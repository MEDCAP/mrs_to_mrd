"""
Convert MR Solutions .MRD raw data to MRD v2.

Three input modes, one per way a scan arrives:

    -t/--tar     Single experiment epsi folder wrapped in a tar file, where each subdirectory folder 
                 represent single repetition. 

    -i/--input   A single .MRD filepath, which is how spectral fid data usually arrives: one file
                 already holds every repetition, so there is nothing to collect or group
    
    -f/--folder  Legacy method for local testing. One experiment folder, walked for its .MRD files

Both -t and -f are given one experiment and take their name from it, so mrs_organize reads the
directories inside only for their scan ids. Which files belong to one scan it decides from the
sequence and the acquisition matrix. All three converge on convert_group_to_mrd, so header choice
and repetition numbering have exactly one implementation.

Nothing here reads a parameter as meaning more than it says: every file converts to acquisitions the
same way, and a file that mrs_organize did not combine with anything is simply a scan of its own.
Whether a scan is calibration data is a question about the scan, answered where that matters rather
than by flagging acquisitions differently on the way in.

Errors are raised only in main(). Past that point a file that cannot be converted is reported on
stderr and skipped, so one unreadable or stray file does not lose the rest of the scan.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path
from typing import BinaryIO, Iterable, List, Optional, Sequence

import numpy as np

# mrd python package
import mrd
import mrs_organize
from MRSreader import MRSdata
from mrs_organize import ScanGroup, sequence_family
from mrs_tar import read_scan_tar

# Time from the start of the excitation pulse to the first sample, in ns. An estimate until the
# sequence reports it: acquisition starts 180us after a 100us pulse begins
TE_NS = np.uint64(1.8E+5)


def keep_mrs(mrs: MRSdata, name: str, family: str) -> bool:
    """
    Whether a parsed file can join a group. Reported rather than raised, so that one bad file does
    not cost the rest of the scan
    Args:
        - mrs: the parsed file
        - name: what to call it on stderr, a path or a tar member name
        - family: the family established for this group, or '' if this is the first file
    Returns:
        - True when the file parsed and belongs to the group's family
    """
    # MRSdata.read_from_file reports a parse failure and returns with rawdata left unset
    if mrs.rawdata is None or mrs.rawdata.size == 0:
        print(f"Skipping {name}: no raw data was read", file=sys.stderr)
        return False
    this_family = sequence_family(mrs.sequence_name)
    # grouping ran off a header-only probe, so this is where the full parse gets to disagree with it
    if family and this_family and this_family != family:
        print(f"Skipping {name}: {this_family} file in a {family} scan", file=sys.stderr)
        return False
    return True


def generate_acquisition(mrs: MRSdata, rep_base: int, rep_count: int) -> Iterable[mrd.StreamItem]:
    """
    Emit one acquisition per (repetition, view) of one MRS file.

    Both families walk the same two loops, bounded by the MRSdata field that names each dimension
    rather than by a position in rawdata.shape. Where the repetitions live is the one thing that
    varies, and counting them off nrepetitions covers every arrangement seen so far:
        - EPSI split one repetition per file:    nrepetitions=1,  nviews=8 ->    8 acquisitions
        - EPSI acquired into a single file:      nrepetitions=N,  nviews=8 ->  N*8 acquisitions
        - spectral, repetitions on the nex axis: nrepetitions=40, nviews=1 ->   40 acquisitions
    rep_base is where this file's repetitions start within the group, so all three number
    identically from the reader's point of view.
    Args:
        - mrs: one parsed MRS file, rawdata indexed
               (nsamples, nviews, nsliceviews, nslices, nechoes, nrepetitions)
        - rep_base: repetition index this file's first repetition maps to
        - rep_count: repetitions in the whole group, for the LAST_IN_REPETITION flag
    Returns:
        - Iterable of mrd.StreamItem.Acquisition
    """
    family = sequence_family(mrs.sequence_name)
    if not family:
        # converted rather than dropped: an unknown sequence still holds acquisitions, and the
        # single view arrangement is the one that assumes least about how they were encoded
        print(f"Unrecognized sequence '{mrs.sequence_name}', converting it as single view data",
              file=sys.stderr)
    TR = np.uint64(mrs.tr * 1.0E+6)                             # mrs.tr in ms, converted to ns
    if family == "epsi":
        # samples hold nswitch echoes, each npoints_per_switch long with a gradient ramp either
        # side. cirrhrat_43_1: 1792 samples / 64 switches = 28, (28 - 12) / 2 = 8 points per ramp.
        # nswitch is a divisor both here and in reconstruction, so a file recording 0 switches is
        # read as 1 rather than being allowed to raise part way through a scan
        nswitch = max(mrs.nswitch, 1)
        points_per_switch = mrs.nsamples // nswitch
        discard = (points_per_switch - mrs.npoints_per_switch) // 2
    for irep in range(mrs.nrepetitions):
        for iview in range(mrs.nviews):
            acq = mrd.Acquisition()
            # MRS acquires on one channel, so add the coil axis: acq.data.shape=(coils=1, samples)
            acq.data = np.expand_dims(mrs.rawdata[:, iview, 0, 0, 0, irep], axis=0)
            # one excitation per (repetition, view): EPSI phase encodes one view per TR and
            # spectral repeats its single view, so counting excitations spaces both correctly
            excitation = np.uint64(irep * mrs.nviews + iview)
            pulse_start = np.uint64(mrs.acquisition_timestamp * 100) + excitation * TR
            acq.head.acquisition_time_stamp_ns = pulse_start + TE_NS
            acq.head.sample_time_ns = mrs.sample_period * 100    # sample_period in units of 100ns
            acq.head.idx.average = mrs.naverages
            acq.head.idx.phase = iview
            repetition = rep_base + irep
            acq.head.idx.repetition = repetition
            acq.head.idx.kspace_encode_step_1 = iview
            acq.head.scan_counter = repetition * mrs.nviews + iview
            if repetition == 0 and iview == 0:
                acq.head.flags |= mrd.AcquisitionFlags.FIRST_IN_REPETITION
            if repetition == rep_count - 1 and iview == mrs.nviews - 1:
                acq.head.flags |= mrd.AcquisitionFlags.LAST_IN_REPETITION
            # both hold on a single view acquisition, so these are two ifs rather than if/elif
            if iview == 0:
                acq.head.flags |= mrd.AcquisitionFlags.FIRST_IN_PHASE
            if iview == mrs.nviews - 1:
                acq.head.flags |= mrd.AcquisitionFlags.LAST_IN_PHASE
            if family == "epsi":
                acq.head.idx.contrast = nswitch                 # echo number in multi-echo
                # recon slices a fixed width off both ends of each echo, so these must stay equal
                acq.head.discard_pre = discard
                acq.head.discard_post = discard
            else:
                acq.head.idx.contrast = 1
                acq.head.discard_pre = 0
                acq.head.discard_post = 0
                acq.phase = np.zeros(mrs.nsamples, dtype=np.float32)   # phase is not recorded
            yield mrd.StreamItem.Acquisition(acq)


def make_header(mrs: MRSdata, meas_id: str, rep_count: int) -> mrd.Header:
    """
    Fill in the MRD header from one file's parameters. Every file in a group was acquired at the
    same matrix, which is what let them be combined, so any of them describes the geometry
    Args:
        - mrs: the file the header describes, chosen by convert_group_to_mrd
        - meas_id: measurement id, e.g. cirrhrat_43_1
        - rep_count: repetitions in the group
    Returns:
        - mrd.Header
    """
    header = mrd.Header()

    subject = mrd.SubjectInformationType()
    subject.patient_id = meas_id            # e.g.) cirrhrat_43_1, KIC_Huh7msps5_08-15-2025.mrs
    header.subject_information = subject

    meas = mrd.MeasurementInformationType()
    meas.sequence_name = mrs.sequence_name
    meas.measurement_id = meas_id
    meas.protocol_name = meas_id.split("_")[0]
    meas.relative_table_position = mrd.ThreeDimensionalFloat(x=mrs.FOVoffset[0] * 1e3,
                                                            y=mrs.FOVoffset[1] * 1e3,
                                                            z=mrs.FOVoffset[2] * 1e3)  # m -> mm
    header.measurement_information = meas

    header.experimental_conditions.h1resonance_frequency_hz = mrs.base_frequency

    encoded_space = mrd.EncodingSpaceType()
    encoded_space.matrix_size = mrd.MatrixSizeType(x=mrs.nsamples, y=mrs.nviews, z=mrs.nslices)
    encoded_space.field_of_view_mm = mrd.FieldOfViewMm(x=mrs.FOV * 1e3, y=mrs.FOV * 1e3, z=0)

    # each limit is the size of one rawdata dimension, as a maximum index
    limits = mrd.EncodingLimitsType()
    limits.kspace_encoding_step_0 = mrd.LimitType(maximum=mrs.nsamples - 1)
    limits.kspace_encoding_step_1 = mrd.LimitType(maximum=mrs.nviews - 1)
    # reconstruction sizes its k-space off the phase limit, so it has to carry the view count too
    limits.phase = mrd.LimitType(maximum=mrs.nviews - 1)
    limits.slice = mrd.LimitType(maximum=mrs.nslices - 1)
    # clamped, so that a file reporting no repetitions still gets a valid header
    limits.repetition = mrd.LimitType(minimum=0, maximum=max(rep_count - 1, 0))

    encoding = mrd.EncodingType()
    encoding.encoded_space = encoded_space
    encoding.encoding_limits = limits
    header.encoding.append(encoding)
    return header


def read_mrs_group(filepaths: Sequence[str], family: str = "") -> List[MRSdata]:
    """
    Parse a group of .MRD files from disk, dropping the ones that cannot contribute
    Args:
        - filepaths: paths in acquisition order
        - family: the family mrs_organize established for the group, so that the full parse is
          checked against the header-only probe grouping ran on
    Returns:
        - parsed MRSdata in the same order, possibly shorter than filepaths
    """
    mrs_list: List[MRSdata] = []
    for filepath in filepaths:
        mrs = MRSdata()                                 # one instance per file, they are all kept
        mrs.read_from_file(filepath)
        if not keep_mrs(mrs, str(filepath), family):
            continue
        family = family or sequence_family(mrs.sequence_name)
        mrs_list.append(mrs)
    return mrs_list


def convert_group_to_mrd(mrs_list: Sequence[MRSdata], meas_id: str, output: BinaryIO) -> bool:
    """
    Write one scan's worth of parsed MRS files as a single MRD v2 stream.

    Repetitions are counted off nrepetitions and summed across the group, which covers EPSI split
    one repetition per file, EPSI acquired with all its repetitions in one file, and spectral
    holding them on the nex axis, without the caller having to say which of the three it has. A
    group of one file is the same thing with one term in the sum.
    Args:
        - mrs_list: parsed files in acquisition order
        - meas_id: measurement id recorded in the header
        - output: writable binary stream. Must be a file object, not a path: BinaryMrdWriter
                  special-cases str only, so a Path would be mistaken for a stream
    Returns:
        - True when a stream was written
    """
    if not mrs_list:
        print(f"Nothing to convert for {meas_id}", file=sys.stderr)
        return False
    rep_count = sum(mrs.nrepetitions for mrs in mrs_list)
    print(f"Converting {meas_id}: {len(mrs_list)} files, {rep_count} repetitions", file=sys.stderr)
    # the writer must be closed to emit the end-of-stream sentinel, hence the with block
    with mrd.BinaryMrdWriter(output) as writer:
        writer.write_header(make_header(mrs_list[0], meas_id, rep_count))
        rep_base = 0
        for mrs in mrs_list:
            writer.write_data(generate_acquisition(mrs, rep_base, rep_count))
            rep_base += mrs.nrepetitions
    return True


def report_warnings(groups: Sequence[ScanGroup]) -> None:
    """
    Print what mrs_organize noticed about the shape of the input: a file from another session, a
    scan id that never arrived, two series in one experiment. These are warnings rather than errors,
    because a scan with a file missing still converts, and saying so beside the converted stream is
    more use than refusing to write one
    Args:
        - groups: every group in this run
    """
    for warning in mrs_organize.check_groups(groups):
        print(f"WARNING {warning}", file=sys.stderr)
    for group in groups:
        for warning in mrs_organize.check_group(group):
            print(f"WARNING {warning}", file=sys.stderr)


def convert_folder_to_mrd(folder: Path, dry_run: bool = False) -> bool:
    """
    Walk one experiment folder for .MRD files and convert each scan to its own stream inside it
        spectral: KIC_Huh7msps5_08-15-2025.mrs/KIC_huh7_5.MRD
            -> KIC_Huh7msps5_08-15-2025.mrs/KIC_Huh7msps5_08-15-2025.mrs_1puls_extrf_KIC.mrd2
        EPSI:     cirrhrat_43_1/epsi/24804/24804_000_0.MRD, one repetition per scan directory
            -> cirrhrat_43_1/cirrhrat_43_1_epsigre_combined.mrd2
        EPSI:     cirrhrat_43_1/epsi/24792/24792_000_0.MRD, acquired at another matrix
            -> cirrhrat_43_1/cirrhrat_43_1_24792.mrd2
    Args:
        - folder: the experiment folder to walk
        - dry_run: report the grouping and what looks wrong with it, converting nothing
    Returns:
        - True when at least one scan was written, or when a dry run found something to convert
    """
    groups = mrs_organize.organize_folder(folder)
    if dry_run:
        mrs_organize.report(groups)
        return bool(groups)
    report_warnings(groups)
    written = False
    for group in groups:
        mrs_list = read_mrs_group(group.files, group.family)
        if not mrs_list:
            print(f"No usable files in {group.meas_id}, skipping", file=sys.stderr)
            continue
        print(f"Writing {group.output_path}", file=sys.stderr)
        with open(group.output_path, "wb") as output:
            written |= convert_group_to_mrd(mrs_list, group.meas_id, output)
    return written


def convert_file_to_mrd(input_path: Path, output_path: Optional[Path] = None) -> bool:
    """
    Convert a single .MRD file. Spectral fid data arrives this way: one file already holds every
    repetition on its nex axis, so there is nothing to collect or group
    Args:
        - input_path: the .MRD file
        - output_path: where to write, defaulting to the name its scan would have been given
    Returns:
        - True when a stream was written
    """
    groups = mrs_organize.organize_folder(input_path)
    if not groups:
        print(f"Nothing to convert in {input_path}", file=sys.stderr)
        return False
    group = groups[0]
    report_warnings(groups)
    destination = output_path or Path(group.output_path)
    mrs_list = read_mrs_group(group.files, group.family)
    if not mrs_list:
        print(f"Nothing to convert in {input_path}", file=sys.stderr)
        return False
    print(f"Converting {input_path} to {destination}", file=sys.stderr)
    with open(destination, "wb") as output:
        return convert_group_to_mrd(mrs_list, group.meas_id, output)


def convert_tar_to_mrd(tar_path: Path, output_path: Path, meas_id_override: str = "") -> bool:
    """
    Convert tar filed directory of single experiment epsi MRS data folder

    Tyger hands a job its input buffer as a named FIFO, which is strictly sequential, so the
    archive is read forward once and held in memory: the parameter block sits after the raw data at
    EOF, the .SPR sidecar can follow the .MRD members in tar order, and the header needs the whole
    group before the first acquisition can be written.

    A Tyger job has one output buffer, so this writes one stream. The members are grouped exactly
    as a folder is, which turns "the caller tarred one scan" from an assumption into a checked
    precondition: more than one scan in the archive is an error, unless --output names a directory
    to write them all into
    Args:
        - tar_path: tar archive of one scan directory, or a FIFO carrying one
        - output_path: where to write the stream, which may also be a FIFO, or a directory when the
          archive holds more than one scan
        - meas_id_override: measurement id to record instead of the tar's root directory name
    Returns:
        - True when a stream was written
    """
    # opening a FIFO for read blocks until the buffer sidecar opens the write end
    with open(tar_path, "rb") as tar_stream:
        tar_meas_id, spr_frequency, members = read_scan_tar(tar_stream)
    fallback = meas_id_override or tar_meas_id or Path(output_path).stem
    if not (meas_id_override or tar_meas_id):
        # nothing to name the scan after: the archive holds no single root directory
        print(f"No single root directory in the tar, calling this scan {fallback}", file=sys.stderr)
    if not spr_frequency:
        print("No base frequency from a .SPR sidecar in the tar", file=sys.stderr)

    payloads = dict(members)
    groups = mrs_organize.organize_members(members, fallback_meas_id=fallback)
    report_warnings(groups)
    into_directory = output_path.is_dir()
    if len(groups) > 1 and not into_directory:
        raise ValueError(f"the tar holds {len(groups)} scans "
                         f"({', '.join(group.output_name for group in groups)}) but --output names "
                         f"a single stream: tar one scan directory, or point --output at a "
                         f"directory")

    written = False
    for group in groups:
        mrs_list: List[MRSdata] = []
        for name in group.files:
            mrs = MRSdata()
            mrs.parse_from_buffer(payloads[name])
            # the .MRD may defer its frequency to the sidecar, as read_from_file does on disk
            mrs.set_base_frequency(spr_frequency)
            if not keep_mrs(mrs, name, group.family):
                continue
            mrs_list.append(mrs)
        # meas_id_override names the scan itself, so it wins over what grouping inferred
        meas_id = meas_id_override or group.meas_id or fallback
        destination = output_path / group.output_name if into_directory else output_path
        # opened even with nothing to write: on Tyger this is the output buffer's FIFO, and closing
        # it empty is what tells the sidecar the job produced nothing, rather than leaving it
        # blocked on a stream that never opens
        with open(destination, "wb") as output:
            written |= convert_group_to_mrd(mrs_list, meas_id, output)
    if not groups and not into_directory:
        with open(output_path, "wb"):       # nothing grouped, but the buffer still has to close
            pass
    return written


def main() -> int:
    """
    Convert MRS data to MRD2, in exactly one of three input modes.

    Every check that can stop the run lives here. Past this point a file that cannot be converted is
    reported and skipped, so the run ends with a status rather than a traceback
    Returns:
        - 0 when at least one stream was written, 1 when nothing was
    """
    parser = argparse.ArgumentParser(description="Convert MR Solutions MRS data to MRD2 format")
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("-t", "--tar", type=Path,
                      help="tar archive of one scan directory, or a FIFO carrying one")
    mode.add_argument("-i", "--input", type=Path,
                      help="single MRS .MRD file")
    mode.add_argument("-f", "--folder", type=Path,
                      help="directory to walk for MRS .MRD files")
    parser.add_argument("-o", "--output", type=Path,
                        help="file or FIFO to write the MRD2 stream to. Required with --tar "
                             "(default: $OUTPUT_PIPE), optional with --input, unused with --folder")
    parser.add_argument("--meas-id", default="",
                        help="measurement id to record instead of the folder's root directory name")
    parser.add_argument("-n", "--dry-run", action="store_true",
                        help="with --folder, report how the files group and what looks wrong with "
                             "them, without converting anything")
    args = parser.parse_args()

    if args.tar:
        output = args.output or os.environ.get("OUTPUT_PIPE")
        if not output:
            parser.error("--output is required with --tar when $OUTPUT_PIPE is unset")
        print(f"Converting tarred scan {args.tar} to {output}", file=sys.stderr)
        written = convert_tar_to_mrd(args.tar, Path(output), args.meas_id)
    elif args.input:
        if not args.input.is_file():
            parser.error(f"{args.input} is not a file")
        written = convert_file_to_mrd(args.input, args.output)
    else:
        if not args.folder.is_dir():
            parser.error(f"{args.folder} is not a directory")
        if args.output:
            parser.error("--output does not apply to --folder: each scan is written beside its "
                         "own files")
        print(f"Converting folder {args.folder}", file=sys.stderr)
        written = convert_folder_to_mrd(args.folder, args.dry_run)

    if not written:
        print("Nothing was converted", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
