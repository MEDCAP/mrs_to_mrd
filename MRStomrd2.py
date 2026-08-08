"""
Convert MR Solutions .MRD raw data to MRD v2.

Three input modes, one per way a scan arrives:

    -t/--tar     one scan directory wrapped in a tar, possibly arriving on a FIFO. This is how
                 Tyger delivers a job's input buffer, and a FIFO is strictly sequential, so nothing
                 on this path seeks
    -i/--input   a single .MRD file, which is how spectral fid data usually arrives: one file
                 already holds every repetition, so there is nothing to collect or group
    -f/--folder  a directory walked for .MRD files, grouped into scans by --unifylevel. This is the
                 local mode the *_recon.sh scripts drive

All three converge on convert_group_to_mrd, so phantom handling, header choice and repetition
numbering have exactly one implementation.

Errors are raised only in main(). Past that point a file that cannot be converted is reported on
stderr and skipped, so one unreadable or stray file does not lose the rest of the scan. A scan of
nothing but phantom data still converts.
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
from MRSreader import MRSdata
from mrs_tar import read_scan_tar

# Sequence name fragments identifying the two acquisition families. EVO1 and EVO2 spell the
# spectral sequence differently, hence two fragments for it: '1puls_extrf_KIC', 'fid...'
EPSI_SEQUENCES = ("epsi",)
SPECTRAL_SEQUENCES = ("1pul", "fid")

# Time from the start of the excitation pulse to the first sample, in ns. An estimate until the
# sequence reports it: acquisition starts 180us after a 100us pulse begins
TE_NS = np.uint64(1.8E+5)

# .MRD files are written one scan per numbered folder, so the name of the file itself carries no
# scan identity. Output is named after the scan directory instead, matching the
# <scan directory>_recon.mrd2 that reconstruction writes beside it
RAW_SUFFIX = "_raw.mrd2"


def sequence_family(mrs: MRSdata) -> str:
    """
    Which acquisition family a file belongs to, read from the sequence it was acquired with rather
    than from the folder depth it was found at, so that all three input modes classify a file the
    same way. --unifylevel is then only about how to group files, not about what is in them
    Args:
        - mrs: a parsed MRS file
    Returns:
        - 'epsi', 'spectral', or '' when the sequence name matches neither family
    """
    name = mrs.sequence_name.lower()
    if any(fragment in name for fragment in EPSI_SEQUENCES):
        return "epsi"
    if any(fragment in name for fragment in SPECTRAL_SEQUENCES):
        return "spectral"
    return ""


def is_phantom(mrs: MRSdata) -> bool:
    """
    Whether a file holds phantom rather than acquired data. A scan averaged more than once is a
    separate calibration scan: it is converted like everything else, but it is flagged
    IS_NAVIGATION_DATA and it consumes no repetition index
    Args:
        - mrs: a parsed MRS file
    Returns:
        - True for phantom data
    """
    return mrs.naverages > 1


def keep_mrs(mrs: MRSdata, name: str, family: str) -> bool:
    """
    Whether a parsed file can join a group. Reported rather than raised, so that one bad file does
    not cost the rest of the scan
    Args:
        - mrs: the parsed file
        - name: what to call it on stderr, a path or a tar member name
        - family: the family already established for this group, or '' if this is the first file
    Returns:
        - True when the file parsed and belongs to the group's family
    """
    # MRSdata.read_from_file reports a parse failure and returns with rawdata left unset
    if mrs.rawdata is None or mrs.rawdata.size == 0:
        print(f"Skipping {name}: no raw data was read", file=sys.stderr)
        return False
    this_family = sequence_family(mrs)
    if not this_family:
        print(f"Skipping {name}: unrecognized sequence '{mrs.sequence_name}'", file=sys.stderr)
        return False
    # an EPSI scan folder can hold a stray fid, and a spectral folder a stray epsi
    if family and this_family != family:
        print(f"Skipping {name}: {this_family} file in a {family} scan", file=sys.stderr)
        return False
    if is_phantom(mrs):
        print(f"Phantom data at {name} ({mrs.naverages} averages)", file=sys.stderr)
    return True


def generate_acquisition(mrs: MRSdata, rep_base: int, rep_count: int) -> Iterable[mrd.StreamItem]:
    """
    Emit one acquisition per (repetition, view) of one MRS file.

    Both families walk the same two loops, bounded by the MRSdata field that names each dimension
    rather than by a position in rawdata.shape. Where the repetitions live is the one thing that
    varies, and counting them off nrepetitions covers every arrangement seen so far:
        - EPSI split one repetition per file:  nrepetitions=1,  nviews=8  ->    8 acquisitions
        - EPSI acquired into a single file:    nrepetitions=N,  nviews=8  ->  N*8 acquisitions
        - spectral, repetitions on the nex axis: nrepetitions=40, nviews=1 ->  40 acquisitions
    rep_base is where this file's repetitions start within the group, so all three number
    identically from the reader's point of view.
    Args:
        - mrs: one parsed MRS file, rawdata indexed
               (nsamples, nviews, nsliceviews, nslices, nechoes, nrepetitions)
        - rep_base: repetition index this file's first repetition maps to
        - rep_count: non-phantom repetitions in the whole group, for the LAST_IN_REPETITION flag
    Returns:
        - Iterable of mrd.StreamItem.Acquisition, empty when the sequence belongs to neither family
    """
    family = sequence_family(mrs)
    if not family:
        print(f"No acquisitions from unrecognized sequence '{mrs.sequence_name}'", file=sys.stderr)
        return
    phantom = is_phantom(mrs)
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
            acq.head.sample_time_ns = mrs.sampleperiod * 100    # sampleperiod in units of 100ns
            acq.head.idx.average = mrs.naverages
            acq.head.idx.phase = iview
            if phantom:
                acq.head.flags |= mrd.AcquisitionFlags.IS_NAVIGATION_DATA
                # deliberately no repetition or kspace_encode_step_1: a phantom is calibration, and
                # its matrix is not the acquired one (2176x12 against 1792x8 in cirrhrat_43_1), so
                # indexing its views would run off the end of a k-space sized from the header
            else:
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


def make_header(mrs: MRSdata,
                meas_id: str,
                rep_count: int,
                phantom_idx_list: Sequence[int]) -> mrd.Header:
    """
    Fill in the MRD header from one file's parameters. Everything describing the acquisition
    geometry comes from acquired rather than phantom data, since a phantom is a separate
    calibration scan whose matrix size is not the one it calibrates
    Args:
        - mrs: the file the header describes, chosen by convert_group_to_mrd
        - meas_id: measurement id, e.g. cirrhrat_43_1
        - rep_count: non-phantom repetitions in the group
        - phantom_idx_list: positions in the group that held phantom data, recorded as provenance
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
    # clamped, so that a group holding nothing but phantoms still gets a valid header
    limits.repetition = mrd.LimitType(minimum=0, maximum=max(rep_count - 1, 0))

    encoding = mrd.EncodingType()
    encoding.encoded_space = encoded_space
    encoding.encoding_limits = limits
    header.encoding.append(encoding)

    # provenance: which files in the group were phantoms. The acquisitions carry the same fact as
    # their IS_NAVIGATION_DATA flag, which is what reconstruction reads
    user_params = mrd.UserParametersType()
    for idx in phantom_idx_list:
        user_params.user_parameter_long.append(mrd.UserParameterLongType(name="phantom_idx",
                                                                        value=idx))
    header.user_parameters = user_params
    return header


def collect_mrd_files(folder: Path) -> List[Path]:
    """
    Find all MRS .MRD files in a folder and its subdirectories recursively
    Args:
        - folder: directory to walk
    Returns:
        - List of Path, sorted, because iterdir order is arbitrary and a file's position in this
          list becomes its acquisition repetition index
    """
    mrd_filepath_list: List[Path] = []
    for entry in sorted(folder.iterdir()):
        if entry.is_dir():
            mrd_filepath_list.extend(collect_mrd_files(entry))
        elif entry.suffix == ".MRD":
            mrd_filepath_list.append(entry)
    return mrd_filepath_list


def group_mrd_files(folder: Path, unifylevel: int) -> List[List[Path]]:
    """
    Group the .MRD files under a folder into scans, by the directory they share unifylevel levels up
        spectral, unifylevel=1: KIC_data/KIC_Huh7msps5_08-15-2025.mrs/KIC_huh7_5.MRD
            grouped on the file's own parent, so one group per .mrs folder
        EPSI, unifylevel=3: cirrhrat_data/cirrhrat_43_1/epsi/24804/24804_000_0.MRD
            grouped on parts[:-2]=(cirrhrat_data, cirrhrat_43_1, epsi), so one group per scan
    Args:
        - folder: directory to walk
        - unifylevel: path components to strip from a file to reach its scan directory
    Returns:
        - List of lists of Path
    """
    mrd_filepath_list = collect_mrd_files(folder)
    mrd_file_groups: List[List[Path]] = []
    # parts[:-0] is the empty tuple rather than the whole path, so unifylevel=1 has to be clamped
    # or every file compares equal and lands in one group
    trim = max(unifylevel - 1, 1)
    for filepath in mrd_filepath_list:
        for group in mrd_file_groups:
            if filepath.parts[:-trim] == group[0].parts[:-trim]:
                group.append(filepath)
                break
        else:
            mrd_file_groups.append([filepath])
    print(f"Grouped {len(mrd_filepath_list)} files into {len(mrd_file_groups)} scans",
          file=sys.stderr)
    return mrd_file_groups


def read_mrs_group(filepaths: Sequence[Path]) -> List[MRSdata]:
    """
    Parse a group of .MRD files from disk, dropping the ones that cannot contribute
    Args:
        - filepaths: paths in acquisition order
    Returns:
        - parsed MRSdata in the same order, possibly shorter than filepaths
    """
    mrs_list: List[MRSdata] = []
    family = ""
    for filepath in filepaths:
        mrs = MRSdata()                                 # one instance per file, they are all kept
        mrs.read_from_file(filepath)
        if not keep_mrs(mrs, str(filepath), family):
            continue
        family = family or sequence_family(mrs)
        mrs_list.append(mrs)
    return mrs_list


def convert_group_to_mrd(mrs_list: Sequence[MRSdata], meas_id: str, output: BinaryIO) -> bool:
    """
    Write one scan's worth of parsed MRS files as a single MRD v2 stream.

    Repetitions are counted off nrepetitions and summed across the group, which covers EPSI split
    one repetition per file, EPSI acquired with all its repetitions in one file, and spectral
    holding them on the nex axis, without the caller having to say which of the three it has.
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
    phantom_idx_list = [i for i, mrs in enumerate(mrs_list) if is_phantom(mrs)]
    # phantoms hold no repetition index, so they do not count towards the total
    rep_count = sum(mrs.nrepetitions for mrs in mrs_list if not is_phantom(mrs))
    # a phantom's matrix size is not the acquired one, so the header has to come from acquired
    # data. Phantoms sort anywhere in a group, the first two positions included, hence the search
    header_mrs = next((mrs for mrs in mrs_list if not is_phantom(mrs)), mrs_list[0])
    if is_phantom(header_mrs):
        print(f"{meas_id} is phantom data only, taking the header from the first file",
              file=sys.stderr)
    print(f"Converting {meas_id}: {len(mrs_list)} files, {len(phantom_idx_list)} phantom, "
          f"{rep_count} repetitions", file=sys.stderr)
    # the writer must be closed to emit the end-of-stream sentinel, hence the with block
    with mrd.BinaryMrdWriter(output) as writer:
        writer.write_header(make_header(header_mrs, meas_id, rep_count, phantom_idx_list))
        rep_base = 0
        for mrs in mrs_list:
            writer.write_data(generate_acquisition(mrs, rep_base, rep_count))
            if not is_phantom(mrs):
                rep_base += mrs.nrepetitions
    return True


def output_filepath(scan_folder: Path, meas_id: str) -> Path:
    """
    Where a converted scan is written, named after the scan directory so that it stays paired with
    the <scan directory>_recon.mrd2 reconstruction writes beside it
    Args:
        - scan_folder: directory the scan is named after
        - meas_id: measurement id
    Returns:
        - Path to write the MRD2 stream to
    """
    return scan_folder / f"{meas_id}{RAW_SUFFIX}"


def convert_tar_to_mrd(tar_path: Path, output_path: Path, meas_id_override: str = "") -> bool:
    """
    Convert one scan directory delivered as a tar.

    Tyger hands a job its input buffer as a named FIFO, which is strictly sequential, so the
    archive is read forward once and held in memory: the parameter block sits after the raw data at
    EOF, the .SPR sidecar can follow the .MRD members in tar order, and the header needs the whole
    group before the first acquisition can be written.
    Args:
        - tar_path: tar archive of one scan directory, or a FIFO carrying one
        - output_path: where to write the stream, which may also be a FIFO
        - meas_id_override: measurement id to record instead of the tar's root directory name
    Returns:
        - True when a stream was written
    """
    # opening a FIFO for read blocks until the buffer sidecar opens the write end
    with open(tar_path, "rb") as tar_stream:
        meas_id, spr_frequency, members = read_scan_tar(tar_stream)
    meas_id = meas_id_override or meas_id
    if not meas_id:
        # nothing to name the scan after: the archive holds no single root directory
        meas_id = Path(output_path).name.replace(RAW_SUFFIX, "")
        print(f"No single root directory in the tar, calling this scan {meas_id}", file=sys.stderr)
    if not spr_frequency:
        print("No base frequency from a .SPR sidecar in the tar", file=sys.stderr)

    mrs_list: List[MRSdata] = []
    family = ""
    for name, filebytes in members:
        mrs = MRSdata()
        mrs.parse_from_buffer(filebytes)
        # the .MRD may defer its frequency to the sidecar, as MRSdata.read_from_file does on disk
        if mrs.base_frequency_in_SPR:
            mrs.set_base_frequency(spr_frequency)
        if not keep_mrs(mrs, name, family):
            continue
        family = family or sequence_family(mrs)
        mrs_list.append(mrs)

    # opened even with nothing to write: on Tyger this is the output buffer's FIFO, and closing it
    # empty is what tells the sidecar the job produced nothing, rather than leaving it blocked
    with open(output_path, "wb") as output:
        return convert_group_to_mrd(mrs_list, meas_id, output)


def convert_file_to_mrd(input_path: Path, output_path: Optional[Path] = None) -> bool:
    """
    Convert a single .MRD file. Spectral fid data arrives this way: one file already holds every
    repetition on its nex axis, so there is nothing to collect or group
    Args:
        - input_path: the .MRD file
        - output_path: where to write, defaulting to <scan directory>_raw.mrd2 beside the input
    Returns:
        - True when a stream was written
    """
    meas_id = input_path.parent.name
    if output_path is None:
        output_path = output_filepath(input_path.parent, meas_id)
    mrs_list = read_mrs_group([input_path])
    if not mrs_list:
        print(f"Nothing to convert in {input_path}", file=sys.stderr)
        return False
    print(f"Converting {input_path} to {output_path}", file=sys.stderr)
    with open(output_path, "wb") as output:
        return convert_group_to_mrd(mrs_list, meas_id, output)


def convert_folder_to_mrd(folder: Path, unifylevel: int) -> bool:
    """
    Walk a folder for .MRD files and convert each scan to its own stream beside its files
        spectral, unifylevel=1: KIC_data/KIC_Huh7msps5_08-15-2025.mrs/KIC_huh7_5.MRD
            -> KIC_data/KIC_Huh7msps5_08-15-2025.mrs/KIC_Huh7msps5_08-15-2025.mrs_raw.mrd2
        EPSI, unifylevel=3: cirrhrat_data/cirrhrat_43_1/epsi/24804/24804_000_0.MRD
            -> cirrhrat_data/cirrhrat_43_1/cirrhrat_43_1_raw.mrd2
    Args:
        - folder: directory to walk
        - unifylevel: path components to strip from a file to reach its scan directory
    Returns:
        - True when at least one scan was written
    """
    written = False
    for group in group_mrd_files(folder, unifylevel):
        meas_id = group[0].parts[-(unifylevel + 1)]             # e.g.) cirrhrat_43_1
        scan_folder = Path(*group[0].parts[:-unifylevel])       # e.g.) cirrhrat_data/cirrhrat_43_1
        mrs_list = read_mrs_group(group)
        if not mrs_list:
            print(f"No usable files in {meas_id}, skipping", file=sys.stderr)
            continue
        with open(output_filepath(scan_folder, meas_id), "wb") as output:
            written |= convert_group_to_mrd(mrs_list, meas_id, output)
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
    parser.add_argument("-u", "--unifylevel", type=int, default=1,
                        help="with --folder, path components to strip from a file to reach its "
                             "scan directory: 1 for spectral, 3 for EPSI (default: 1)")
    parser.add_argument("-o", "--output", type=Path,
                        help="file or FIFO to write the MRD2 stream to. Required with --tar "
                             "(default: $OUTPUT_PIPE), optional with --input, unused with --folder")
    parser.add_argument("--meas-id", default="",
                        help="with --tar, measurement id to record instead of the archive's root "
                             "directory name")
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
        print(f"Converting folder {args.folder} with unify level {args.unifylevel}",
              file=sys.stderr)
        written = convert_folder_to_mrd(args.folder, args.unifylevel)

    if not written:
        print("Nothing was converted", file=sys.stderr)
        return 1
    return 0


# -f with -u 3 consolidates the files as appropriate for EPSI, spectral data uses -u 1
if __name__ == "__main__":
    raise SystemExit(main())
