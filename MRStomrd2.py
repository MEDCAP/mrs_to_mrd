"""
Convert MR Solutions .MRD raw data to MRD v2.

Three input modes, one per way a scan arrives:

    -t/--tar     Single experiment epsi folder wrapped in a tar file, where each subdirectory folder 
                 represent single repetition and pre-scan data that do not have the same dimension 
                 as rawdata. 

    -i/--input   A single .MRD filepath, or a directory holding exactly one, which is how spectral
                 fid data usually arrives: one file already holds every repetition, so there is
                 nothing to collect or group

    -f/--folder  Legacy method for local testing. One experiment folder, walked for its .MRD files

    -w/--window  Plot the echo position and sampling window

Both -t and -f represent one experiment and take their name from it, so mrs_organize reads the
directories inside only for their scan ids. Which files belong to one scan it decides from the
sequence and the acquisition matrix, and a file's position in the group decides where its
repetitions start, so a series split one file per repetition and a series already carried on one
file's repetition axis number identically in the stream.
"""

from __future__ import annotations

import argparse
import os
import sys
from itertools import product
from pathlib import Path
from typing import Iterable, Optional

import numpy as np

# mrd python package
import mrd
from MRSreader import MRSdata
import mrs_organize
from mrs_organize import read_scan_tar
from epsi_window import check_peak_position


def generate_acquisition(mrs: MRSdata,
                         rep_count: int = 0) -> Iterable[mrd.StreamItem]:
    """
    Emit one acquisition per adc event from one MRS file.

    Only this file's own repetitions are walked, and rep_idx is where they start in the group, so
    the two ways a series arrives number identically in the stream
        - EPSI split one repetition per file: nrepetitions=1, rep_idx counting up per file
        - EPSI or spectral in a single file:  nrepetitions=N, rep_idx=0

    Each axis takes the MRD index that means the same thing
        nviews       -> kspace_encode_step_1
        nsliceviews  -> kspace_encode_step_2
        nslices      -> slice
        nechoes      -> contrast
        nrepetitions -> repetition
    Args:
        - mrs: one parsed MRS file, rawdata indexed
               (nsamples, nviews, nsliceviews, nslices, nechoes, nrepetitions)
        - rep_count: repetitions in the whole group, for the LAST_IN_REPETITION flag. Defaults to
          this file being the whole group
    Returns:
        - Iterable of mrd.StreamItem.Acquisition
    """
    # pre_scan data does not encode rep_count
    if rep_count is None:
        rep_count = mrs.nrepetitions
    grid = product(range(rep_count), range(mrs.nechoes), range(mrs.nslices),
                   range(mrs.nsliceviews), range(mrs.nviews))
    for counter, (irep, iecho, islice, isliceview, iview) in enumerate(grid):
        acq = mrd.Acquisition()
        acq.head.acquisition_time_stamp_ns = np.uint64(mrs.acquisition_timestamp * 100) # mrs timestamp is in 100ns
        acq.head.sample_time_ns = mrs.sample_period * 100       # sample_period in units of 100ns
        acq.head.idx.average = mrs.naverages                    # number of averages, not index
        acq.head.idx.repetition = irep
        if irep == 0:
            acq.head.flags |= mrd.AcquisitionFlags.FIRST_IN_REPETITION
        if irep == rep_count - 1:
            acq.head.flags |= mrd.AcquisitionFlags.LAST_IN_REPETITION
        acq.head.idx.kspace_encode_step_1 = iview
        acq.head.idx.kspace_encode_step_2 = isliceview
        acq.head.idx.slice = islice
        acq.head.idx.contrast = iecho                           # index of echoes
        # unique and increasing across the whole group, since rep_idx places this file's block
        acq.head.scan_counter = counter
        if iview == 0:
            acq.head.flags |= mrd.AcquisitionFlags.FIRST_IN_PHASE
        if iview == mrs.nviews - 1:
            acq.head.flags |= mrd.AcquisitionFlags.LAST_IN_PHASE
        # epsi sequence acquire in the order of rampup -> readout -> rampdown -> rephasing, so one
        # switch period is ramp + npoints_per_switch + 3 ramps. tramp is in us and sample_period in
        # units of 100ns, hence the /10 that puts them both in samples
        # example cirrhrat_43_1: tramp=112us / 28 = 4, and 4 + 12 + 12 = 1792 samples / 64 switches
        if mrs.nswitches > 1:
            ramp = mrs.tramp // (mrs.sample_period // 10)
            acq.head.user_int = [mrs.nswitches, mrs.npoints_per_switch]
            acq.head.discard_pre = ramp
            acq.head.discard_post = ramp * 3    # rampdown + rephasing, which is two ramps
        if mrs.naverages > 1 and mrs.nrepetitions == 1 and 'epsi' in mrs.sequence_name:
            acq.head.flags |= mrd.AcquisitionFlags.IS_NAVIGATION_DATA
        # MRS acquires on one channel, so add the coil axis: acq.data.shape=(coils=1, samples)
        if mrs.nrepetitions == 1:
            irep = 0
        acq.data = np.expand_dims(mrs.rawdata[:, iview, isliceview, islice, iecho, irep], axis=0)
        yield mrd.StreamItem.Acquisition(acq)


def generate_header(mrs: MRSdata, meas_id: str, rep_count: Optional[int]=None) -> mrd.Header:
    """
    Fill in the MRD header from one file's parameters. Every file in a group was acquired at the
    same matrix, which is what let them be combined, so any of them describes the geometry
    Args:
        - mrs: the file the header describes, chosen by convert_group_to_mrd
        - meas_id: measurement id, e.g. cirrhrat_43_1
        - rep_count: total count of repetitions in the group
    Returns:
        - mrd.Header
    """
    header = mrd.Header()

    subject = mrd.SubjectInformationType()
    subject.patient_id = meas_id            # e.g.) cirrhrat_43_1, KIC_Huh7msps5_08-15-2025.mrs
    header.subject_information = subject

    seq = mrd.SequenceParametersType()
    seq.t_r = [mrs.tr]
    seq.t_e = [mrs.te]
    seq.flip_angle_deg = [mrs.flip_angle]
    header.sequence_parameters = seq 

    meas = mrd.MeasurementInformationType()
    meas.sequence_name = mrs.sequence_name
    meas.measurement_id = meas_id
    meas.protocol_name = meas_id.split("_")[0]
    meas.relative_table_position = mrd.ThreeDimensionalFloat(x=mrs.FOVoffset[0] * 1e3,
                                                             y=mrs.FOVoffset[1] * 1e3,
                                                             z=mrs.FOVoffset[2] * 1e3)  # m -> mm
    header.measurement_information = meas
    
    header.experimental_conditions.h1resonance_frequency_hz = mrs.base_frequency
    
    user_param = mrd.UserParametersType()
    user_param.user_parameter_long.append(
            mrd.UserParameterLongType(name="tramp", value=int(mrs.tramp)))
    header.user_parameters = user_param

    encoded_space = mrd.EncodingSpaceType()
    encoded_space.matrix_size = mrd.MatrixSizeType(x=mrs.nsamples, y=mrs.nviews, z=mrs.nsliceviews)
    encoded_space.field_of_view_mm = mrd.FieldOfViewMm(x=mrs.FOV * 1e3, y=mrs.FOV * 1e3, z=0)

    # each limit is the size of one rawdata dimension, as a maximum index, and every dimension
    limits = mrd.EncodingLimitsType()
    limits.kspace_encoding_step_0 = mrd.LimitType(maximum=mrs.nsamples - 1)
    limits.kspace_encoding_step_1 = mrd.LimitType(maximum=mrs.nviews - 1)
    limits.kspace_encoding_step_2 = mrd.LimitType(maximum=mrs.nsliceviews - 1)
    # reconstruction sizes its k-space off the phase limit, so it has to carry the view count too
    limits.phase = mrd.LimitType(maximum=mrs.nviews - 1)
    limits.slice = mrd.LimitType(maximum=mrs.nslices - 1)
    limits.contrast = mrd.LimitType(maximum=mrs.nechoes - 1)     
    if rep_count is None:
        rep_count = mrs.nrepetitions
    limits.repetition = mrd.LimitType(minimum=0, maximum=rep_count - 1)
    encoding = mrd.EncodingType()
    encoding.encoded_space = encoded_space
    encoding.encoding_limits = limits
    header.encoding.append(encoding)
    return header

def convert_folder_to_mrd(folder: Path,
                          dry_run: bool = False,
                          check_window: bool = False) -> bool:

    """
    Walk one experiment folder for .MRD files and convert each scan to its own stream inside it
        spectral: {experiment_folder}/{KIC_Huh7msps5_08-15-2025.MRD}
            -> KIC_Huh7msps5_08-15-2025.mrs/KIC_Huh7msps5_08-15-2025.mrs_1puls_extrf_KIC.mrd2
        EPSI:     {experiment_folder}/{modal}/{scan_id}/{24804_000_0.MRD}, one repetition per scan directory
            -> cirrhrat_43_1/cirrhrat_43_1_epsigre_combined.mrd2
        EPSI:     a subdirectory among those has navg>1 and nrep=1 instead - an averaged prescan,
                  converted into the same stream as the data it calibrates rather than a file of
                  its own, its acquisitions flagged IS_NAVIGATION_DATA
    Args:
        - folder: the experiment folder to walk
        - dry_run: report the grouping and what looks wrong with it, converting nothing
        - check_window: report sampling window and try shifting the echo position
    Returns:
        - True when at least one scan was written, or when a dry run found something to convert
    """
    # scan subfolders of provided folder and group them into rawdata and prescan data as object ScanGroup
    grouped_files = mrs_organize.organize_folder(folder)
    if dry_run:
        mrs_organize.report(grouped_files)
        return bool(grouped_files)
    if grouped_files is None:
        print(f"No data to convert in {folder}", file=sys.stderr)
        return False
    # plot signals for each switch
    if check_window:
        check_peak_position(grouped_files)

    # If mrs.nrepetitions>1, single file carries multiple repetitions and rep_count=mrs.nrepetitions
    # If mrs.nrepetitions==1, a group of files represent entire repetitions and rep_count={number of rawdata files} 
    rep_count = grouped_files.nrepetitions
    writer: Optional[mrd.BinaryMrdWriter] = None
    try:
        # first convert raw data files before phantom files in the group
        for filepath in grouped_files.rawdata_file_list:
            mrs = MRSdata()                     # one at a time, released once written
            mrs.read_from_file(filepath)
            if writer is None:
                print(f'Writing file at {grouped_files.output_path}', file=sys.stderr)
                writer = mrd.BinaryMrdWriter(grouped_files.output_path)
                writer.write_header(generate_header(mrs, grouped_files.meas_id, rep_count))
            writer.write_data(generate_acquisition(mrs, rep_count))
        # next, convert the phantom files in the group navg>1 if they exist. They only ever join a
        # stream the rawdata already opened, since the header is never built from a prescan
        if writer is not None:
            for filepath in grouped_files.prescan_file_list:
                mrs = MRSdata()
                mrs.read_from_file(filepath)
                writer.write_data(generate_acquisition(mrs, None))
    finally:
        if writer is not None:
            writer.close()
    if writer is None:
        print(f"No data to convert for {grouped_files.meas_id}", file=sys.stderr)
        return False

    return True


def convert_tar_to_mrd(tar_path: Path, output_path: Path) -> bool:
    """
    Convert tar filed directory of single experiment epsi MRS data folder, exactly as --folder
    converts that directory unpacked

    Tyger hands a job its input buffer as a named FIFO, which is strictly sequential, so the
    archive is read forward once and held in memory: the parameter block sits after the raw data at
    EOF, the .SPR sidecar can follow the .MRD members in tar order, and the header needs the whole
    group before the first acquisition can be written.

    A Tyger job has one output buffer, so the tar holds one experiment and this writes one stream
    Args:
        - tar_path: tar archive of one experiment directory, or a FIFO carrying one
        - output_path: where to write the stream, which may also be a FIFO
    Returns:
        - True when a stream was written
    """
    # opening a FIFO for read blocks until the buffer sidecar opens the write end
    with open(tar_path, "rb") as tar_stream:
        tar_meas_id, spr_frequency, members = read_scan_tar(tar_stream)
    fallback = tar_meas_id or Path(output_path).stem
    if not tar_meas_id:
        # nothing to name the scan after: the archive holds no single root directory
        print(f"No single root directory in the tar, calling this scan {fallback}", file=sys.stderr)
    if not spr_frequency:
        print("No base frequency from a .SPR sidecar in the tar", file=sys.stderr)

    payloads = dict(members)
    grouped_files = mrs_organize.organize_members(members, fallback_meas_id=fallback)
    if grouped_files is None:
        # nothing grouped, but on Tyger the output buffer's FIFO still has to open and close, or
        # the sidecar is left blocked on a stream that never opens
        print(f"No data to convert in {tar_path}", file=sys.stderr)
        with open(output_path, "wb"):
            pass
        return False

    rep_count = grouped_files.nrepetitions
    writer: Optional[mrd.BinaryMrdWriter] = None
    rep_idx = 0
    try:
        for name in grouped_files.rawdata_file_list:
            mrs = MRSdata()
            mrs.parse_from_buffer(payloads[name])
            mrs.set_base_frequency(spr_frequency)  # the .MRD may defer its frequency to the sidecar
            if writer is None:
                print(f'Writing file at {output_path}', file=sys.stderr)
                writer = mrd.BinaryMrdWriter(str(output_path))
                writer.write_header(generate_header(mrs, grouped_files.meas_id, rep_count))
            writer.write_data(generate_acquisition(mrs, rep_idx, rep_count))
            rep_idx += mrs.nrepetitions
        if writer is not None:
            # a prescan is not a repetition of the acquisition it calibrates and is not counted in
            # rep_count, so it cannot carry a repetition index past the limit the header declared.
            # It sits at repetition 0 and is told apart by IS_NAVIGATION_DATA instead
            for name in grouped_files.prescan_file_list:
                mrs = MRSdata()
                mrs.parse_from_buffer(payloads[name])
                mrs.set_base_frequency(spr_frequency)
                writer.write_data(generate_acquisition(mrs, 0, rep_count))
    finally:
        if writer is not None:
            writer.close()
    if writer is None:
        print(f"No data to convert for {grouped_files.meas_id}", file=sys.stderr)
        with open(output_path, "wb"):
            pass
        return False
    return True


def main() -> int:
    """
    Convert MRS data to MRD2, in exactly one of three input modes, or with -w report where the echo
    peaks inside each gradient switch of what those modes resolve to and convert nothing.

    Every check that can stop the run lives here. Past this point a file that cannot be converted is
    reported and skipped, so the run ends with a status rather than a traceback
    Returns:
        - 0 when at least one stream was written, or one scan was reported on with -w. 1 when
          nothing was
    """
    parser = argparse.ArgumentParser(description="Convert MR Solutions MRS data to MRD2 format")
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("-t", "--tar", type=Path,
                      help="tar archive of one scan directory, or a FIFO carrying one")
    mode.add_argument("-f", "--folder", type=Path,
                      help="directory to walk for MRS .MRD files")
    parser.add_argument("-o", "--output", type=Path,
                        help="file or FIFO to write the MRD2 stream to. Required with --tar "
                             "(default: $OUTPUT_PIPE), optional with --input, unused with --folder")
    parser.add_argument("-n", "--dry-run", action="store_true",
                        help="with --folder only: report how the files group them, without converting")
    parser.add_argument("-w", "--window", action="store_true",
                        help="with --folder only: plot where the echo peaks in each gradient switch, "
                             "then take the drift along the switch train out and plot it again, "
                             "converting nothing")
    args = parser.parse_args()

    if args.tar and not args.tar.exists():
        parser.error(f"{args.tar} does not exist")
    if args.folder and not args.folder.is_dir():
        parser.error(f"{args.folder} is not a directory")
    if args.tar and (args.window or args.dry_run):
        parser.error(f"-t can only us to convert file")

    if args.folder:
        # convert folder allows
        print(f"Converting folder of single experiment folder {args.folder}", file=sys.stderr)
        written = convert_folder_to_mrd(args.folder, args.dry_run, args.window)
    elif args.tar:
        output = args.output or Path(os.environ.get("OUTPUT_PIPE", ""))
        if not str(output):
            parser.error("--tar needs --output, or $OUTPUT_PIPE set")
        print(f"Converting tar of single experiment folder {args.tar}", file=sys.stderr)
        written = convert_tar_to_mrd(args.tar, output)
    if not written:
        print("Nothing was converted", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
