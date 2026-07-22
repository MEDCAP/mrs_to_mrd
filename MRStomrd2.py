"""
Convert a list of MRS files to MRD2 format
EPSI:
- raw folder can be one experiment folder with multiple scans or one parent folder with multiple experiment folders
- skip scans with more than one measurement as phantom data
- upload raw folder to s3 bucket with filtered files by the frontend script based on file names > list of files or list of list of files
- based on file list count, dynamically create a tyger codespec with n inputs and 1 output
- codespec python script can have optional length of list as argument input
- upload buffer outputs to s3 bucket
- [FRONT] show completion on mrd

Spectral:
- raw folder can be one experiment folder or parent folder with multiple experiments
- upload raw folder to s3 bucket
- use the same codespec yml file for conversion

Dockerfile
- bash script to take command line arguments
- Dockerfile is all I need for mrd fork
- convers]

Edge cases:
- group of files can include spectral and epsi files, failing the choice of epsi_recon and spectral_recon functions. Upon upload, check for intended sequence
and prompt for parameters for reconstruction. For local run, use unifylevel number to distinguish epsi or fid
- 
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path
from typing import BinaryIO, Iterable, List, Union
import numpy as np

# mrd python package
import mrd
from MRSreader import MRSdata

def append_pulseq(mrs: MRSdata, pulseq_mrd_file: BinaryIO):
    """
    Read pulse   and sequence from mrd file
    TODO: Add modification of pulseq file depending on parameters on MRSData
    @param
        - mrs: instance of class MRSData that extracts parameters from raw MRS file
        - pulseq_mrd_file: mrd file that contain pulseq 
    @return
        - 
    """
    with mrd.BinaryMrdReader(pulseq_mrd_file) as reader:
        head = reader.header()  # ignore header
        for item in reader.read_data():
            if isinstance(item, mrd.StreamItem.PulseqDefinitions):
                # add fov from MRSData
                yield mrd.StreamItem.PulseqDefinitions(item)
            if isinstance(item, mrd.StreamItem.Blocks): # StreamItem.Blocks is a vector of mrd.Block
                yield mrd.StreamItem.Blocks(item)   
            if isinstance(item, mrd.StreamItem.Rf): # StreamItem.Rf is a mrd.RFEvent
                yield mrd.StreamItem.Rf(item)
            if isinstance(item, mrd.StreamItem.ArbitraryGradient): # StreamItem.ArbitraryGradient is a mrd.ArbitraryGradient
                yield mrd.StreamItem.ArbitraryGradient(item)
            if isinstance(item, mrd.StreamItem.TrapezoidalGradient): # StreamItem.TrapezoidalGradient is a mrd.TrapezoidalGradient
                yield mrd.StreamItem.TrapezoidalGradient(item)
            if isinstance(item, mrd.StreamItem.Adc): # StreamItem.Adc is a mrd.ADCEvent
                yield mrd.StreamItem.Adc(item)
            if isinstance(item, mrd.StreamItem.Shape):
                yield mrd.StreamItem.Shape(item)

    definitions = mrd.PulseqDefinitions()  
    # epsi base pulseq is wip and need to adapt based on slice selection parameters gz
    if mrs.pplfile.find("epsi") != -1:
        definitions.gradient_raster_time_ns = 1e5
        definitions.radiofrequency_raster_time_ns = mrs.sampleperiod * 100
        definitions.adc_raster_time_ns = 1e5
        definitions.block_duration_raster_ns = 1e5
        definitions.name = "MRS epsi"
        definitions.fov = mrd.ThreeDimensionalFloat(x=mrs.FOV, y=mrs.FOV, z=mrs.FOV)    # 3d share the same FOV
    return

def generate_pulseq(mrs: MRSdata) -> Iterable[mrd.StreamItem]:
    """
    Generate pulseq stream items based on MRSdata parameters
    Args:
        - mrs: MRSdata object containing parameters extracted from .MRD file
    Returns:
        - Iterable of mrd.StreamItem containing pulseq definitions, shapes, and events for the
    """
    TE = np.uint64(1.8E+5)              # 180us in ns: just an estimate for now, start acquiring 180us after 100us pulse start
    TR = np.uint64(mrs.tr * 1.0E+6)     # mrs.tr from ms to ns
    pulse_length = np.uint64(1.0E+5)    # 100us in ns: just an estimate for now
    definitions = mrd.PulseqDefinitions()
    # start, end, and duration are defined in multiples of raster times(seconds) as float
    definitions.gradient_raster_time = 1e-05                                        # typical value is 1e-05s=10us
    definitions.radiofrequency_raster_time = 1e-06                                  # typical value is 1e-06s=1us. mrs.sample period=100ns * 400
    definitions.adc_raster_time = 1e-07                                             # typical value is 1e-07s=100ns 
    definitions.block_duration_raster = 1e-05                                       # typical value is 1e-05s=10us
    definitions.name = "epsi"
    definitions.fov = mrd.ThreeDimensionalFloat(x=mrs.FOV, y=mrs.FOV, z=mrs.FOV)    # fov in m
    definitions.custom['TE_ms'] = 0.18                                                  # TE in ms, just an estimate for now, start acquiring 180us after pulse
    definitions.custom['TR_ms'] = mrs.tr
    definitions.custom['acq_start_time_ns'] = str(mrs.acqstarttime * 100)  # acqstarttime in units of 100ns converted to ns
    yield mrd.StreamItem.PulseqDefinitions(definitions)

    # define shape of the RF pulse uncompressed as currently pulseq-mrd conversion does not support compression
    rf_amp_shape = mrd.Shape()
    rf_amp_shape.id = 1
    # pulse length=100us / dt=10us = 10samples
    rf_amp_shape.num_samples = pulse_length // definitions.radiofrequency_raster_time_ns
    # shape data is normalized to [-1, 1] and amplitude is set in RFPulseEvent field
    rf_amp_shape.data = np.ones(rf_amp_shape.num_samples, dtype=np.float64)
    yield mrd.StreamItem.Shape(rf_amp_shape)

    # define shape of the RF pulse uncompressed as currently pulseq-mrd conversion does not support compression
    rf_phase_shape = mrd.Shape()
    rf_phase_shape.id = 2
    # pulse length=100us / dt=10us = 10samples
    rf_phase_shape.num_samples = pulse_length // definitions.radiofrequency_raster_time_ns
    # phase is unknown so set to zeros
    rf_phase_shape.data = np.zeros(rf_amp_shape.num_samples, dtype=np.float64)
    yield mrd.StreamItem.Shape(rf_phase_shape)

    # define RF event to specify amplitude and phase, offsets
    rf = mrd.RFEvent()
    rf.id = 1                       # correspond to block.rf=1
    rf.amp = float(1E+5)            # peak amplitude in (Hz) in float value is guessed
    rf.mag_id = rf_amp_shape.id
    rf.phase_id = rf_phase_shape.id
    rf.time_id = 0                  # time_id=0 to use radiofrequency_raster_time_ns
    rf.center_ns = pulse_length // 2 # center of the pulse in ns
    rf.delay_ns = 0                 # delay before rf pulse start
    rf.freq_ppm = 0                 # freq offset in ppm relative to main system's freq
    rf.phase_ppm = 0                # phase offset in rad/MHz proportional to main system's freq
    rf.freq_offset = 0              # freq offset in Hz
    rf.phase_offset = 0             # phase offset in rad
    rf.use = mrd.RFPulseUse.EXCITATION
    yield mrd.StreamItem.Rf(rf)

def generate_acquisition(mrs: MRSdata, head: mrd.Header, idx: int) -> Iterable[mrd.StreamItem]:
    """
    Extract acquisition from MRS data file
    Args:
        - mrs: MRSdata object containing rawdata and parameters extracted from .MRD file
               mrs.rawdata shape=(samples, views, sliceviews, slices, echoes, nex)
               EPSI: (1280, 8, 1, 1, 1, 1)
        - head: mrd.Header object containing header information
        - idx: index of the current file in the grouped EPSI files
    Returns:
        - Iterable of mrd.StreamItem.Acquisition containing acquisition data and parameters for each acquisition
    """
    TE = np.uint64(1.8E+5)           # in ns: just an estimate for now, start acquiring 180us after pulse
    TR = np.uint64(mrs.tr * 1.0E+6)  # convert mrs.tr from ms to ns
    # encode acquisition for EPSI sequence
    if "epsi" in mrs.pplfile:
        if mrs.rawdata.shape[2:] != (1, 1, 1, 1):
            raise ValueError(f"Expected rawdata shape (samples, views, 1, 1, 1, 1) but got {mrs.rawdata.shape}")
        elif mrs.rawdata.shape[0] == 1 or mrs.rawdata.shape[1] == 1:
            raise ValueError(f"Expected rawdata shape (samples, views, 1, 1, 1, 1) but got {mrs.rawdata.shape}")
        for iview in range(mrs.rawdata.shape[1]):                                    # rawdata.shape[1]=views
            acq = mrd.Acquisition()
            if mrs.navg > 1:                                                         # mrs.avg>1 is a phantom data
                acq.head.flags |= mrd.AcquisitionFlags.IS_NAVIGATION_DATA
            acq.head.idx.average = mrs.navg                                          # phantom scans have avg>1
            # skip phantom data of idx for repetition and scan_counter 
            phantom_idx_list = [x.value for x in head.user_parameters.user_parameter_long]
            # repetition and scan_counter encoded only for nonphantom data
            if idx not in phantom_idx_list:
                acq.head.idx.repetition = idx - len([x for x in phantom_idx_list if x<idx])
                acq.head.scan_counter = acq.head.idx.repetition * mrs.rawdata.shape[1] + iview
                acq.head.idx.kspace_encode_step_1 = iview
            if iview == 0:
                acq.head.flags |= mrd.AcquisitionFlags.FIRST_IN_PHASE
                if idx == 0:           
                    acq.head.flags |= mrd.AcquisitionFlags.FIRST_IN_REPETITION
            elif iview == mrs.rawdata.shape[1] - 1:                                  # rawdata.shape[1]=views
                acq.head.flags = mrd.AcquisitionFlags.LAST_IN_PHASE
                if acq.head.idx.repetition == head.encoding[0].encoding_limits.repetition.maximum:
                    acq.head.flags |= mrd.AcquisitionFlags.LAST_IN_REPETITION        # last non-phantom acq in repetition
            acq.head.idx.phase = iview
            acq.data = np.expand_dims(np.squeeze(mrs.rawdata)[:, iview], axis=0)     # mrs.rawdata.shape=(samples,views) acq.data.shape=(coils=1, samples)
            pulse_start = np.uint64(mrs.acqstarttime * 100) + np.uint64(iview * TR)  # mrs.acqstarttime in units of 100ns converted to ns
            acq.head.acquisition_time_stamp_ns = pulse_start + TE
            acq.head.idx.contrast = mrs.nswitch                                      # contrast=echo number in multi_echo=nswitch
            totalppswitch = round(mrs.rawdata.shape[0] / mrs.nswitch)                # samples/nswitch=npts_per_echo + 2 * npts_per_ramp
            acq.head.discard_pre = int((totalppswitch - mrs.nppswitch) / 2)          # npts_per_ramp to discard
            acq.head.discard_post = acq.head.discard_pre
            acq.head.sample_time_ns = mrs.sampleperiod * 100                         # sampleperiod in units of 100ns
            yield mrd.StreamItem.Acquisition(acq)

    # MRS->mrd conversion for spectral sequence
    elif "1pul" in mrs.pplfile or "fid" in mrs.pplfile: # spectral sequence ppl file have two names for EVO1 and EVO2
        nrep = mrs.rawdata.shape[5]
        for iview in range(nrep):
            # encode acquisition field
            acq = mrd.Acquisition()
            # MRS is single channel acquisition. Add new empty coil dimension since acq.data.shape=(coils,samples)
            acq.data = np.expand_dims(mrs.rawdata[:, 0, 0, 0, 0, iview], (0))
            pulse_start = np.uint64(mrs.acqstarttime * 100) + np.uint64(iview * TR) # mrs.acqstarttime in units of 100ns converted to ns
            acq.head.idx.repetition = iview
            acq.head.acquisition_time_stamp_ns = pulse_start + TE
            acq.head.idx.contrast = 1
            acq.head.discard_pre = 0
            acq.head.discard_post = 0
            acq.head.sample_time_ns = mrs.sampleperiod * 100                  # convert to ns (sampleperiod is in units of 100ns)
            acq.phase = np.zeros((mrs.rawdata.shape[0]), dtype=np.float32)    # acquisition phase array set to zeros
            yield mrd.StreamItem.Acquisition(acq)
    else:
        raise Exception(f"Invalid sequence protocol {mrs.pplfile}")

def make_header(mrs: MRSdata, meas_id: str, rep_count: int, phantom_idx_list: List[int]) -> mrd.Header:
    """
    Fill in MRD header based on mrs parameters fields and meas_id from filename
    Args:
        - mrs: MRSdata object containing parameters extracted from .MRD file
        - meas_id: string extracted from filename to use as measurement id in header
        - rep_count: total count of data in this repetition
        - phantom_idx_list: list of idx where phantom data with navg>1 show up in repetition
    Returns:
        - mrd.Header object with filled in fields based on mrs parameters and meas_id
    """
    # make mrd2 header. For now only filling in sequence name but some day should do more
    header = mrd.Header()
    s = mrd.SubjectInformationType()
    s.patient_id = meas_id                  # e.g.) cirrhrat_1_1, PYR_HepG22_00_00_0000.mrs
    header.subject_information = s
    
    meas = mrd.MeasurementInformationType()
    meas.sequence_name = mrs.pplfile
    meas.relative_table_position = mrd.ThreeDimensionalFloat(x=mrs.FOVoff[0]*1e3, y=mrs.FOVoff[1]*1e3, z=mrs.FOVoff[2]*1e3)
    meas.measurement_id = meas_id
    meas.protocol_name = meas_id.split("_")[0]
    header.measurement_information = meas

    header.experimental_conditions.h1resonance_frequency_hz = mrs.basefreq

    e = mrd.EncodingSpaceType()
    e.matrix_size = mrd.MatrixSizeType(x=mrs.rawdata.shape[0], y=mrs.rawdata.shape[1], z=mrs.rawdata.shape[3])  # samples, views, slices
    e.field_of_view_mm = mrd.FieldOfViewMm(x=mrs.FOV*1e3, y=mrs.FOV*1e3, z=0)                                   # mrs.FOV in m converted to mm

    limits_rep = mrd.LimitType()
    limits_rep.minimum = 0
    limits_rep.maximum = rep_count - len(phantom_idx_list) - 1  # max repetition as index
    
    # encode size of each dimension of rawdata
    limits_kspace_encode_step_0 = mrd.LimitType()
    limits_kspace_encode_step_0.maximum = mrs.rawdata.shape[0] - 1

    limits_kspace_encode_step_1 = mrd.LimitType()
    limits_kspace_encode_step_1.maximum = mrs.rawdata.shape[1] - 1     # max views as index

    limits_slice = mrd.LimitType()
    limits_slice.maximum = mrs.rawdata.shape[3] - 1     # max slices as index

    limits = mrd.EncodingLimitsType()
    limits.kspace_encoding_step_0 = limits_kspace_encode_step_0
    limits.kspace_encoding_step_1 = limits_kspace_encode_step_1
    limits.repetition = limits_rep
    limits.slice = limits_slice
    enc = mrd.EncodingType()
    enc.encoded_space = e
    enc.encoding_limits = limits
    header.encoding.append(enc)

    user_params = mrd.UserParametersType()
    for idx in phantom_idx_list:
        user_params.user_parameter_long.append(mrd.UserParameterLongType(name="phantom_idx", value=idx))
    header.user_parameters = user_params
    return header

def collect_mrd_files(folder: Path) -> List[Path]:
    """
    Find all mrs .MRD files in the input folder and its subdirectories recursively
    Args:
        folder: Path object to the root directory
    Returns:
        List of Path objects
    """
    mrd_filepath_list: List[Path] = []
    file_extension = ".MRD"
    # recursively find all file paths with .MRD extension in the rootdir
    for entry in folder.iterdir():
        if entry.is_dir():
            mrd_filepath_list.extend(collect_mrd_files(entry))
            continue
        if file_extension in entry.name:
            mrd_filepath_list.append(entry)
        else:
            continue
    return mrd_filepath_list

def group_mrd_files(folder: Path, unifylevel: int) -> List[List[Path]]:
    """
    Group mrd files in the folder based on unify level
        spectral data=unify_level=1
            example folder: /pyruvate_data/PYR_HepG22.mrs/12345_000_0.MRD
            single file for each group
        epsi data=unify_level=3
            example folder: /cirrhrat_data/cirrhrat_0_1/epsi/12345/12345_000_0.MRD
            extract paths of 3 directories above: file.parts[:-(unifylevel-1)]=(cirrhrat_data, cirrhrat_0_1, epsi)
            group files with the same paths
    Args:
        folder: Path object to the input folder
        unifylevel: integer specifying the number of levels to unify
    Returns:
        List of lists of Path objects
    """
    mrd_filepath_list = collect_mrd_files(folder)   # list of all MRS filepaths in the folder   
    mrd_file_groups: List[List[Path]] = []          # list of list as grouped mrd files
    for file in mrd_filepath_list:
        is_grouped = False
        for group in mrd_file_groups:
            if file.parts[:-(unifylevel-1)] == group[0].parts[:-(unifylevel-1)]:
                group.append(file)
                is_grouped = True
        if not is_grouped:
            mrd_file_groups.append([file])
    print(f"Grouped {len(mrd_filepath_list)} files into {len(mrd_file_groups)} groups", file=sys.stderr)
    return mrd_file_groups

def convert_mrs_folder_to_mrd(folder: Path, unifylevel: int) -> None:
    """
    Convert a list of MRS files grouped by the same EPSI measurements into MRD format with header and acquisition
    Edge case: there could be an error fid file among the epsi files. Use unifylevel and mrs.pplfile to skip it 
    Spectral data: unify_level=1 for /pyruvate_data/PYR_HepG22.mrs/12345_000_0.MRD
        meas_id = PYR_HepG22.mrs
        output_dir = /pyruvate_data/PYR_HepG22.mrs/raw.mrd2
    EPSI data: unify_level=3 for /cirrhrat_data/cirrhrat_0_1/epsi/12345/12345_000_0.MRD
        meas_id = cirrhrat_0_1
        output_dir = /cirrhrat_data/cirrhrat_0_1/raw.mrd2
    Args:
        folder: Path object to the input folder containing MRS files
        unifylevel: integer specifying the number of layers from the edge of filepath to group files
    Returns:
        None
    """
    mrs = MRSdata()
    mrd_file_groups = group_mrd_files(folder, unifylevel)
    for group in mrd_file_groups:
        meas_id = group[0].parts[-(unifylevel+1)]                                       # e.g.) meas_id=cirrhrat_0_1
        raw_filepath = os.path.join(Path(*group[0].parts[:-unifylevel]), "raw.mrd2")    # e.g.) cirrhrat_data/cirrhrat_0_1/raw.mrd2
        with mrd.BinaryMrdWriter(raw_filepath) as writer:
            # iterate through group to count the number of phantoms
            phantom_idx_list: List[int] = []                                            # index of phantom data in each file group
            for j, filepath in enumerate(group):
                mrs.mread3d(filepath)
                # if files in the group are inconsistent between epsi and fid with unifylevel, remove it
                if unifylevel == 3 and "1pul" in mrs.pplfile or "fid" in mrs.pplfile:
                    group.remove(filepath)
                if unifylevel == 1 and "epsi" in mrs.pplfile:
                    group.remove(filepath)
                # if navg>1 the file is phantom data
                if mrs.navg > 1:                                                # navg>1 is phantom
                    phantom_idx_list.append(j)
                    print(f"Phantom data found at {filepath}", file=sys.stderr)
            print(f"Found {len(phantom_idx_list)} phantom data among {len(group)} files", file=sys.stderr)
            for idx, filepath in enumerate(group):
                mrs.mread3d(filepath)
                if idx==0:
                    if mrs.navg == 1:
                        header = make_header(mrs, meas_id, len(group), phantom_idx_list)
                        writer.write_header(header)                                 # write header for the first non-phantom data in the group
                    else:
                        raise ValueError(f"First file in group {filepath} is phantom data")            
                writer.write_data(generate_acquisition(mrs, header, idx))
                # writer.write_data(generate_pulseq(mrs))                       # generate pulseq field, currently not writing pulseq to save time as it is not used for reconstruction
    
def convert_mrs_file_to_mrd(input: Path,
                            output: Path):
    """
    Convert single MRS .MRD file to mrd2 format
    """
    mrs = MRSdata()
    with mrd.BinaryMrdWriter(output) as writer:
        mrs.mread3d(input)
        meas_id = input.parent.name
        rep_count = mrs.rawdata.shape[5]
        header = make_header(mrs, meas_id, rep_count, [])
        writer.write_header(header)
        writer.write_data(generate_acquisition(mrs, header, 0))

def main() -> int:
    """
    In the specified folder, listen for new file and incrementally write to mrd2 file
    - f/--folder: folder containing MRS data files
    - u/--unifylevel: directory levels to unify when grouping files (default: 1)
    - i/--input: single MRS .MRD file input for conversion to MRD2 format, if folder is not passed
    """

    parser = argparse.ArgumentParser(description='Convert MRS data folder to MRD2 format')
    parser.add_argument("-f", "--folder", type=Path, required=False,
                        help="Directory containing MRS data files")
    parser.add_argument("-u", "--unifylevel", type=int, required=False, default=1,
                        help="Directory levels to unify when grouping files (default: 1)")
    parser.add_argument("-i", "--input", type=Path, required=False,
                        help="Single MRS .MRD file input for conversion to MRD2 format")
    args = parser.parse_args()

    # if folder is passed, collect files. Otherwise, convert take MRS file input
    if args.folder and args.input:
        raise ValueError("Cannot specify both --folder and --input")
    elif args.folder:
        if not args.folder.is_dir():
            raise ValueError(f"{args.folder} is not a directory")
        print(f"Convert folder {args.folder} with unify level {args.unifylevel}", file=sys.stderr)
        convert_mrs_folder_to_mrd(args.folder, args.unifylevel)
    elif args.input:
        if not args.input.is_file():
            raise ValueError(f"{args.input} is not a file")
        print(f"Convert single input file {args.input} to mrd2", file=sys.stderr)
        output = os.path.join(args.input.parent,"raw.mrd")
        convert_mrs_file_to_mrd(args.input, output)
    else:
        raise ValueError("Either --folder or --input must be specified")
    return 0

# -u 3 consolidates the files as appropriate for EPSI. 
# spectral data like Bukola's and David's uses -u 1
if __name__ == "__main__":
    raise SystemExit(main())