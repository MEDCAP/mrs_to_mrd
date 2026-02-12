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
- convers
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Iterable, List, Sequence, BinaryIO, Type

import numpy as np

# mrd python package
# https://github.com/MEDCAP/mrd-fork/tree/dev# is added at root of this repository as a git submodule
# path to mrd python package is 'root/mrd-fork/python'
MRD_PYTHON_PATH = Path(__file__).resolve().parent / "mrd-fork" / "python"
sys.path.insert(0, str(MRD_PYTHON_PATH))
import mrd
from MRSreader import MRSdata

MRSCONVERT_DEBUG = True

def append_pulseq(mrs: MRSData, pulseq_mrd_file: binaryIO):
    """
    Append mrs pulse amplitude and phase shape on mrd dataset as a separate file
    TODO: Add modification of pulseq file depending on parameters on MRSData
    @param
        - mrs: instance of class MRSData that extracts parameters from raw MRS file
        - amp_file: npy file containing pulse shape in 
        - phase_file: npy file containing pulse phase in radians
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
    if mrs.pplfile.find("epsi") >= 0:
        definitions.gradient_raster_time_ns = 1e5
        definitions.radiofrequency_raster_time_ns = mrs.sampleperiod * 100
        definitions.adc_raster_time_ns = 1e5
        definitions.block_duration_raster_ns = 1e5
        definitions.name = "MRS epsi"
        definitions.fov = mrd.ThreeDimensionalFloat(x=mrs.FOV, y=mrs.FOV, )


def generate_pulseq_acquisition(mrs: MRSdata, measurement_index: int) -> Iterable[mrd.StreamItem]:
    """
    Define pulseq fields and extract acquisition from MRSdata structure.
    """
    # define pulseq definitions with pseudo rf pulses
    # all time units are aligned to be in ns
    pulse_length = np.uint64(1.0E+5)   # 100us in ns: guess pulelength to 100us
    TE = np.uint64(1.8E+5)         # 180us in ns: just an estimate for now, start acquiring 180us after 100us pulse start
    TR = np.uint64(mrs.tr * 1.0E+6)  # mrs.tr from ms to ns
    definitions = mrd.PulseqDefinitions()
    # start, end, and duration are measured in multiples of raster times
    definitions.gradient_raster_time_ns = 1                             # typical pulseq value is 1e-05s=10us
    definitions.radiofrequency_raster_time_ns = mrs.sampleperiod * 100    # sample period is units of 100ns
    definitions.adc_raster_time_ns = 1                                  # typical pulseq value is 1e-07s=100ns 
    definitions.block_duration_raster_ns = 1                            # typical pulseq value is 1e-05s=10us
    definitions.name = "MRS epsi"                       
    definitions.fov = mrd.ThreeDimensionalFloat(x=mrs.FOV, y=mrs.FOV, z=mrs.FOV)    # optional three dimensionak float
    definitions.custom['TE_ns'] = str(TE)
    definitions.custom['TR_ns'] = str(TR)
    definitions.custom['acq_start_time_ns'] = str(mrs.acqstarttime * 100)  # acqstarttime in units of 100ns converted to ns
    definitions.custom['pulse_length_ns'] = str(pulse_length)               
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

    # encode pulse events and acquisition in a time-series for EPSI sequence
    if mrs.pplfile.find("epsi") >= 0:
        nPE = mrs.rawdata.shape[1]  # number of phase encoding
        blocks = []
        for iacq in range(nPE):
            # encode pulse events and acquisition in a time-series for each phase-encoding
            block1 = mrd.Block()
            block1.id = iacq*2 + 1      # id is non-zero unique values, starting at 1 for rf pulse, 2 for adc period, 3 for next pulse
            block1.duration = pulse_length  # pulse length in units of definitions.block_duration_raster_ns
            block1.rf = 1                # rf pulse id
            block1.gx = 0                # no gradients for now
            block1.gy = 0
            block1.gz = 0
            block1.adc = 0
            block1.ext = 0               # extension to save LABEL for counters and flags, TRIGGERS
            blocks.append(block1)
            
            # to next rf pulse
            block2 = mrd.Block()
            block2.id = iacq*2 + 2               # increment id from the previous value
            block2.duration = TR - pulse_length
            block2.rf = 0
            block2.gx = 0
            block2.gy = 0
            block2.gz = 0
            block2.adc = 0
            block2.ext = 0
            blocks.append(block2)
        
        # Yield all blocks at once
        yield mrd.StreamItem.Blocks(blocks)
        
        # Now yield acquisitions
        for iacq in range(nPE):
            # encode acquisition field
            acq = mrd.Acquisition()
            if(iacq == 0):
                acq.head.flags = mrd.AcquisitionFlags.FIRST_IN_PHASE
            if(iacq == mrs.rawdata.shape[1] - 1):
                acq.head.flags = mrd.AcquisitionFlags.LAST_IN_PHASE
            acq.data = np.transpose(np.expand_dims(mrs.rawdata[:, iacq, 0, 0, 0, 0], (0)))
            acq.head.acquisition_center_frequency = mrs.basefreq
            acq.head.idx.phase = iacq
            # pulseq only defines duration
            # timestamp of each pulse start is recalculated here
            pulse_start = np.uint64(definitions.custom['acq_start_time_ns']) + np.uint64(iacq * TR)
            acq.head.acquisition_time_stamp_ns = pulse_start + TE
            acq.head.idx.contrast = mrs.nswitch
            totalppswitch = int(mrs.rawdata.shape[0] / mrs.nswitch + 0.7)
            acq.head.discard_pre = int((totalppswitch - mrs.nppswitch) / 2)
            acq.head.discard_post = acq.head.discard_pre
            acq.head.sample_time_ns = mrs.sampleperiod * 100                  # sampleperiod is in units of 100ns
            acq.phase = np.zeros((mrs.rawdata.shape[0]), dtype=np.float32)    # acquisition phase array set to zeros
            yield mrd.StreamItem.Acquisition(acq)

    # MRS->mrd conversion for spectral sequence
    elif mrs.pplfile.find("1pul") >= 0:
        nrep = mrs.rawdata.shape[5]  # number of acquisitions per image file (phase-encodes)
        blocks = []
        for iacq in range(nrep):
            # encode pulse events and acquisition in a time-series for 1pul sequence
            block1 = mrd.Block()
            block1.id = iacq*2 + 1      # id is non-zero unique values, starting at 1 for rf pulse, 2 for adc period, 3 for next pulse
            block1.duration = pulse_length  # pulse length in units of definitions.block_duration_raster_ns
            block1.rf = 1                # rf pulse id
            block1.gx = 0                # no gradients for now
            block1.gy = 0
            block1.gz = 0
            block1.adc = 0
            block1.ext = 0               # extension to save LABEL for counters and flags, TRIGGERS
            blocks.append(block1)
            
            # to next rf pulse
            block2 = mrd.Block()
            block2.id = iacq*2 + 2               # increment id from the previous value
            block2.duration = TR - pulse_length
            block2.rf = 0
            block2.gx = 0
            block2.gy = 0
            block2.gz = 0
            block2.adc = 0
            block2.ext = 0
            blocks.append(block2)
        
        # Blocks stream is a list of mrd.Block() 
        yield mrd.StreamItem.Blocks(blocks)
        
        # Now yield acquisitions
        for iacq in range(nrep):
            # encode acquisition field
            acq = mrd.Acquisition()
            acq.data = np.transpose(np.expand_dims(mrs.rawdata[:, 0, 0, 0, 0, iacq], (0)))
            acq.head.acquisition_center_frequency = mrs.basefreq
            pulse_start = np.uint64(definitions.custom['acq_start_time_ns']) + np.uint64(iacq * TR)
            acq.head.idx.repetition = iacq
            acq.head.acquisition_time_stamp_ns = pulse_start + TE
            acq.head.idx.contrast = 1
            acq.head.discard_pre = 0
            acq.head.discard_post = 0
            acq.head.sample_time_ns = mrs.sampleperiod * 100                  # convert to ns (sampleperiod is in units of 100ns)
            acq.phase = np.zeros((mrs.rawdata.shape[0]), dtype=np.float32)    # acquisition phase array set to zeros
            yield mrd.StreamItem.Acquisition(acq)

def make_header(mrs: MRSdata, meas_id: str) -> mrd.Header:
    # make mrd2 header. For now only filling in sequence name but some day should do more
    h = mrd.Header()
    meas = mrd.MeasurementInformationType()
    meas.sequence_name = mrs.pplfile
    # this may be a misuse of the relative table position but I couldn't find FOV offset anywhere
    meas.relative_table_position = mrd.ThreeDimensionalFloat(x=mrs.FOVoff[0], y=mrs.FOVoff[1], z=mrs.FOVoff[2])
    meas.measurement_id = meas_id
    h.measurement_information = meas
    # This is definitely a misuse of h1 resonance frequency for non-1H acquisitions also saved as acq center freq
    exp = mrd.ExperimentalConditionsType()
    exp.h1_resonance_frequency_hz = mrs.basefreq
    h.experimental_condition = exp
    e = mrd.EncodingType()
    e.encoded_space.field_of_view_mm.x = mrs.FOV
    e.encoded_space.field_of_view_mm.y = mrs.FOV
    e.encoded_space.field_of_view_mm.z = mrs.FOV
    h.encoding.append(e)    # encoding as a list
    return(h)

def collect_mrd_files(rootdir: Path, mrs: MRSdata) -> List[Path]:
    """
    Find all .MRD files in the root directory and its subdirectories.
    Assumption: Skip .MRD file with more than 1 average as it is considered as phantom data.
    Args:
        rootdir: Path object to the root directory
        mrs: MRSdata object
    Returns:
        List of Path objects
    """
    mrd_filepath_list: List[Path] = []
    # recursively find all file paths with .MRD extension in the rootdir
    for entry in rootdir.iterdir():
        if entry.is_dir():
            mrd_filepath_list.extend(collect_mrd_files(entry, mrs))
            continue
        if ".MRD" in entry.name:
            mrs.mread3d(str(entry)) # extract MRS data from .MRD file into mrs object to read number of averages
            print(f"Filepath {entry} with {mrs.navg} avg", file=sys.stderr)
            # skip the file with more than 1 average as phantom data
            if mrs.navg > 1:
                print(f"Skipping {entry} with {mrs.navg} avg as phantom data", file=sys.stderr)
                continue
            mrd_filepath_list.append(entry)
        else:
            # if the entry is not a .MRD file, continue
            continue
    return mrd_filepath_list

def group_mrd_files(rootdir: Path, unifylevel: int, mrs: MRSdata) -> List[List[Path]]:
    """
    In the rootdir, make a list of MRD files that have the same parental paths by unify level.
    For example, if unifylevel is 1, and the rootdir is /A/B/C/protocolA/, the list of lists will be:
    rootdir: epsi_kidney_data
    If unify level is 3, 
    condition: 3 dir above ischemia_27_1
    If unify level is 1, 1 dir above is 8262
    Args:
        rootdir: Path object to the root directory
        unifylevel: integer specifying the number of levels to unify
        mrs: MRSdata object
    Returns:
        List of lists of Path objects
    """
    if not rootdir.is_dir():
        raise ValueError(f"Root directory {rootdir} is not a directory")
    mrd_filepath_list = collect_mrd_files(rootdir, mrs)
    mrd_file_groups: List[List[Path]] = []
    for path in mrd_filepath_list:
        parts = path.parts
        added_to_group = False
        # add same folder name filepath into the same group
        for group in mrd_file_groups:
            group_parts = group[0].parts
            if all(parts[i] == group_parts[i] for i in range(len(parts) - unifylevel)):
                group.append(path)
                added_to_group = True
                break
        # if new folder name is found, append a new group
        if not added_to_group:
            mrd_file_groups.append([path])
    return mrd_file_groups

def convert_mrs_folder_to_mrd(basedir: Path, unifylevel: int) -> None:
    """
    Frontend typescript organizes the files into groups of same protocol
    Presigned URL sent to upload a folder of MRD files combined together into same protocol
    S3 bucket 
    """
    mrs = MRSdata()
    mrd_file_groups = group_mrd_files(basedir, unifylevel, mrs)
    for group in mrd_file_groups:
        meas_id = group[0].parts[-(unifylevel + 1)] # protocol folder name e.g.) ischemia_27_1
        output_dir = Path(*group[0].parts[: -unifylevel])
        print(f"grouping {len(group)} files into {output_dir}", file=sys.stderr)
        writer = mrd.BinaryMrdWriter(str(output_dir / "raw.mrd2"))
        for index, file_path in enumerate(group):
            # write the header once for the entire group after continue
            # re-read the data to set the correct header for each group
            mrs.mread3d(str(file_path))
            if index == 0:
                header = make_header(mrs, meas_id)
                writer.write_header(header)
            writer.write_data(generate_pulseq_acquisition(mrs, index))
        writer.close()


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description='Convert MRS data folder to MRD2 format')
    # required directory path with MRS data files as Path type
    parser.add_argument('-f', '--folder', type=Path, required=False,
                        help='Base directory containing MRS data files')
    parser.add_argument('-u', '--unifylevel', type=int, required=False, default=1,
                        help='Directory levels to unify when grouping files (default: 1)')
    args = parser.parse_args()
    # if folder is passed, collect files. Otherwise, take single file input
    print(f"running with base director {args.folder}", file=sys.stderr)
    print(f"unify level set to {args.unifylevel}", file=sys.stderr)
    convert_mrs_to_mrd(args.folder, args.unifylevel)
    return 0

# -u 3 consolidates the files as appropriate for EPSI. 
# spectral data like Bukola's and David's uses -u 1
if __name__ == "__main__":
    raise SystemExit(main())