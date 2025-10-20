"""
python script to convert MRS folder to MRD files in place
"""

import numpy as np
import sys
import os
import argparse
from statistics import mode

# mrd python package
# https://github.com/MEDCAP/mrd-fork/tree/dev# is added at root of this repository as a git submodule
# path to mrd python package is 'root/mrd-fork/python' 
sys.path.insert(0, 'mrd-fork/python')
import mrd
from MRSreader import MRSdata

mrs = MRSdata()

def generate_pulseq_acquisition(g, ide):
    '''
    Define pulseq fields and extract acquisition from MRSdata structure g 
    '''
    # define pulseq definitions with pseudo rf pulses
    # all time units are aligned to be in ns
    pulse_length = np.uint64(1.0E+5)   # 100us in ns: guess pulelength to 100us
    TE = np.uint64(1.8E+5)         # 180us in ns: just an estimate for now, start acquiring 180us after 100us pulse start
    TR = np.uint64(g.tr * 1.0E+6)  # g.tr from ms to ns
    definitions = mrd.PulseqDefinitions()
    # start, end, and duration are measured in multiples of raster times
    definitions.gradient_raster_time_ns = 1                             # typical pulseq value is 1e-05s=10us
    definitions.radiofrequency_raster_time_ns = g.sampleperiod * 100    # sample period is units of 100ns
    definitions.adc_raster_time_ns = 1                                  # typical pulseq value is 1e-07s=100ns 
    definitions.block_duration_raster_ns = 1                            # typical pulseq value is 1e-05s=10us
    definitions.name = "MRS epsi"                       
    definitions.fov = g.FOV
    definitions.custom['TE_ns'] = str(TE)
    definitions.custom['TR_ns'] = str(TR)
    definitions.custom['acq_start_time_ns'] = str(g.acqstarttime * 100)  # acqstarttime in units of 100ns converted to ns
    definitions.custom['pulse_length_ns'] = str(pulse_length)               
    yield mrd.StreamItem.PulseqDefinitions(definitions)

    # define shape of the RF pulse uncompressed as currently pulseq-mrd conversion does not support compression
    rf_amp_shape = mrd.shape()
    rf_amp_shape.id = 1
    # pulse length=100us / dt=10us = 10samples
    rf_amp_shape.num_samples = pulse_length // definitions.radiofrequency_raster_time_ns
    # shape data is normalized to [-1, 1] and amplitude is set in RFPulseEvent field
    rf_amp_shape.data = np.ones(rf_amp_shape.num_samples, dtype=np.float64)
    yield mrd.StreamItem.Shape(rf_amp_shape)

    # define shape of the RF pulse uncompressed as currently pulseq-mrd conversion does not support compression
    rf_phase_shape = mrd.shape()
    rf_phase_shape.id = 2
    # pulse length=100us / dt=10us = 10samples
    rf_phase_shape.num_samples = pulse_length // definitions.radiofrequency_raster_time_ns
    # phase is unknown so set to zeros
    rf_phase_shape.data = np.zeros(rf_amp_shape.num_samples, dtype=np.float64)
    yield mrd.StreamItem.Shape(rf_phase_shape)

    # define RF event to specify amplitude and phase, offsets
    rf = mrd.RFPulseEvent()
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
    yield mrd.StreamItem.RFPulseEvent(rf)

    # encode pulse events and acquisition in a time-series for EPSI sequence
    if(g.pplfile.find('epsi') >= 0):
        [npts, nPE] = g.rawdata.shape   # number of points per acquisition and phase-encoding
        for iacq in range(nPE):
            # encode pulse events and acquisition in a time-series for each phase-encoding
            block = mrd.Block()
            block.id = iacq*2 + 1      # id is non-zero unique values, starting at 1 for rf pulse, 2 for adc period, 3 for next pulse
            block.duration = pulse_length  # pulse length in units of definitions.block_duration_raster_ns
            block.rf = 1                # rf pulse id
            block.gx = 0                # no gradients for now
            block.gy = 0
            block.gz = 0
            block.adc = 0
            block.ext = 0               # extension to save LABEL for counters and flags, TRIGGERS
            yield mrd.StreamItem.Block(block)
            # to next rf pulse
            block = mrd.Block
            block.id += 1               # increment id from the previous value
            block.duration = TR - pulse_length
            block.rf = 0
            block.gx = 0
            block.gy = 0
            block.gz = 0
            block.adc = 0
            block.ext = 0
            yield mrd.StreamItem.Block(block)

            # encode acquisition field
            acq = mrd.Acquisition()
            if(iacq == 0):
                acq.head.flags = mrd.AcquisitionFlags.FIRST_IN_PHASE
            if(iacq == g.rawdata.shape[1] - 1):
                acq.head.flags = mrd.AcquisitionFlags.LAST_IN_PHASE
            acq.data = np.transpose(np.expand_dims(g.rawdata[:, iacq, 0, 0, 0, 0], (0)))
            acq.head.acquisition_center_frequency = g.basefreq
            acq.head.idx.phase = iacq
            # pulseq only defines duration
            # timestamp of each pulse start is recalculated here
            pulse_start = np.uint64(definitions.custom['acq_start_time_ns']) + np.uint64(iacq * TR)
            acq.head.acquisition_time_stamp_ns = pulse_start + TE
            acq.head.idx.contrast = g.nswitch
            totalppswitch = int(g.rawdata.shape[0] / g.nswitch + 0.7)
            acq.head.discard_pre = int((totalppswitch - g.nppswitch) / 2)
            acq.head.discard_post = acq.head.discard_pre
            acq.head.sample_time_ns = g.sampleperiod * 100                  # sampleperiod is in units of 100ns
            acq.phase = np.zeros((g.rawdata.shape[0]), dtype=np.float32)    # acquisition phase array set to zeros
            yield mrd.StreamItem.Acquisition(acq)

    # MRS->mrd conversion for spectral sequence
    elif(g.pplfile.find('1pul') >= 0):
        nrep = g.rawdata.shape[5]  # number of acquisitions per image file (phase-encodes)
        for iacq in range(nrep):
            # encode pulse events and acquisition in a time-series for 1pul sequence
            block = mrd.Block()
            block.id = iacq*2 + 1      # id is non-zero unique values, starting at 1 for rf pulse, 2 for adc period, 3 for next pulse
            block.duration = pulse_length  # pulse length in units of definitions.block_duration_raster_ns
            block.rf = 1                # rf pulse id
            block.gx = 0                # no gradients for now
            block.gy = 0
            block.gz = 0
            block.adc = 0
            block.ext = 0               # extension to save LABEL for counters and flags, TRIGGERS
            yield mrd.StreamItem.Block(block)
            # to next rf pulse
            block = mrd.Block
            block.id += 1               # increment id from the previous value
            block.duration = TR - pulse_length
            block.rf = 0
            block.gx = 0
            block.gy = 0
            block.gz = 0
            block.adc = 0
            block.ext = 0
            yield mrd.StreamItem.Block(block)

            # encode acquisition field
            acq = mrd.Acquisition()
            acq.data = np.transpose(np.expand_dims(g.rawdata[:, 0, 0, 0, 0, iacq], (0)))
            acq.head.acquisition_center_frequency = g.basefreq
            pulse_start = np.uint64(definitions.custom['acq_start_time_ns']) + np.uint64(iacq * TR)
            acq.head.idx.repetition = iacq
            acq.head.acquisition_time_stamp_ns = pulse_start + TE
            acq.head.idx.contrast = 1
            acq.head.discard_pre = 0
            acq.head.discard_post = 0
            acq.head.sample_time_ns = g.sampleperiod * 100                  # convert to ns (sampleperiod is in units of 100ns)
            acq.phase = np.zeros((g.rawdata.shape[0]), dtype=np.float32)    # acquisition phase array set to zeros
            yield mrd.StreamItem.Acquisition(acq)

def groupMRDfiles_collect(rootdir):
    '''
    Find all .MRD files in this directory and subdirectories. Skip files with more than 1 average as they are phantom data.
    Args:
        rootdir: path to the root directory as a string specified by -f command line argument
    Returns
        List of paths of .MRD files as strings with / as path separator
    '''
    l = []
    d = os.listdir(rootdir)
    for f in d:
        if(os.path.isdir(rootdir + '/' + f)):
            l.extend(groupMRDfiles_collect(rootdir + '/' + f))
        elif(f.find('.MRD') > 0):
            mrsdata_filepath = rootdir + '/' + f
            # mrs class is filled when reading navg
            mrs.mread3d(mrsdata_filepath)
            if mrs.navg > 1:
                print(f'Skipping phantom data with {mrs.navg} avg as phantom data', file=sys.stderr)
                continue
            l.append(mrsdata_filepath)
    return(l)

def groupMRDfiles(rootdir, unifylevel):
    '''
    Group .MRD files into groups of files that have the same path up to the unifylevel
    Args:
        rootdir: path to the root directory as a string specified by -f command line argument
        unifylevel: 3 to concatenate each image file into single MRD file
                    1 to keep each image file as a separate MRD file
                    as an integer specified by -u command line argument
    Returns
        List of paths of .MRD files to be grouped together
    '''
    if(not os.path.isdir(rootdir)):
        return([])
    l = groupMRDfiles_collect(rootdir)
    groups = []
    # for each .MRD files in the list,
    for f in l:
        fs = f.split('/')
        addedtogroup = False
        for g in groups:
            gs = g[0].split('/')
            if(len(fs) == len(gs)):
                issame = True
                for i in range(len(fs) - unifylevel):
                    issame = issame and (fs[i] == gs[i])
                if(issame):
                    g.append(f)
                    addedtogroup = True
        if(not addedtogroup):
            groups.append([f])
    return(groups)

def make_header(mrs, measID):
    # make mrd2 header. For now only filling in sequence name but some day should do more
    h = mrd.Header()
    h.measurement_information = mrd.MeasurementInformationType()
    h.measurement_information.sequence_name = mrs.pplfile
    table_pos = mrd.ThreeDimensionalFloat()
    table_pos.x = mrs.FOVoff[0]
    table_pos.y = mrs.FOVoff[1]
    table_pos.z = mrs.FOVoff[2]
    h.measurement_information.relative_table_position = table_pos
    # this may be a misuse of the relative table position but I couldn't find FOV offset anywhere
    # This is definitely a misuse of h1 resonance frequency for non-1H acquisitions
    h.experimental_conditions.h1_resonance_frequency_hz = mrs.basefreq
    h.measurement_information.measurement_id = measID
    et = mrd.EncodingType()
    et.encoded_space.field_of_view_mm.x = mrs.FOV
    et.encoded_space.field_of_view_mm.y = mrs.FOV
    et.encoded_space.field_of_view_mm.z = mrs.FOV
    h.encoding = [et]
    return(h)

# check to see if a directory is specified
def convert_mrs_to_mrd(basedir, unifylevel):
    # assemble groups of MRD files
    groups = groupMRDfiles(basedir, unifylevel)
    # groups are lists of MRS files to be grouped together in a single .MRD2 file format
    for g in groups:
        # get the measurement ID from the first file in the group
        measID = g[0].split('/')[-(unifylevel + 1)]
        basedir = '/'.join(g[0].split('/')[:(-unifylevel)])
        print(f'grouping {len(g)} files into {basedir}', file=sys.stderr)
        w = mrd.BinaryMrdWriter(basedir + '/' + 'raw.mrd2')
        for ig in range(len(g)):
            # write the header once for the entire group after continue
            # re-read the data to set the correct header for each group
            mrs.mread3d(g[ig])
            if ig == 0:
                h = make_header(mrs, measID)
                w.write_header(h)
            w.write_data(generate_pulseq_acquisition(mrs, ig))
        w.close()

# -u 3 consolidates the files as appropriate for EPSI. 
# spectral data like Bukola's and David's uses -u 1
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert MRS data folder to MRD2 format')
    parser.add_argument('-f', '--folder', type=str, required=True,
                        help='Base directory containing MRS data files')
    parser.add_argument('-u', '--unifylevel', type=int, required=False, default=1,
                        help='Directory levels to unify when grouping files (default: 1)')
    args = parser.parse_args()
    print(f'running with base director {args.folder}', file=sys.stderr)
    print(f'unify level set to {args.unifylevel}', file=sys.stderr)
    convert_mrs_to_mrd(args.folder, args.unifylevel)