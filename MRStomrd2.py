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
import mrd
from MRSreader import MRSdata

mrs = MRSdata()

def generate_pulseq_acquisition(mrs, ide):
    '''
    Define pulseq fields and extract acquisition from MRSdata structure g 
    '''
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
    if(mrs.pplfile.find('epsi') >= 0):
        import pdb; pdb.set_trace()
        [npts, nPE] = mrs.rawdata.shape   # number of points per acquisition and phase-encoding
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
    elif(mrs.pplfile.find('1pul') >= 0):
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
    meas = mrd.MeasurementInformationType()
    meas.sequence_name = mrs.pplfile
    # this may be a misuse of the relative table position but I couldn't find FOV offset anywhere
    meas.relative_table_position = mrd.ThreeDimensionalFloat(x=mrs.FOVoff[0], y=mrs.FOVoff[1], z=mrs.FOVoff[2])
    meas.measurement_id = measID
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