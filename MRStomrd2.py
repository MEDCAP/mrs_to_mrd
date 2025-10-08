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
from pathlib import Path
sys.path.append(os.path.join(Path(__file__).parent, 'mrd-fork','python'))
import mrd
from MRSreader import MRSdata

mrs = MRSdata()

def generate_acquisition(g, ide):
    # make the pulse, gradient, and (raw data) acquisitions. g is the MRSdata structure read in from an individual
    # file, and ide is the index of this file in the image series
    npts = g.rawdata.shape[0] # number of points per acquisition
    pulse_len = 100E-6        # guess 100us, I don't know what it is
    pulsedt = 10.0E-6         # specify pulse in steps of 10us
    TE = 180E-6               # just an estimate for now, start acquiring 180us after beginning of 100us pulse
    if(g.pplfile.find('epsi') >= 0):
        nPE = g.rawdata.shape[1]  # number of acquisitions per image file (phase-encodes)
        for iacq in range(nPE):
            # make the pulse
            p = mrd.Pulse()
            CO_freq = g.basefreq
            pulse_start = np.uint64(g.acqstarttime * 100) + np.uint64(iacq * g.tr * 1.0E+6)
            pulse_end = pulse_start + np.uint64(pulse_len * 1.0E+9)
            p.head.pulse_time_stamp_ns = pulse_start
            p.head.sample_time_ns = g.sampleperiod * 100  # convert to ns (sampleperiod is in units of 100ns)
            npulsepoints = int(pulse_len / pulsedt)
            p.amplitude = np.zeros((1, npulsepoints)).astype(np.float32)
            p.phase = np.zeros(npulsepoints).astype(np.float32)
            for idt in range(npulsepoints):
                p.amplitude[0, idt] = 1.0    # we don't know the pulse amplitude or phase yet
                p.phase[idt] = 0.0
            yield mrd.StreamItem.Pulse(p)
            # make the gradients
            gr = mrd.Gradient()
            gr.head.gradient_time_stamp_ns = pulse_start + np.uint64(TE * 1.0E+9)
            gr.head.gradient_sample_time_ns = 10000 # 10us
            # at some point we'll have to reproduce the whole flyback gradient waveform but not now
            yield mrd.StreamItem.Gradient(gr)
            # make the acquisition
            a = mrd.Acquisition()
            if(iacq == 0):
                a.head.flags = mrd.AcquisitionFlags.FIRST_IN_PHASE
            if(iacq == g.rawdata.shape[1] - 1):
                a.head.flags = mrd.AcquisitionFlags.LAST_IN_PHASE
            a.data = np.transpose(np.expand_dims(g.rawdata[:, iacq, 0, 0, 0, 0], (0)))
            a.head.acquisition_center_frequency = g.basefreq
            a.head.idx.phase = iacq
            a.head.acquisition_time_stamp_ns = pulse_start + np.uint64(TE * 1.0E+9)
            a.head.idx.contrast = g.nswitch
            totalppswitch = int(g.rawdata.shape[0] / g.nswitch + 0.7)
            a.head.discard_pre = int((totalppswitch - g.nppswitch) / 2)
            a.head.discard_post = a.head.discard_pre
            a.head.sample_time_ns = g.sampleperiod * 100  # convert to ns (sampleperiod is in units of 100ns)
            a.phase = np.zeros((g.rawdata.shape[0]), dtype=np.float32)
            yield mrd.StreamItem.Acquisition(a)
    elif(g.pplfile.find('1pul') >= 0):
        nrep = g.rawdata.shape[5]  # number of acquisitions per image file (phase-encodes)
        print(nrep)
        for iacq in range(nrep):
            print('converting measurement #: ', iacq)
            # make the pulse
            p = mrd.Pulse()
            CO_freq = g.basefreq
            pulse_start = np.uint64(g.acqstarttime * 100) + np.uint64(iacq * g.tr * 1.0E+6)
            pulse_end = pulse_start + np.uint64(pulse_len * 1.0E+9)
            p.head.pulse_time_stamp_ns = pulse_start
            p.head.sample_time_ns = g.sampleperiod * 100  # convert to ns (sampleperiod is in units of 100ns)
            npulsepoints = int(pulse_len / pulsedt)
            p.amplitude = np.zeros((1, npulsepoints)).astype(np.float32)
            p.phase = np.zeros(npulsepoints).astype(np.float32)
            for idt in range(npulsepoints):
                p.amplitude[0, idt] = 1.0    # we don't know the pulse amplitude or phase yet
                p.phase[idt] = 0.0
            yield mrd.StreamItem.Pulse(p)
            # no gradients
            # make the acquisition
            a = mrd.Acquisition()
            a.data = np.transpose(np.expand_dims(g.rawdata[:, 0, 0, 0, 0, iacq], (0)))
            a.head.acquisition_center_frequency = g.basefreq
            # print('writing center freq = ', a.head.acquisition_center_frequency)
            a.head.idx.repetition = iacq
            a.head.acquisition_time_stamp_ns = pulse_start + np.uint64(TE * 1.0E+9)
            a.head.idx.contrast = 1
            a.head.discard_pre = 0
            a.head.discard_post = 0
            a.head.sample_time_ns = g.sampleperiod * 100  # convert to ns (sampleperiod is in units of 100ns)
            a.phase = np.zeros((g.rawdata.shape[0]), dtype=np.float32)
            yield mrd.StreamItem.Acquisition(a)

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
            mrs.mread3d(mrsdata_filepath)
            if mrs.navg > 1:
                print(f'Skipping {mrsdata_filepath} with {mrs.navg} avg as phantom data')
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
    h.measurement_information.relative_table_position = mrd.ThreeDimensionalFloat()
    # this may be a misuse of the relative table position but I couldn't find FOV offset anywhere
    h.measurement_information.relative_table_position.x = mrs.FOVoff[0]
    h.measurement_information.relative_table_position.y = mrs.FOVoff[1]
    h.measurement_information.relative_table_position.z = mrs.FOVoff[2]
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
        print('grouping', len(g), 'files into', basedir)
        w = mrd.BinaryMrdWriter(basedir + '/' + 'raw.mrd2')
        for ig in range(len(g)):
            # write the header once for the entire group after continue
            if ig == 0:
                h = make_header(mrs, measID)
                w.write_header(h)
            w.write_data(generate_acquisition(mrs, ig))
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
    print('running with base director =', args.folder)
    print('unify level set to', args.unifylevel)
    convert_mrs_to_mrd(args.folder, args.unifylevel)