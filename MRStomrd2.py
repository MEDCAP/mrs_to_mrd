import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from statistics import mode
sys.path.insert(0, 'mrd-fork/python')
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
            a.head.num_echoes = g.nswitch
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
            print(iacq)
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
            a.head.idx.repetition = iacq
            a.head.acquisition_time_stamp_ns = pulse_start + np.uint64(TE * 1.0E+9)
            a.head.num_echoes = 1
            a.head.discard_pre = 0
            a.head.discard_post = 0
            a.head.sample_time_ns = g.sampleperiod * 100  # convert to ns (sampleperiod is in units of 100ns)
            a.phase = np.zeros((g.rawdata.shape[0]), dtype=np.float32)
            yield mrd.StreamItem.Acquisition(a)

def groupMRDfiles_collect(rootdir):
    l = []
    d = os.listdir(rootdir)
    for f in d:
        if(os.path.isdir(rootdir + '/' + f)):
            l.extend(groupMRDfiles_collect(rootdir + '/' + f))
        elif(f.find('.MRD') > 0):
            l.append(rootdir + '/' + f)
    return(l)

def groupMRDfiles(rootdir, unifylevel):
    if(not os.path.isdir(rootdir)):
        return([])
    l = groupMRDfiles_collect(basedir)
    groups = []
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

# MAIN SCRIPT
# check to see if a directory is specified
basedir = '.'
for iarg in range(len(sys.argv)):
    if(sys.argv[iarg] == '-f' and iarg < len(sys.argv) - 1):
        basedir = sys.argv[iarg + 1]
unifylevel = 1
for iarg in range(len(sys.argv)):
    if(sys.argv[iarg] == '-u' and iarg < len(sys.argv) - 1):
        try:
            unifylevel = int(sys.argv[iarg + 1])
        except:
            print('bad unify level set format')
print('running with base director =', basedir)
print('unify level set to', unifylevel)

# assemble groups of MRD files
groups = groupMRDfiles(basedir, unifylevel)
# clean groups so they all have the same length
groups = [[f for f in g if os.path.getsize(f) == mode([os.path.getsize(f) for f in g])] for g in groups]

for g in groups:
    measID = g[0].split('/')[-(unifylevel + 1)]
    basedir = '/'.join(g[0].split('/')[:(-unifylevel)])
    print('grouping', len(g), 'files into', basedir)
    w = mrd.BinaryMrdWriter(basedir + '/' + 'raw.mrd2')
    # make and write output file header
    for ig in range(len(g)):
        mrs.mread3d(g[ig])
        if(ig == 0):
            h = make_header(mrs, measID)
            w.write_header(h)
        w.write_data(generate_acquisition(mrs, ig))
    w.close()
