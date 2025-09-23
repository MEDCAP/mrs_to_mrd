import numpy as np
import matplotlib.pyplot as plt
from statistics import mode
import os
from scipy.optimize import minimize, Bounds
from MRSreader import MRSdata
from lorn import lornfit, lor1fit, lorneval, lor1plot, lornputspect, lornpackx0, lornunpackx0, \
        lorngetpeakparams, lornputpeakparams
import sys
sys.path.insert(0, 'python')
import mrd

debugphasing = False
debuglorn = True


basedir = '.'
fidpad = 4 
lb = 42 # Hz

def findmrd2files(basedir, targetfiletype):
    mrd2list = []
    if(not os.path.isdir(basedir) and basedir.find('raw.mrd2') > -1 and basedir.find(targetfiletype) > -1):
        return([basedir])
    if(os.path.isdir(basedir)):
        d = os.listdir(basedir)
        for f in d:
            mrd2list.extend(findmrd2files(basedir + '/' + f, targetfiletype))
    return(mrd2list)

def generate_aux_images(imglist):
    for iimg in range(len(imglist)):
        imghead = mrd.ImageHeader(image_type=mrd.ImageType.BITMAP)
        img = mrd.Image(head=imghead, data=np.expand_dims(imglist[iimg], (2, 3, 4)))
        yield(mrd.StreamItem.ImageUint32(img))

def get_clarg(clarg, arg, ty):
    clargarr = clarg.split(' ')
    for iarg in range(len(clargarr) - 1):
        print(iarg, clargarr[iarg])
        if(clargarr[iarg] == arg):
            return(ty(clargarr[iarg + 1]))

def generate_epsi_images(h, m):
    # turn the metabolite array (metabolite x image number x rows x columns) into streamable images
    time_between_images = 3   # approximately 3s between images
    nmet = m.shape[0]
    nimg = m.shape[1]
    measfreq = h.experimental_conditions.h1resonance_frequency_hz
    fov = np.array([h.encoding[0].encoded_space.field_of_view_mm.x, \
            h.encoding[0].encoded_space.field_of_view_mm.y, \
            h.encoding[0].encoded_space.field_of_view_mm.z], dtype=np.float32)
    pos = h.measurement_information.relative_table_position
    for ide in range(nimg):
        imghead = mrd.ImageHeader(image_type=mrd.ImageType.MAGNITUDE)
        if(ide == 0):
            imghead.flags = mrd.ImageFlags.FIRST_IN_SET
        elif(ide == nimg - 1):
            imghead.flags = mrd.ImageFlags.LAST_IN_SET
        imghead.measurement_uid = ide
        imghead.measurement_freq = measfreq + np.uint32(measfreq * peakoffsets / 1E+6 + 0.5)
        measfreq = h.experimental_conditions.h1resonance_frequency_hz
        imghead.measurement_freq_label = np.array(peaknames, dtype=np.dtype(np.object_))
        imghead.field_of_view = fov
        imghead.position = np.array([pos.x, pos.y, pos.z], dtype=np.float32)
        imghead.col_dir = np.zeros((3,), dtype=np.dtype(np.float32))
        imghead.line_dir = np.zeros((3,), dtype=np.dtype(np.float32))
        imghead.slice_dir = np.zeros((3,), dtype=np.dtype(np.float32))
        imghead.patient_table_position = np.zeros((3,), dtype=np.dtype(np.float32))
        imghead.average = 0
        imghead.slice = 0    # only one slice
        imghead.contrast = 0 # not sure what to do with this yet
        imghead.phase = 0    # or this
        imghead.repetition = ide
        imghead.set = 0      # or this
        imghead.acquisition_time_stamp_ns = ide * time_between_images * 1000000000
        imghead.physiology_time_stamp_ns = []
        imghead.image_type = mrd.ImageType.MAGNITUDE | mrd.ImageType.SPIN_DENSITY_MAP
        imghead.image_index = ide
        imghead.image_series_index = ide
        imghead.user_int = []
        imghead.user_float = []
        img = mrd.Image(head=imghead, data=np.expand_dims(np.transpose(m[:, ide, :, :]), (0, 1)))
        yield(mrd.StreamItem.ImageDouble(img))

def epsi_recon():
    auximages = []
    numimages = sum([a.head.flags & mrd.AcquisitionFlags.LAST_IN_PHASE == \
            mrd.AcquisitionFlags.LAST_IN_PHASE for a in raw_acquisition_list])
    a = raw_acquisition_list[0]
    sampletime = a.head.sample_time_ns / 1.0E+9
    centerfreq = a.head.acquisition_center_frequency
    totalppswitch = int(a.data.shape[0] / a.head.num_echoes + 0.1)
    nro = totalppswitch - a.head.discard_pre - a.head.discard_post
    npe = int(len(raw_acquisition_list) / numimages + 0.1)
    kspace = np.zeros((numimages, npe, nro, a.head.num_echoes * fidpad), dtype = 'complex')
    ia = 0
    for iimg in range(kspace.shape[0]):
        for ipe in range(kspace.shape[1]):
            a = raw_acquisition_list[ia]
            for iecho in range(a.head.num_echoes):
                kspace[iimg, ipe, :, iecho] = a.data[(iecho * totalppswitch + \
                    a.head.discard_pre):(iecho * totalppswitch + a.head.discard_pre + kspace.shape[2]), 0]
            ia += 1
    print('doing fft')
    img = np.fft.fftshift(np.fft.fftn(kspace, axes = (1, 2, 3)), axes = (1, 2, 3))
    print('finding biggest voxel')
    currmax = 0
    for ide in range(numimages):
        for j in range(npe):
            for k in range(nro):
                thismax = np.max(np.abs(img[ide, j, k, :]))
                if(thismax > currmax):
                    currmax = thismax
                    maxj = j
                    maxk = k
                    maxide = ide
    # estimate noise level by looking at the last image in the series
    noise = np.mean(np.abs(img[-1, :, :, :]))
    print('noise =', noise)
    maxspect = img[maxide,maxj,maxk,:].copy()
    globalspect = np.zeros((a.head.num_echoes * fidpad), dtype = 'complex')
    for ide in range(numimages):
        print('phasing img', ide)
        # go through all the voxels and minimize the phase difference
        for j in range(npe):
            for k in range(nro):
                thisspect = img[ide, j, k, :]
                if(np.max(np.abs(thisspect)) < noise * 3):
                    continue
                bestoverlap = 0
                for r in range(-15,16):
                    thisrollspect = np.roll(thisspect,r)
                    S0 = np.sum(np.real(thisrollspect * np.conj(maxspect)))
                    Spi2 = np.sum(np.real(thisrollspect * 1j * np.conj(maxspect)))
                    overlap = S0**2 + Spi2**2
                    if(overlap > bestoverlap):
                        bestr = r
                        bestoverlap = overlap
                        th0 = np.pi / 2 - np.arctan2(S0, Spi2)
                img[ide,j,k,:] = np.roll(img[ide,j,k,:], bestr) * np.exp(1j * th0)
                if(debugphasing):
                    plt.clf()
                    plt.plot(np.real(img[ide,j,k,:]))
                    plt.plot(np.real(maxspect))
                    plt.plot([0, 256], [100*bestr, 100*bestr], 'k')
                    plt.draw()
                    plt.pause(.1)
                if(ide > 1):
                    globalspect += img[ide, j, k, :]
    print('end phasing')
    BW = 1 / sampletime / totalppswitch
    xscale = np.array(range(len(globalspect))) / len(globalspect) * BW / centerfreq * 1E+6
    globalspect /= np.max(np.abs(globalspect))
    # estimate peak widths using the FWHM of the largest peak
    maxpeakidx = np.argmax(np.abs(globalspect))
    leftidx = -1
    rightidx = -1
    for isp in range(len(globalspect)):
        if(np.abs(globalspect[(maxpeakidx - isp) % len(globalspect)]) < 0.5 and leftidx == -1):
            leftidx = -isp
        if(np.abs(globalspect[(maxpeakidx + isp) % len(globalspect)]) < 0.5 and rightidx == -1):
            rightidx = isp
    widthguess = (rightidx - leftidx) * (xscale[1] - xscale[0]) / 2
    lornputspect(xscale, globalspect, widthguess, 1.0, debuglorn)
    # fit global spectrum to 6 lorentzians. Biggest peak has got to be either pyruvate or urea. 
    # Maybe some day this will fail if it's lactate
    x0 = np.zeros(((4 * npeaks) + 2))
    x1 = np.zeros((len(biggestpeaklist), len(x0)))
    diff = np.zeros((len(biggestpeaklist)))
    for icg in range(len(biggestpeaklist)):
        centers = (xscale[np.argmax(np.abs(globalspect))] - (peakoffsets - \
                peakoffsets[biggestpeaklist[icg]])) % (BW / centerfreq * 1E+6)
        for ip in range(npeaks):
            x0[3 * npeaks + ip] = np.abs(globalspect[np.argmin(np.abs(xscale - centers[ip]))])
            x0[2 * npeaks + ip] = np.angle(globalspect[np.argmin(np.abs(xscale - centers[ip]))])
        lornputpeakparams(centers, np.ones((npeaks)) * widthguess, x0[(2 * npeaks):(3 * npeaks)], debuglorn)
        print('begin minimize', icg)
        x1[icg, :] = minimize(lornfit, x0).x
        for ip in range(npeaks):
            if(x1[icg, 3 * npeaks + ip] < 0):
                x1[icg, 3 * npeaks + ip] *= -1
                x1[icg, 2 * npeaks + ip] += np.pi
        diff[icg] = np.sum(np.abs(globalspect - lorneval(x1[icg, :])))
#        plt.plot(xscale, np.real(globalspect), 'r')
#        plt.plot(xscale, np.imag(globalspect), 'g')
#        plt.plot(xscale, np.real(lorneval(x1[icg, :])), 'k')
#        plt.plot(xscale, np.imag(lorneval(x1[icg, :])), 'k')
#        plt.draw()
#        plt.pause(.1)
        print("icg ", icg, "diff = ", diff[icg])
    centers = (xscale[np.argmax(np.abs(globalspect))] - (peakoffsets - \
            peakoffsets[biggestpeaklist[np.argmin(diff)]])) % (BW / centerfreq * 1E+6)
    lornputpeakparams(centers, np.ones((npeaks)) * widthguess, x0[(2 * npeaks):(3 * npeaks)], debuglorn)
    centers, widths, phases, amplitudes, baseline = lornunpackx0(x1[np.argmin(diff)], debuglorn)
    specteval = lorneval(x1[np.argmin(diff)])
    plt.clf()
    plt.plot(xscale, np.real(globalspect), 'r')
    plt.plot(xscale, np.imag(globalspect), 'g')
    plt.plot(xscale, np.real(specteval), 'k')
    plt.plot(xscale, np.imag(specteval), 'k')
    plt.gcf().canvas.draw()
    thisbmp = np.array(plt.gcf().canvas.renderer._renderer, dtype=np.uint32)
    thisbmp[:,:,0] = thisbmp[:,:,0]*2**12 + thisbmp[:,:,1]*2**8 + thisbmp[:,:,2]*2**4 + thisbmp[:,:,3]
    auximages.append(thisbmp[:,:,0])
    for ip in range(0, npeaks):
        plt.plot([centers[ip], centers[ip]], [-1, 1], 'k')
        plt.text(centers[ip], .95-ip*.07, str(centers[ip]))
#    plt.savefig(basedir + '/' + fn + '/epsi/globalfit', dpi=300, bbox_inches='tight')
#    np.save(basedir + '/' + fn + '/epsi/xscale', xscale)
#    np.save(basedir + '/' + fn + '/epsi/globalfit_centers', centers)
#    np.save(basedir + '/' + fn + '/epsi/globalfit_widths', widths)
#    np.save(basedir + '/' + fn + '/epsi/globalfit_amplitudes', amplitudes)
#    np.save(basedir + '/' + fn + '/epsi/globalfit_phases', phases)
#    np.save(basedir + '/' + fn + '/epsi/gloablfit_baseline', baseline)
#    np.save(basedir + '/' + fn + '/epsi/globalfit_globalspect', globalspect)
    # now do voxel fits
    lornputpeakparams(centers, widths, phases, debuglorn)
    metabolites = np.zeros((npeaks, numimages, npe, nro))
    for ide in range(numimages):
        print('voxel fit img', ide)
        for j in range(npe):
            for k in range(nro):
                thisspect = img[ide,j,k,:]
                scaling = np.max(np.abs(thisspect))
                if(np.max(np.abs(thisspect)) < noise * 3):
                    continue
                thisspect /= scaling
                lornputspect(xscale, thisspect, widths, 1.0, False)
                x0 = np.zeros((npeaks + 2))
                for ip in range(npeaks):
                    x0[ip] = np.abs(thisspect[np.argmin(np.abs(xscale - centers[ip]))])
                bounds = Bounds(np.concatenate((np.zeros((npeaks)), [-.1, -.1])), \
                        np.concatenate((x0[:npeaks] * 1.5, [.1, .1])))
                x1 = minimize(lor1fit, x0, bounds=bounds)
                metabolites[:, ide, j, k] = x1.x[:npeaks] * scaling
    return([metabolites, auximages])

def generate_spectra(h, measurementtimes_ns, peakfrequencies, peakamplitudes, spectra):
    # turn the metabolite array (metabolite x image number x rows x columns) into streamable images
    nspect = spectra.shape[0]
    ntimepoints = spectra.shape[1]
    measfreq = h.experimental_conditions.h1resonance_frequency_hz
    for ispect in range(nspect):
        # 'image' that is the measured spectrum at this time point
        imghead = mrd.ImageHeader(image_type=mrd.ImageType.COMPLEX)
        if(ispect == 0):
            imghead.flags = mrd.ImageFlags.FIRST_IN_SET
        elif(ispect == nspect - 1):
            imghead.flags = mrd.ImageFlags.LAST_IN_SET
        imghead.measurement_uid = ispect
        imghead.measurement_freq = measfreq
        imghead.repetition = ispect
        imghead.acquisition_time_stamp_ns = measurementtimes_ns[ispect]
        imghead.image_index = ispect
        imghead.image_series_index = ispect
        spect = mrd.Image(head=imghead, data=np.expand_dims(np.transpose(spectra[ispect, :]), (0, 1, 2, 3)))
        yield(mrd.StreamItem.ImageComplexDouble(spect))
        # 'image' quantifying the amplitudes for this time point
        imghead = mrd.ImageHeader(image_type=mrd.ImageType.MAGNITUDE)
        if(ispect == 0):
            imghead.flags = mrd.ImageFlags.FIRST_IN_SET
        elif(ispect == ntimepoints - 1):
            imghead.flags = mrd.ImageFlags.LAST_IN_SET
        imghead.measurement_uid = ispect
        imghead.measurement_freq = peakfrequencies
        imghead.measurement_freq_label = np.array(peaknames, dtype=np.dtype(np.object_))
        imghead.repetition = ispect
        imghead.acquisition_time_stamp_ns = measurementtimes_ns[ispect]
        imghead.image_index = ispect
        imghead.image_series_index = ispect
        spect = mrd.Image(head=imghead, data=np.expand_dims(np.transpose(peakamplitudes[:, ispect]), (0, 1, 2, 3)))
        yield(mrd.StreamItem.ImageDouble(spect))

def spectra_recon(h):
    auximages = []
    numspectra = len(raw_acquisition_list)
    a = raw_acquisition_list[0]
    sampletime = a.head.sample_time_ns / 1.0E+9
    centerfreq = a.head.acquisition_center_frequency
    sampletime = a.head.sample_time_ns / 1.0E+9
    npts = len(a.data)
    kspace = np.zeros((numspectra, len(a.data)), dtype = 'complex')
    ia = 0
    for ispect in range(kspace.shape[0]):
        kspace[ispect, :] = raw_acquisition_list[ispect].data[:, 0]
        for ipt in range(kspace.shape[1]):
             tk = ipt * sampletime
             kspace[ispect, ipt] *= np.exp(-tk * lb)
    spectra = np.fft.fftshift(np.fft.fft(kspace, axis = (1)), axes = (1))
    currmax = 0
    for ispect in range(numspectra):
        thismax = np.max(np.abs(spectra[ispect, :]))
        if(thismax > currmax):
            currmax = thismax
            maxispect = ispect
    # estimate noise level by looking at the last image in the series
    noise = np.mean(np.abs(spectra[-1, :]))
    print('noise =', noise)
    maxspect = spectra[maxispect,:].copy()
    maxpt = np.argmax(np.abs(maxspect))
    BW = 1 / sampletime
    # take points within 25 ppm of highest signal
    lowpt = int(maxpt - 25E-6 / (BW / centerfreq) * npts)
    hipt = int(maxpt + 25E-6 / (BW / centerfreq) * npts)
    newBW = BW * (hipt - lowpt) / npts
    spectra = spectra[:, lowpt:hipt]
    globalspect = np.zeros(hipt - lowpt, dtype = 'complex')
    xscale = np.array(range(len(globalspect))) / len(globalspect) * newBW / centerfreq * 1E+6
    for ispect in range(numspectra):
        if(np.max(np.abs(spectra[ispect, :])) > noise * 5):
            globalspect += spectra[ispect, :]
    globalspect /= np.max(np.abs(globalspect))
    # estimate peak widths using the FWHM of the largest peak
    maxpeakidx = np.argmax(np.abs(globalspect))
    leftidx = -1
    rightidx = -1
    for isp in range(len(globalspect)):
        if(np.abs(globalspect[(maxpeakidx - isp) % len(globalspect)]) < 0.5 and leftidx == -1):
             leftidx = maxpeakidx - isp
        if(np.abs(globalspect[(maxpeakidx + isp) % len(globalspect)]) < 0.5 and rightidx == -1):
             rightidx = maxpeakidx + isp
    widthguess = (rightidx - leftidx) * (xscale[1] - xscale[0]) / 4
    lornputspect(xscale, globalspect, widthguess, wigglefactor, debuglorn)
    # fit global spectrum to npeaks lorentzians.
    x0 = np.zeros(((4 * npeaks) + 2))
    x1 = np.zeros((len(biggestpeaklist), len(x0)))
    diff = np.zeros((len(biggestpeaklist)))
    for icg in range(len(biggestpeaklist)):
        centers = (xscale[np.argmax(np.abs(globalspect))] - (peakoffsets - \
                peakoffsets[biggestpeaklist[icg]])) % (BW / centerfreq * 1E+6)
        for ip in range(npeaks):
             x0[3 * npeaks + ip] = np.abs(globalspect[np.argmin(np.abs(xscale - centers[ip]))])
             x0[2 * npeaks + ip] = np.angle(globalspect[np.argmin(np.abs(xscale - centers[ip]))])
        lornputpeakparams(centers, np.ones((npeaks)) * widthguess, x0[(2 * npeaks):(3 * npeaks)], debuglorn)
        print('begin minimize', icg)
        x1[icg, :] = minimize(lornfit, x0).x
        for ip in range(npeaks):
             if(x1[icg, 3 * npeaks + ip] < 0):
                 x1[icg, 3 * npeaks + ip] *= -1
                 x1[icg, 2 * npeaks + ip] += np.pi
        diff[icg] = np.sum(np.abs(globalspect - lorneval(x1[icg, :])))
    centers = (xscale[np.argmax(np.abs(globalspect))] - (peakoffsets - \
                peakoffsets[biggestpeaklist[np.argmin(diff)]])) % (BW / centerfreq * 1E+6)
    lornputpeakparams(centers, np.ones((npeaks)) * widthguess, x0[(2 * npeaks):(3 * npeaks)], debuglorn)
    thex1 = x1[np.argmin(diff)]
    centers, widths, phases, amplitudes, baseline = lornunpackx0(thex1, debuglorn)
    specteval = lorneval(thex1)
    plt.clf()
    plt.plot(xscale, np.real(globalspect), 'r')
    plt.plot(xscale, np.imag(globalspect), 'g')
    plt.plot(xscale, np.real(specteval), 'k')
    plt.plot(xscale, np.imag(specteval), 'k')
    for ip in range(0, npeaks):
        plt.plot([centers[ip], centers[ip]], [-1, 1], 'k')
        plt.text(centers[ip], .95-ip*.07, str(centers[ip]))
    plt.gcf().canvas.draw()
    thisbmp = np.array(plt.gcf().canvas.renderer._renderer, dtype=np.uint32)
    thisbmp[:,:,0] = thisbmp[:,:,0]*2**12 + thisbmp[:,:,1]*2**8 + thisbmp[:,:,2]*2**4 + thisbmp[:,:,3]
    auximages.append(thisbmp[:,:,0])
    # now do voxel fits
    lornputpeakparams(centers, widths, phases, debuglorn)
    peakamplitudes = np.zeros((npeaks, numspectra))
    measurementtimes_ns = [a.head.acquisition_time_stamp_ns - \
            raw_acquisition_list[0].head.acquisition_time_stamp_ns for a in raw_acquisition_list]
    for ispect in range(numspectra):
        print('voxel fit img', ispect)
        thisspect = spectra[ispect,:]
        scaling = np.max(np.abs(thisspect))
        thisspect /= scaling
        lornputspect(xscale, thisspect, widths, wigglefactor, False)
        x0 = np.zeros((npeaks + 2))
        for ip in range(npeaks):
            x0[ip] = np.abs(thisspect[np.argmin(np.abs(xscale - centers[ip]))])
        bounds = Bounds(np.concatenate((np.zeros((npeaks)), [-.1, -.1])), np.concatenate((x0[:npeaks] * 1.5, [.1, .1])))
        x1 = minimize(lor1fit, x0, bounds=bounds)
        print(x1.x)
        peakamplitudes[:, ispect] = x1.x[:npeaks] * scaling
    plt.clf()
    for ipeak in range(npeaks):
        plt.plot(np.array(measurementtimes_ns) * 1.0E-9, peakamplitudes[ipeak,:])
    plt.gcf().canvas.draw()
    thisbmp = np.array(plt.gcf().canvas.renderer._renderer, dtype=np.uint32)
    thisbmp[:,:,0] = thisbmp[:,:,0]*2**12 + thisbmp[:,:,1]*2**8 + thisbmp[:,:,2]*2**4 + thisbmp[:,:,3]
    auximages.append(thisbmp[:,:,0])
    return(measurementtimes_ns, spectra, centerfreq + np.uint32(centers * centerfreq / 1.0E+6), \
            peakamplitudes, auximages)

# read command line arguments
wigglefactor = 1.0
peakoffsets = []
peaknames = []
biggestpeaklist = []
targetfiletype = ''
for iarg in range(len(sys.argv) - 1):
    try:
        floatarg = float(sys.argv[iarg + 1])
    except:
        floatarg = np.nan
    if(sys.argv[iarg] == '-f'):
        basedir = sys.argv[iarg + 1]
        print('setting base dir to', basedir)
        continue
    if(sys.argv[iarg] == '-w' and not np.isnan(floatarg)):
        wigglefactor = floatarg
        print('setting wiggle factor to', wigglefactor)
        continue
    if(sys.argv[iarg] == '-n'):
        targetfiletype = sys.argv[iarg + 1]
        print('setting target file type to', targetfiletype)
        continue
    if(sys.argv[iarg][0] == '-' and not np.isnan(floatarg)):
        if(sys.argv[iarg][-1] == '*'):
            # this peak is certified small, do not include it in the 'biggest peak' list
            peaknames.append(sys.argv[iarg][1:-1])
        else:
            peaknames.append(sys.argv[iarg][1:])
            biggestpeaklist.append(len(peaknames) - 1)
        peakoffsets.append(floatarg)
fnames = findmrd2files(basedir, targetfiletype)
# now look for specification of metabolite peaks
# BA's cirrhrat is -bic* 0.0 -urea 2.3 -pyr 9.7 -ala* 15.2 -hyd* 18.1 -lac 21.8
# SZ's mouse kidney is -bic 0.0 -urea 2.3 -pyr 9.7 -ala 15.2 -poop 15.9 -hyd 18.1 -lac 21.8
# BA's spectra is -bic* -0.4 -urea 2.1 -urea2* 2.3 -pyr 9.7 -ala* 15.2 -hyd* 18.1 -lac 21.8 -w 0.5
# DT's spectra is -urea 0.0 -KIC 8.6 -leu* 13.0 -hyd* 18.1 -?* 21.8 -w 1.5
peakoffsets = np.array(peakoffsets)
npeaks = len(peakoffsets)

# read raw data mrd file
for f in fnames:
    print('doing', f)
    with mrd.BinaryMrdReader(f) as w:
        raw_header = w.read_header()
        raw_streamables_list = list(w.read_data())
        raw_pulse_list = [x.value for x in raw_streamables_list if type(x.value) == mrd.Pulse]
        raw_pulse_list.sort(key = lambda x: x.head.pulse_time_stamp_ns)
        raw_gradient_list = [x.value for x in raw_streamables_list if type(x.value) == mrd.Gradient]
        raw_gradient_list.sort(key = lambda x: x.head.gradient_time_stamp_ns)
        raw_acquisition_list = [x.value for x in raw_streamables_list if type(x.value) == mrd.Acquisition]
        raw_acquisition_list.sort(key = lambda x: x.head.acquisition_time_stamp_ns)
        img_list = [x.value for x in raw_streamables_list if type(x.value) == mrd.Image]

    with mrd.BinaryMrdWriter(f.replace('raw.mrd2', 'recon.mrd2')) as w:
        # write output file header
        w.write_header(raw_header)
        if(raw_header.measurement_information.sequence_name.find('epsi') > -1):
            [metabolites, auximages] = epsi_recon()
            w.write_data(generate_epsi_images(raw_header, metabolites))
        elif(raw_header.measurement_information.sequence_name.find('1pul') > -1):
            [measurementtimes_ns, spectra, peakfrequencies, peakamplitudes, auximages] = spectra_recon(raw_header)
            w.write_data(generate_spectra(raw_header, measurementtimes_ns, peakfrequencies, peakamplitudes, spectra))
        if(len(auximages) > 0):
            w.write_data(generate_aux_images(auximages))
