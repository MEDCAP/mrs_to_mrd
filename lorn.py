import numpy as np
import matplotlib.pyplot as plt

centers = []
widths = []
phases = []
spect = []
xscale = []
BW = 0
wigglefactor = 1

def lornunpackx0(x0, debug):
    # takes an array of 4 x number of peaks + 2 (complex offset) numbers and unpacks them into 
    npeaks = int(len(x0) / 4)

    #c = centers + widths * np.arctan(x0[:npeaks]) / np.pi
    #w = widths * (1 + np.arctan(x0[(npeaks):(2 * npeaks)]) * 2 / np.pi)
    # changed from above to below for Bukola's bicarbonate
    c = centers + np.arctan(x0[:npeaks]) / np.pi * wigglefactor
    w = widths * (1 + np.arctan(x0[(npeaks):(2 * npeaks)]) * 1.8 / np.pi)

    ph = x0[(2 * npeaks):(3 * npeaks)]
    A = x0[(3 * npeaks):(4 * npeaks)]
    b = x0[4 * npeaks] + 1j * x0[4 * npeaks + 1]
    if(debug):
        print('unpacking x0 =', x0)
        print('  --> centers =', c)
        print('  --> widths =', w)
        print('  --> phases =', ph)
        print('  --> amplitudes =', A)
        print('  --> baseline =', b)
    return(c, w, ph, A, b)

def lornpackx0(c, w, ph, a, b, debug):
    npeaks = len(c)
    x0 = np.zeros((4 * len(c) + 2))

    #x0[:npeaks] = np.tan((c - centers) * np.pi / widths)
    #x0[npeaks:(2 * npeaks)] = np.tan((w / widths - 1) * np.pi / 2)
    # changed from above to below for Bukola's bicarbonate
    x0[:npeaks] = np.tan((c - centers) / wigglefactor * np.pi / 1.8)
    x0[npeaks:(2 * npeaks)] = np.tan((w / widths - 1) * np.pi)

    x0[(2 * npeaks):(3 * npeaks)] = ph
    x0[(3 * npeaks):(4 * npeaks)] = a
    x0[4 * npeaks] = np.real(b)
    x0[4 * npeaks + 1] = np.imag(b)
    if(debug):
        print('packing centers =', c)
        print('  widths =', w)
        print('  phases =', ph)
        print('  amplitudes =', a)
        print('  baseline =', b)
        print('  --> x0 =', x0)
    return(x0)

def lornputspect(x, g, w, wf, debug):
    global spect
    global xscale
    global BW
    global widths
    if(debug):
        print('putting spect')
    xscale = x
    BW = max(xscale) - min(xscale) + (xscale[1] - xscale[0])
    spect = g
    widths = w
    wigglefactor = wf

def lorngetpeakparams(debug):
    global centers, widths, phases
    if(debug):
        print('getting centers =', centers)
        print('getting widths =', widths)
        print('getting phases =', phases)
    return(centers, widths, phases)

def lornputpeakparams(c, w, ph, debug):
    global centers, widths, phases
    if(debug):
        print('putting centers =', c)
        print('putting widths =', w)
        print('putting phases =', ph)
    centers = c
    widths = w
    phases = ph

def lorneval(x0):
    global spect
    c, w, ph, A, b = lornunpackx0(x0, False)
    y = np.zeros(len(spect), dtype='complex') + b

    for j in range(len(c)):
        y += A[j] * np.exp(1j * ph[j]) / (1 + 1j * (xscale - c[j]) / w[j])
        y += A[j] * np.exp(1j * ph[j]) / (1 + 1j * (xscale - c[j] - BW) / w[j])
        y += A[j] * np.exp(1j * ph[j]) / (1 + 1j * (xscale - c[j] + BW) / w[j])

#    plt.clf()
#    plt.plot(np.real(spect), 'r')
#    plt.plot(np.imag(spect), 'g')
#    plt.plot(np.real(y), 'b')
#    plt.plot(np.imag(y), 'c')
#    plt.draw()
#    plt.pause(.001)
    return(y)

def lornfit(x0):
    global spect
    y = lorneval(x0)
    return(np.sum(abs(y - spect)))


def lor1fit(x0):
    global spect, centers, widths, phases, xscale
    amplitudes = x0[:-2]
    baseline = x0[-2] + 1j * x0[-1]
    y = np.zeros(len(spect), dtype='complex') + baseline
    diff = 0

    for ip in range(len(centers)):
        c = amplitudes[ip] * np.exp(1j * phases[ip])
        y += c / (1 + 1j * (xscale - centers[ip]) / widths[ip])
        y += c / (1 + 1j * (xscale - centers[ip] - BW) / widths[ip])
        y += c / (1 + 1j * (xscale - centers[ip] + BW) / widths[ip])
    diff = np.sum(abs(y-spect))

#    for isp in range(len(spect)):
#        for ip in range(len(centers)):
#            c = amplitudes[ip] * np.exp(1j * phases[ip])
#            y[isp] += c / (1 + 1j * (xscale[isp] - centers[ip]) / widths[ip])
#            y[isp] += c / (1 + 1j * (xscale[isp] - centers[ip] - BW) / widths[ip])
#            y[isp] += c / (1 + 1j * (xscale[isp] - centers[ip] + BW) / widths[ip])
#            diff += abs(y[isp]-spect[isp])
    return(diff)

def lor1plot(x0):
    global spect, centers, widths, phases, xscale
    amplitudes = x0[:-2]
    baseline = x0[-2] + 1j * x0[-1]
    y = np.zeros(len(spect), dtype='complex') + baseline
    for isp in range(len(spect)):
        for ip in range(len(centers)):
            c = amplitudes[ip] * np.exp(1j * phases[ip])
            y[isp] += c / (1 + 1j * (xscale[isp] - centers[ip]) / widths[ip])
            y[isp] += c / (1 + 1j * (xscale[isp] - centers[ip] - BW) / widths[ip])
            y[isp] += c / (1 + 1j * (xscale[isp] - centers[ip] + BW) / widths[ip])
    plt.plot(xscale, np.real(spect), '--b')
    plt.plot(xscale, np.imag(spect), '--r')
    plt.plot(xscale, np.real(y), '-k')
    plt.plot(xscale, np.imag(y), '-k')
