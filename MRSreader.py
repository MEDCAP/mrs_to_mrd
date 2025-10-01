import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.optimize import minimize

MRSdatadebug = False

class MRSdata:
    def __init__(self):
        self.samples = 0
        self.views = 0
        self.sliceviews = 0
        self.slices = 0
        self.echoes = 0
        self.nex = 0
        self.basefreq = 0
        self.pplfile = ''
        self.sampleperiod = 0
        self.nslc = 0
        self.acqstarttime = 0
        self.alpha = 0
        self.navg = 0
        self.nswitch = 0
        self.npswitch = 0
        self.FOVoff = [0.0, 0.0, 0.0]
        self.FOVaspect = 0.0
        self.FOV = 0.0
        self.tr = 0.0
    def mread3d(self, f):
        if(MRSdatadebug):
            print('reading MRS file ', f)
        containingfolder = '/'.join(f.split('/')[:-1])
        d = os.listdir(containingfolder)
        # first get base frequency from SPR file
        for auxf in d:
            if(auxf.find('.SPR') > 0):
                fd = open(containingfolder + '/' + auxf, 'rb')
                fdstr = str(fd.read())
                freqidx = fdstr.find('FREQ')
                if(freqidx > 0):
                    self.basefreq = int(float(fdstr[(freqidx + 5):(freqidx + 19)]) * 1.0E+6 + 0.5)
                    if(MRSdatadebug):
                        print('   setting base frequency to ', self.basefreq, ' Hz')
            if self.basefreq == 0:
                self.basefreq = 74941736
        # now read MRS file
        fd = open(f, 'rb')
        fdbytes = fd.read()
        self.samples = np.frombuffer(fdbytes[0:4], dtype = 'int32')[0]
        self.views = np.frombuffer(fdbytes[4:8], dtype = 'int32')[0]
        self.sliceviews = np.frombuffer(fdbytes[8:12], dtype = 'int32')[0]
        self.slices = np.frombuffer(fdbytes[12:16], dtype = 'int32')[0]
        self.type = np.frombuffer(fdbytes[18:20], dtype = 'int16')[0] 
        self.echoes = np.frombuffer(fdbytes[152:156], dtype = 'int32')[0]
        self.nex = np.frombuffer(fdbytes[156:160], dtype = 'int32')[0]
        totalpts = self.samples * self.views * self.sliceviews * self.slices * self.echoes * self.nex
        dstart = 512
        if self.type == 3:
            dend = dstart + totalpts * 2
            rawdata = np.frombuffer(fdbytes[dstart:dend], dtype = 'int16')
        elif self.type == 16:
            dend = dstart + totalpts * 2
            d = np.frombuffer(fdbytes[dstart:dend], dtype = 'uint8')
            rawdata = d[::2] + 1j * d[1::2]
        elif self.type == 17:
            dend = dstart + totalpts * 2
            d = np.frombuffer(fdbytes[dstart:dend], dtype = 'int8')
            rawdata = d[::2] + 1j * d[1::2]
        elif self.type == 18 or self.type == 19:
            dend = dstart + totalpts * 4
            d = np.frombuffer(fdbytes[dstart:dend], dtype = 'int16')
            rawdata = d[::2] + 1j * d[1::2]
        elif self.type == 20:
            dend = dstart + totalpts * 8
            rawdata = d[::2] + 1j * d[1::2]
            d = np.frombuffer(fdbytes[dstart:dend], dtype = 'int32')
        elif self.type == 21:
            dend = dstart + totalpts * 8
            d = np.frombuffer(fdbytes[dstart:dend], dtype = 'float32')
            rawdata = d[::2] + 1j * d[1::2]
        elif self.type == 22:
            dend = dstart + totalpts * 16
            d = np.frombuffer(fdbytes[dstart:dend], dtype = 'float64')
            rawdata = d[::2] + 1j * d[1::2]
        else:
            print('unknown data format')
            return
        # parameters describes settings, appended to the end of the file
        self.rawdata = np.reshape(rawdata, (self.samples, self.views, self.sliceviews, self.slices, \
                self.echoes, self.nex), order = 'F')
        if(MRSdatadebug):
            print('   reading data ', self.rawdata.shape, ' nsamp x nview x nslcview x nslc x necho x nex')
        self.parameters = str(fdbytes[dend:])
        # set ppm file name
        endidx = self.parameters.find('.ppl')
        if endidx == -1:
            return("")
        for beginidx in range(endidx, 0, -1):
            if(self.parameters[beginidx] == '/' or self.parameters[beginidx] == '\\'):
                 break
        self.pplfile = self.parameters[(beginidx + 1):endidx]
        if(MRSdatadebug):
            print('   setting ppl file name =', self.pplfile)
        # set sample period in 1/10ths of a microsecond
        beginidx = self.parameters.find('SAMPLE_PERIOD')
        beginidx += self.parameters[beginidx:].find(',') + 1
        endidx =  min(beginidx + self.parameters[(beginidx + 1):].find(',') + 1, \
                beginidx + self.parameters[(beginidx + 1):].find('\\r') + 1)
        self.sampleperiod = int(self.parameters[beginidx:endidx])
        if(MRSdatadebug):
            print('   setting sample period to ', self.sampleperiod, ' tenths of a microsecond')
        # set number of slices
        beginidx = self.parameters.find('NO_SLICES')
        beginidx += self.parameters[beginidx:].find(',') + 1
        endidx =  min(beginidx + self.parameters[(beginidx + 1):].find(',') + 1, \
                beginidx + self.parameters[(beginidx + 1):].find('\\r') + 1)
        try:
            self.nslc = int(self.parameters[beginidx:endidx])
            if(MRSdatadebug):
                print('   setting num slices to ', self.nslc)
        except:
            print('   nslc not specified')
        # set number of averages
        beginidx = self.parameters.find('NO_AVERAGES')
        beginidx += self.parameters[beginidx:].find(',') + 1
        endidx =  min(beginidx + self.parameters[(beginidx + 1):].find(',') + 1, \
                beginidx + self.parameters[(beginidx + 1):].find('\\r') + 1)
        self.navg = int(self.parameters[beginidx:endidx])
        if(MRSdatadebug):
            print('   setting num averages to ', self.navg)
        # set alpha
        beginidx = self.parameters.find('alpha')
        beginidx += self.parameters[beginidx:].find(',') + 1
        endidx =  min(beginidx + self.parameters[(beginidx + 1):].find(',') + 1, \
                beginidx + self.parameters[(beginidx + 1):].find('\\r') + 1)
        try:
            self.alpha = int(self.parameters[beginidx:endidx])
            if(MRSdatadebug):
                print('   setting flip angle to ', self.alpha)
        except:
            print('   alpha not specified')
        # set tr
        beginidx = self.parameters.find('tr,')
        beginidx += self.parameters[beginidx:].find(',') + 1
        endidx =  min(beginidx + self.parameters[(beginidx + 1):].find(',') + 1, \
                beginidx + self.parameters[(beginidx + 1):].find('\\r') + 1)
        try:
            self.tr = int(self.parameters[beginidx:endidx])
            if(MRSdatadebug):
                print('   setting tr to ', self.tr)
        except:
            print('   tr not specified')
        # set number of switches
        beginidx = self.parameters.find('no_switches')
        beginidx += self.parameters[beginidx:].find(',') + 1
        endidx =  min(beginidx + self.parameters[(beginidx + 1):].find(',') + 1, \
                beginidx + self.parameters[(beginidx + 1):].find('\\r') + 1)
        try:
            self.nswitch = int(self.parameters[beginidx:endidx])
            if(MRSdatadebug):
                print('   setting number of switches ', self.nswitch)
        except:
            print('   nswitches not specified')
        # set number of points per switch
        beginidx = self.parameters.find('no_pts_switch')
        beginidx += self.parameters[beginidx:].find(',') + 1
        endidx =  min(beginidx + self.parameters[(beginidx + 1):].find(',') + 1, \
                beginidx + self.parameters[(beginidx + 1):].find('\\r') + 1)
        try:
            self.nppswitch = int(self.parameters[beginidx:endidx])
            if(MRSdatadebug):
                print('   setting flip angle to ', self.alpha)
        except:
            print('   points per switch not specified')
        # set FOV offsets
        beginidx = self.parameters.find('FOV_OFFSETS')
        beginidx += self.parameters[beginidx:].find(',') + 1
        endidx =  min(beginidx + self.parameters[(beginidx + 1):].find(',') + 1, \
                beginidx + self.parameters[(beginidx + 1):].find('\\r') + 1)
        try:
            self.FOVoff[0] = float(self.parameters[beginidx:endidx]) / 1000.0
            beginidx = endidx + 1
            endidx =  min(beginidx + self.parameters[(beginidx + 1):].find(',') + 1, \
                    beginidx + self.parameters[(beginidx + 1):].find('\\r') + 1)
            self.FOVoff[1] = float(self.parameters[beginidx:endidx]) / 1000.0
            beginidx = endidx + 1
            endidx =  min(beginidx + self.parameters[(beginidx + 1):].find(',') + 1, \
                    beginidx + self.parameters[(beginidx + 1):].find('\\r') + 1)
            self.FOVoff[2] = float(self.parameters[beginidx:endidx]) / 1000.0
            if(MRSdatadebug):
                print('   setting FOV offsets (m) to ', self.FOVoff[0], self.FOVoff[1], self.FOVoff[2])
        except:
            print('   FOV offsets not specified')
        # set FOV
        beginidx = self.parameters.find('FOV') + 3
        if(beginidx >  3):
            endidx =  min(beginidx + self.parameters[(beginidx + 1):].find(',') + 1, \
                    beginidx + self.parameters[(beginidx + 1):].find('\\r') + 1)
            try:
                self.FOV = float(self.parameters[beginidx:endidx]) / 1000.0
                if(MRSdatadebug):
                    print('   setting FOV (m)', self.FOV)
            except:
                print('   FOV not specified')
        # set aspect ratio
        beginidx = self.parameters.find('aspect_ratio')
        beginidx += self.parameters[beginidx:].find(',') + 1
        endidx =  min(beginidx + self.parameters[(beginidx + 1):].find(',') + 1, \
                beginidx + self.parameters[(beginidx + 1):].find('\\r') + 1)
        try:
            self.FOVaspect = float(self.parameters[beginidx:endidx])
            if(MRSdatadebug):
                print('   setting FOV aspect ratio to ', self.FOVaspect)
        except:
            print('   FOV aspect not specified')
        # set acquisition start time
        beginidx = self.parameters.find('AcquisitionStartTime') + len('AcquisitionStartTime')
        endidx =  beginidx + self.parameters[(beginidx + 1):].find('\\r') + 1
        self.acqstarttime = int(self.parameters[beginidx:endidx])
        if(MRSdatadebug):
            print('   setting acq start time to ', self.acqstarttime, ' since some time in units of 100ns')
