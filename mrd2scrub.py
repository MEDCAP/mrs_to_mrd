import numpy as np
import matplotlib.pyplot as plt
from statistics import mode
import os
from scipy.optimize import minimize, Bounds
from MRSreader import MRSdata
from lorn import lornfit, lor1fit, lorneval, lor1plot, lornputspect, lornpackx0, lornunpackx0, \
        lorngetpeakparams, lornputpeakparams
import sys
sys.path.insert(0, '../../mrd-fork/python')
import mrd

debugphasing = False
debuglorn = True

g = MRSdata()

basedir = "C:/Users/steph/Desktop/data/shurik_all_ischemia_data_081425"
#basedir = "C:/Users/steph/Desktop/data/cirrhrat data"
fidpad = 4 
lb = 42 # Hz

def findmrd2files(basedir):
    mrd2list = []
    if(not os.path.isdir(basedir) and basedir.find('.mrd2') > 0):
        return([basedir])
    if(os.path.isdir(basedir)):
        d = os.listdir(basedir)
        for f in d:
            mrd2list.extend(findmrd2files(basedir + '/' + f))
    return(mrd2list)

# check to see if a directory is specified
for iarg in range(len(sys.argv)):
    if(sys.argv[iarg] == '-f' and iarg < len(sys.argv) - 1):
        basedir = sys.argv[iarg + 1]
fnames = findmrd2files(basedir)
for f in fnames:
    print('removing', f)
    os.remove(f)
