'''
Python script to take a 
- Can run locally by passing command line arguments
- Can run with tyger by stdin.buffer by piping
  cat ${input} | python mrd2recon.py ${parameter arguments}
'''

import numpy as np
from typing import BinaryIO
import argparse
from mrd2recon import reconstruct_mrs

from MRSreader import MRSdata
import sys
sys.path.insert(0, 'python')
import mrd


if __name__ == "__main__":
    # read command line arguments to set spectral fitting parameters for each metabolite
    parser = argparse.ArgumentParser(description='Reconstruct MRS data from raw MRD file to processed MRD2 file')
    parser.add_argument('-i', '--input', type=str, required=False, help='Input file, defaults to stdin')
    parser.add_argument('-o', '--output', type=str, required=False, help='Output file, defaults to stdout')
    parser.add_argument('-w', '--wigglefactor', type=float, required=False, default=1.0, help='Wiggle factor for spectral fitting, default 1.0')
    parser.add_argument('-p', '--peakoffsets', type=str, required=False, help='List of peak offsets in ppm. Suffix with * to indicate small peak')
    
    args = parser.parse_args()
    input = open(args.input, 'rb') if args.input is not None else sys.stdin.buffer
    output = open(args.output, 'wb') if args.output is not None else sys.stdout.buffer
    # peakoffsets, peaknames, biggestpeaklist, npeaks
    peakoffsets = []
    biggestpeakidx = []
    peaknames = []
    peak_args = args.peakoffsets.split(' ')
    for i, item in enumerate(peak_args):
        if item.startswith('-'):
            peakoffsets.append(float(peak_args[i + 1]))
            if item.endswith('*'):  # small peak with suffix *
                peaknames.append(item[1:-1])
            else:                   # big peak without suffix
                peaknames.append(item[1:])
                biggestpeakidx.append(len(peaknames) - 1)
    peakoffsets = np.array(peakoffsets)
    

    # now look for specification of metabolite peaks
    # BA's cirrhrat is -bic* 0.0 -urea 2.3 -pyr 9.7 -ala* 15.2 -hyd* 18.1 -lac 21.8
    # SZ's mouse kidney is -bic 0.0 -urea 2.3 -pyr 9.7 -ala 15.2 -poop 15.9 -hyd 18.1 -lac 21.8
    # BA's spectra is -bic* -0.4 -urea 2.1 -urea2* 2.3 -pyr 9.7 -ala* 15.2 -hyd* 18.1 -lac 21.8 -w 0.5
    # DT's spectra is -urea 0.0 -KIC 8.6 -leu* 13.0 -hyd* 18.1 -?* 21.8 -w 1.5
    reconstruct_mrs(input, output, biggestpeakidx, peakoffsets, peaknames)