from fractions import Fraction
import sys
import numpy as np
import matplotlib

import matplotlib.pyplot as plt
import os
from pathlib import Path
from typing import BinaryIO
import argparse

# path to mrd python package is 'root/mrd-fork/python' 
import mrd
# from watchdog.observers import Observer
# from watchdog.events import FileSystemEventHandler
from scipy.ndimage import zoom

def plot_mrd(input: BinaryIO):
    current_streamables_list = []
    current_acquisition_list = []
    current_image_list = []
    current_phantom_list = []
    if phantom:
        with mrd.BinaryMrdReader(phantom) as ph_reader:
            ph_reader.read_header()
            current_phantom_list = [x.value for x in list(ph_reader.read_data()) if type(x.value) == mrd.Image]

    # read all streamable items and sort them into types
    if input:
        with mrd.BinaryMrdReader(input) as reader:
            current_header = reader.read_header()
            current_streamables_list = list(reader.read_data())
            current_acquisition_list = [x.value for x in current_streamables_list if type(x.value) == mrd.Acquisition]
            current_acquisition_list.sort(key = lambda x: x.head.acquisition_time_stamp_ns)
            current_image_list = [x.value for x in current_streamables_list if type(x.value) == mrd.Image]
    nimg = 0
    for i in current_image_list:
        if(i.head.image_type == mrd.ImageType.COMPLEX):
            # show each bitmap figure in a separate window
            bmp = np.squeeze(i.data)
            if bmp.ndim < 2:
                continue
            unpacked_bmp = np.zeros([bmp.shape[0], bmp.shape[1], 4], dtype=np.uint8)
            unpacked_bmp[:,:,0] = ((bmp & 0x000000FF) / 2**0).astype(np.uint8)
            unpacked_bmp[:,:,1] = ((bmp & 0x0000FF00) / 2**8).astype(np.uint8)
            unpacked_bmp[:,:,2] = ((bmp & 0x00FF0000) / 2**16).astype(np.uint8)
            unpacked_bmp[:,:,3] = ((bmp & 0xFFFF0000) / 2**24).astype(np.uint8)
            fig2 = plt.figure()
            fig2.suptitle(f'File: {input.name}')
            plt.xticks([])
            plt.yticks([])
            plt.imshow(unpacked_bmp)
            plt.show()
        else:
            # this is a metabolite image (for epsi)
            nimg += 1
            metshape = np.shape(i.data)[-3:] + (nimg,)  # metshape=(y, x, frequency, rep)

#    for i in current_phantom_list:
#        phantom_data = i.data[0, 0, :, :, 0]    # phantom_data shape=(y, x, frequency=0)

    if metshape[0] * metshape[1] > 1:
        zf = 2
        plotfig = np.zeros((1, 1), dtype=np.float32)
        met = np.zeros(metshape, dtype=np.float32)
        normalized_met = np.zeros_like(met)
        if(metshape[0] * metshape[1] > 1):
            [height, width, nmet, nimg] = metshape
            print('h-w-met-img', metshape)
            plotfig = np.zeros((height * nmet * zf, width * nimg * zf))
    #        phantomplotfig = np.zeros_like(plotfig)
            iimg = 0
            for i in current_image_list:
                if(i.head.image_type == mrd.ImageType.COMPLEX):
                    continue
                # put this image's data into the global metabolite array (for plotting and saving)
                met[: ,: ,:, iimg] = i.data[0, 0, ...]
                for imet in range(metshape[2]):
                    normalized_met[:, :, imet, iimg] = i.data[0, 0, :, :, imet]# / zoom(phantom_data, 2/3)
                iimg += 1
            for imet in range(nmet):
                # find maximum (positive or negative) for scaling 
                themax = np.max(met[:, :, imet, :])
                phmax = np.max(normalized_met[:, :, imet, :])
                for iimg in range(nimg):
                    thisimg = met[:, :, imet, iimg]
                    plotfig[(imet * height * zf):((imet + 1) * height * zf), (iimg * width * zf):((iimg + 1) * \
                            width * zf)] = zoom(thisimg, zf, order=2) / themax
                    plotfig[imet * height * zf, :] = 1

#                phantomplotfig[(imet * height * zf):((imet + 1) * height * zf), (iimg * width * zf):((iimg + 1) * \
#                        width * zf)] = zoom(thisimg, zf, order=2) / phmax
#                phantomplotfig[imet * height * zf, :] = 1e-16*8

        plt.imshow(plotfig)
        plt.title(f'Filename: {input.name}')
        plt.show()
        return plotfig

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot MRD file contents')
    parser.add_argument('-i', '--input', type=Path, required=False, help='Input file, defaults to stdin')
    parser.add_argument('-p', '--phantom', type=str, required=False, help='Phantom file')
    args = parser.parse_args()

    if args.input:
        input = open(str(args.input), 'rb')
    else:
        input = ""    
    if args.phantom:
        phantom = open(str(args.phantom), 'rb')
    else:
        phantom = ""
    ratdata = plot_mrd(input)
    if ratdata is not None:
        plt.figure(1, figsize=(3, 10))
        plt.imshow(ratdata, cmap='gray')
        plt.title(f'Filename: {args.input.name} data rows=bic, urea, pyr, ala, hyd, lac')
        plt.colorbar(fraction=0.03, pad=0.02, shrink=0.6)
        plt.show()
