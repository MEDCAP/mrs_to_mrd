import sys
import numpy as np
import matplotlib.pyplot as plt
import os
from pathlib import Path
from typing import BinaryIO
import argparse

# path to mrd python package is 'root/mrd-fork/python' 
sys.path.append(os.path.join(Path(__file__).parent, 'mrd-fork','python'))
import mrd
# from watchdog.observers import Observer
# from watchdog.events import FileSystemEventHandler
from scipy.ndimage import zoom

def plot_mrd(input: BinaryIO, filename: str):
    current_header = []
    current_streamables_list = []
    current_pulse_list = []
    current_gradient_list = []
    current_acquisition_list = []
    current_image_list = []

    # read all streamable items and sort them into types
    with mrd.BinaryMrdReader(input) as w:
        current_header = w.read_header()
        current_streamables_list = list(w.read_data())
        current_pulse_list = [x.value for x in current_streamables_list if type(x.value) == mrd.Pulse]
        current_pulse_list.sort(key = lambda x: x.head.pulse_time_stamp_ns)
        current_gradient_list = [x.value for x in current_streamables_list if type(x.value) == mrd.Gradient]
        current_gradient_list.sort(key = lambda x: x.head.gradient_time_stamp_ns)
        current_acquisition_list = [x.value for x in current_streamables_list if type(x.value) == mrd.Acquisition]
        current_acquisition_list.sort(key = lambda x: x.head.acquisition_time_stamp_ns)
        current_image_list = [x.value for x in current_streamables_list if type(x.value) == mrd.Image]
        # images do not necessarily have a time, so don't sort them
        print('found', len(current_pulse_list), 'pulses,', len(current_gradient_list), 'gradients,', \
            len(current_acquisition_list), 'acquisitions,', len(current_image_list), 'images')
    fig1 = plt.figure()
    fig1.suptitle('Grad, Acq, Pulse File: ' + filename)
    plt.tight_layout()
    if(current_pulse_list):
        # draw pulses
        plt.subplot(3,1,1)
        for p in current_pulse_list:
            for q in range(p.amplitude.shape[0]):
                plt.plot(((np.array(range(p.amplitude.shape[1] + 2)) - 1) * \
                        p.head.sample_time_ns + p.head.pulse_time_stamp_ns) * \
                        1.0E-9, np.concatenate((np.array([0]), p.amplitude[q, :], np.array([0]))))
        plt.xlabel('time (s)')
        plt.ylabel('pulse amplitude (V)')
        plt.title('RF Pulses')
    if(current_gradient_list):
        # draw gradients
        plt.subplot(3,1,2)
        for g in current_gradient_list:
            plt.plot((np.array(range(len(g.data))) * g.head.gradient_sample_time_ns + g.head.gradient_time_stamp_ns) * \
                    1.0E-9, g.data)
            # plt.plot((np.array(range(len(g.data))) * g.head.gradient_sample_time_ns + g.head.gradient_time_stamp_ns) * \
            #         1.0E-9, g.data[:,1])
            # plt.plot((np.array(range(len(g.data))) * g.head.gradient_sample_time_ns + g.head.gradient_time_stamp_ns) * \
            #         1.0E-9, g.data[:,2])
        plt.xlabel('time (s)')
        plt.ylabel('gradient amplitude (mT/m)')
        plt.title('Gradients')
    if(current_acquisition_list):
        # draw acquisition raw data
        plt.subplot(3,1,3)
        for a in current_acquisition_list:
            plt.plot(np.array(range(len(a.data))) * a.head.sample_time_ns + a.head.acquisition_time_stamp_ns, \
                    np.real(a.data[:,0]))
            plt.plot(np.array(range(len(a.data))) * a.head.sample_time_ns + a.head.acquisition_time_stamp_ns, \
                    np.imag(a.data[:,0]))
        plt.xlabel('time (s)')
        plt.ylabel('signal (uV)')
        plt.title('Acquisitions')
    plt.xticks([])
    plt.yticks([])
    if(current_image_list):
        nimg = 0
        for i in current_image_list:
            if(i.head.image_type == mrd.ImageType.BITMAP):
                # show each bitmap figure in a separate window
                bmp = np.squeeze(i.data)
                unpacked_bmp = np.zeros([bmp.shape[0], bmp.shape[1], 4], dtype=np.uint8)
                unpacked_bmp[:,:,0] = ((bmp & 0x000000FF) / 2**0).astype(np.uint8)
                unpacked_bmp[:,:,1] = ((bmp & 0x0000FF00) / 2**8).astype(np.uint8)
                unpacked_bmp[:,:,2] = ((bmp & 0x00FF0000) / 2**16).astype(np.uint8)
                unpacked_bmp[:,:,3] = ((bmp & 0xFFFF0000) / 2**24).astype(np.uint8)
                fig2 = plt.figure()
                fig2.suptitle('File: ' + filename)
                plt.xticks([])
                plt.yticks([])
                plt.imshow(unpacked_bmp)
            else:
                # this is a metabolite image (for epsi)
                nimg += 1
                metshape = np.shape(i.data)[-3:] + (nimg,)
        zf = 2 # zoom factor for image interpolation
        met = np.zeros(metshape, dtype=np.float32)
        if(metshape[0] * metshape[1] > 1):
            [height, width, nmet, nimg] = metshape
            plotfig = np.zeros((height * nmet * zf, width * nimg * zf))
            iimg = 0
            for i in current_image_list:
                if(i.head.image_type == mrd.ImageType.BITMAP):
                    continue
                # put this image's data into the global metabolite array (for plotting and saving)
                met[: ,: ,:, iimg] = np.squeeze(i.data)
                iimg += 1
            for imet in range(nmet):
                # find maximum (positive or negative) for scaling 
                themax = np.max(met[:, :, imet, :])
                for iimg in range(nimg):
                    thisimg = met[:, :, imet, iimg]
                    plotfig[(imet * height * zf):((imet + 1) * height * zf), (iimg * width * zf):((iimg + 1) * \
                            width * zf)] = zoom(np.rot90(thisimg), zf, order=2) / themax
                    plotfig[imet * height * zf, :] = 1
            fig3 = plt.figure()
            fig3.suptitle('metabolic images, File: ' + filename)
            plt.imshow(plotfig, cmap='gray')
            plt.xticks([])
            plt.yticks([])
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot MRD file contents')
    parser.add_argument('-i', '--input', type=str, required=False, help='Input file, defaults to stdin')
    args = parser.parse_args()
    if args.input is None:
        input = sys.stdin.buffer
        filename = ''
    else:
        input = open(args.input, 'rb')
        filename = Path(args.input).parent.stem
    plot_mrd(input, filename)
