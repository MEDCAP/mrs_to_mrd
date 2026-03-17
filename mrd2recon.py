"""
As new file comes in, convert to mrd2 file 

tyger args:
    - python mrd2recon.py
    - -i
    - $(INPUT_PIPE)
    - -o
    - $(OUTPUT_PIPE)

Across multiple raw.mrd2 files to reconstruct, have multiple jobs
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
from statistics import mode
import os
from scipy.optimize import minimize, Bounds
import sys
from pathlib import Path
from typing import Any, BinaryIO, Iterable, List, Union

import mrd

def apply_line_broadening(acq: mrd.Acquisition, line_broadening: float) -> np.ndarray:
    """
    Apply line broadening apodization A(t) = exp(-pi * LB * t) to each echo in an acquisition.
    Args:
        - acq: mrd.Acquisition object
        - line_broadening: float, line broadening factor in Hz
    Returns:
        - np.ndarray of shape (nsamples, nechoes) with apodization applied
    """
    nechoes = acq.head.idx.contrast
    totalppswitch = round(acq.samples() / nechoes)
    nsamples = totalppswitch - acq.head.discard_pre - acq.head.discard_post
    result = np.zeros((nsamples, nechoes), dtype='complex')
    for iecho in range(nechoes):
        tk = iecho * acq.head.sample_time_ns * totalppswitch / 1.0e+9
        extracted_samples = acq.data[0, (iecho * totalppswitch + acq.head.discard_pre):(iecho * totalppswitch + acq.head.discard_post + nsamples)]
        result[:, iecho] = extracted_samples * np.exp(-tk * line_broadening)
    return result

def reconstruct_epsi(head: mrd.Header, input: Iterable[mrd.Acquisition], line_broadening: float) -> Iterable[mrd.NdArrayDouble]:
    """
    Generator to reconstruct from acquisition for epsi without spectral fitting
    Args:
        - head: mrd.Header object
        - input: Iterable[mrd.Acquisition] object
        - line_broadening: float, line broadening factor in Hz
    Returns:
        - Iterable[mrd.NdArrayDouble] object
    acquisitio rawdata=[samples=1280, views=8, slcview=1, slice=1, echoes=1, nex=1], avg=1
    """
    def kspace_to_ndarray(kspace: np.ndarray, ref_acq: mrd.Acquisition) -> Iterable[Union[mrd.NdArrayDouble, mrd.NdArrayComplexDouble]]:
        """
        Reconstruct kspace(nviews, nsamples, nechoes) to ndarray(nviews, nsamples, nechoes) as an intermediate array
        Args:
            - kspace: kspace array in the shape of (views, samples, echoes)
            - ref_acq: 
        """
        nonlocal current_max, max_spect
        print(f"Reconstructing kspace for repetition={ref_acq.head.idx.repetition}", file=sys.stderr)
        dim = range(kspace.ndim)
        img = np.fft.fftshift(np.fft.fftn(kspace, axes=dim), axes=dim)
        print(head)
        for iview in range(kspace.shape[0]):
            for isample in range(kspace.shape[1]):
                # find the max abs across multi-echoes
                this_max = np.max(np.abs(img[iview, isample, :]))
                if this_max > current_max:
                    current_max = this_max
                    max_spect = np.copy(img[iview, isample, :])

        # at the last repetition, store max, global variable
        if ref_acq.head.idx.repetition == head.encoding[0].encoding_limits.repetition.maximum:
            # store max spectral
            yield mrd.NdArrayComplexDouble(
                head=mrd.NdArrayHeader(
                    dimension_labels=[mrd.ArrayDimension.CONTRAST],
                    array_type=mrd.ArrayType.USER_MAP
                ),
                meta=mrd.ArrayMeta({"description": [mrd.ArrayMetaValue.String("max spectral")]}),
                data=max_spect
            )

            # store last repetition as noise data
            yield mrd.NdArrayDouble(
                head=mrd.NdArrayHeader(
                    array_type=mrd.ArrayType.NOISE
                ),
                meta=mrd.ArrayMeta({"description": [mrd.ArrayMetaValue.String("noise")]}),
                data=np.mean(np.abs(img))
            )

        # flag IS_NAVIGATION_DATA is a phantom acquisition
        if ref_acq.head.flags & mrd.AcquisitionFlags.IS_NAVIGATION_DATA:
            print(f"Phantom data found at repetition={ref_acq.head.idx.repetition}", file=sys.stderr)
            yield mrd.NdArrayComplexDouble(
                head=mrd.NdArrayHeader(
                    dimension_labels=[mrd.ArrayDimension.Y, mrd.ArrayDimension.X, mrd.ArrayDimension.CONTRAST],
                    array_type=mrd.ArrayType.PHANTOM
                ),
                meta=mrd.ArrayMeta({"description": [mrd.ArrayMetaValue.String("urea phantoms")]}),
                data=img
            )
        # store acquisition including last repetition with receiver bandwidth as meta value
        else:
            dwell_time_ns = ref_acq.head.sample_time_ns / ref_acq.samples()
            receiver_bandwidth = 1 / (dwell_time_ns / 1e+9)
            yield mrd.NdArrayComplexDouble(
                head=mrd.NdArrayHeader(
                    dimension_labels=[mrd.ArrayDimension.Y, mrd.ArrayDimension.X, mrd.ArrayDimension.CONTRAST],
                    image_type=mrd.ArrayImageType.COMPLEX,
                    measurement_uid = ref_acq.head.measurement_uid,
                    average = ref_acq.head.idx.average,
                    repetition = ref_acq.head.idx.repetition,
                    acquisition_time_stamp_ns = ref_acq.head.acquisition_time_stamp_ns,
                ),
                meta=mrd.ArrayMeta(
                    {"receiver bandwidth(Hz)": [mrd.ArrayMetaValue.Float64(receiver_bandwidth)],
                     "line broadening(Hz)": [mrd.ArrayMetaValue.Float64(line_broadening)]}),
                data=img
            )
            # keep acquisition data in streamitem
    # store item.data for each kspace shape
    current_rep = -1
    kspace = None
    reference_acq = None
    enc = head.encoding[0]
    FIDPAD = 1
    current_max = -np.inf
    max_spect = None
    
    if enc.encoding_limits.phase != None:
        nviews = enc.encoding_limits.phase.maximum + 1
    if enc.encoding_limits.repetition != None:
        nreps = enc.encoding_limits.repetition.maximum + 1

    for acq in input:
        nechoes = acq.head.idx.contrast
        totalppswitch = round(acq.samples() / nechoes)
        nsamples = totalppswitch - acq.head.discard_pre - acq.head.discard_post     # samples = (npts_per_echo + 2 * discard) * (contrast=nechoes)
        # first acq of new repetition
        if acq.head.idx.repetition != current_rep:
            if kspace is not None and reference_acq is not None:
                yield from kspace_to_ndarray(kspace, reference_acq)
            kspace = np.zeros((nviews, nsamples, nechoes * FIDPAD), dtype='complex')
            kspace_apodized = np.zeros((nviews, nsamples, nechoes * FIDPAD), dtype='complex')
            reference_acq = acq
            current_rep = acq.head.idx.repetition
        # while idx.repetition=current_rep incrementally fill in kspace for each line
        if kspace is not None:
            view = acq.head.idx.kspace_encode_step_1 if acq.head.idx.kspace_encode_step_1 is not None else 0
            kspace[view, :, :] = apply_line_broadening(acq, line_broadening)
        # keep acquisition data in streamitem
        yield mrd.StreamItem.Acquisition(acq)
    if kspace is not None and reference_acq is not None:
        yield from kspace_to_ndarray(kspace, reference_acq)
        kspace = None
        reference_acq = None



def reconstruct_spectral():
    pass

def _ndarray_to_stream_item(arr: mrd.NdArray) -> mrd.StreamItem:
    """Map NdArray payload to StreamItem; dtype selects union arm (not the generic alias)."""
    dt = arr.data.dtype
    if dt == np.uint16:
        return mrd.StreamItem.NdArrayUint16(arr)
    if dt == np.int16:
        return mrd.StreamItem.NdArrayInt16(arr)
    if dt == np.uint32:
        return mrd.StreamItem.NdArrayUint32(arr)
    if dt == np.int32:
        return mrd.StreamItem.NdArrayInt32(arr)
    if dt == np.float32:
        return mrd.StreamItem.NdArrayFloat(arr)
    if dt == np.float64:
        return mrd.StreamItem.NdArrayDouble(arr)
    if dt == np.complex64:
        return mrd.StreamItem.NdArrayComplexFloat(arr)
    if dt == np.complex128:
        return mrd.StreamItem.NdArrayComplexDouble(arr)
    raise TypeError(f"Unsupported NdArray dtype for stream: {dt}")


# convert iterable mrd objects to mrd stream
def generate_stream(input: Iterable[Any]) -> Iterable[mrd.StreamItem]:
    for item in input:
        if isinstance(item, mrd.Acquisition):
            yield mrd.StreamItem.Acquisition(item)
        elif isinstance(item, mrd.NdArray):
            yield _ndarray_to_stream_item(item)
        else:
            continue

def acquisition_reader(input: Iterable[mrd.StreamItem]) -> Iterable[mrd.Acquisition]:
    """
    Generator to yield acquisition as iterable
    Assign NOISE_MEASUREMENT flag to avg > 1 as phantom data
    """
    for item in input:
        if not isinstance(item, mrd.StreamItem.Acquisition):
            continue
        # to ignore phantom acquisition with navg>1, uncomment below
        # if item.value.head.flags & mrd.AcqusitionFlags.IS_NAVIGATION_DATA:
        #     continue
        yield item.value

def append_header(header: mrd.Header, line_broadening: float):
    header.user_parameters.user_parameter_double.append(
        mrd.UserParameterDouble(name="line_broadening_factor", value=line_broadening)
    )
    return header

def reconstruct_mrs(input: BinaryIO, output: BinaryIO, line_broadening: float):
    with mrd.BinaryMrdReader(input) as reader:
        with mrd.BinaryMrdWriter(output) as writer:
            header = reader.read_header()
            append_header(header, line_broadening)
            writer.write_header(header)
            writer.write_data(
                generate_stream(
                    reconstruct_epsi(header,
                        acquisition_reader(reader.read_data())
                    )
                )
            )

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Reconstruct MRS data from mrd2 file")
    parser.add_argument("-f", "--folder", type=Path, required=False, help="Folder with mrs raw.mrd2 file")
    parser.add_argument("-i", "--input", type=Path, required=False, help="Input mrd2 file")
    parser.add_argument("-o", "--output", type=Path, required=False, help="Output mrd2 file")
    parser.add_argument("-lb", "--line-broadening", type=float, default=42, required=False, help="Line broadening factor in Hz")
    args = parser.parse_args()

    if args.folder and args.input:
        raise ValueError("Cannot specify both --folder and --input")
    elif args.folder:
        if not args.folder.is_dir():
            raise ValueError(f"{args.folder} is not a directory")
        else:
            # search for raw.mrd2 or targetfilename and turn into a list
            mrd2_filepaths: List[Path] = []
            for root, dirnames, filenames in os.walk(args.folder):
                if "raw.mrd2" in filenames:
                    mrd2_filepaths.append(Path(os.path.join(root, "raw.mrd2")))
            if len(mrd2_filepaths) > 0:
                # for raw filepath, run reconstruct with parameters from cmd line arguments
                for i, input_filepath in enumerate(mrd2_filepaths):
                    recon_filepath = input_filepath.with_name(input_filepath.parent.name + "_recon.mrd2")
                    print(f"Reconstructing {i+1}/{len(mrd2_filepaths)} at: {recon_filepath}", file=sys.stderr)
                    input = open(input_filepath, "rb")
                    output = open(recon_filepath, "wb")
                reconstruct_mrs(input, output, args.line_broadening)
            else:
                raise ValueError(f"No mrd2 files found in {args.folder}")
    elif args.input:
        if not args.input.is_file():
            raise ValueError(f"{args.input} is not a file")
        if not args.output.is_file():
            raise ValueError(f"{args.output} is not a file")
        input = open(args.input, "rb")
        output = open(args.output, "wb")
        reconstruct_mrs(input, output, args.line_broadening)
    else:
        raise ValueError("Either --folder or --input must be specified")