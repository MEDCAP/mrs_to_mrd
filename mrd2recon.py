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

def denoise_svd(matrix: np.ndarray, rank: int = None) -> tuple:
    """
    Truncated-SVD (low-rank) denoising of a complex signal matrix, following
    Francischello et al., NMR Biomed. 2021;34:e4285.

    The columns of `matrix` are the complex time-domain FIDs of a time series of
    spectra. In the absence of noise a series of metabolite spectra with fixed line
    shapes has rank ~= number of metabolites; additive Gaussian noise makes it full
    rank. Keeping only the first `rank` singular values (M_hat = U @ diag(S_r) @ Vh)
    recovers the low-rank signal subspace and discards the noise subspace.

    Args:
        - matrix: complex ndarray of shape (nsamples, nspectra); each column is a FID
        - rank: number of singular values to retain. If None, it is estimated with the
          Gavish-Donoho optimal hard threshold for a matrix with unknown noise level
          (thresh = omega(beta) * median(S), beta = min(m, n) / max(m, n)).
    Returns:
        - (matrix_hat, S, rank): denoised matrix, the singular values, and the rank used
    """
    U, S, Vh = np.linalg.svd(matrix, full_matrices=False)
    if rank is None:
        m, n = matrix.shape
        beta = min(m, n) / max(m, n)
        omega = 0.56 * beta ** 3 - 0.95 * beta ** 2 + 1.82 * beta + 1.43
        thresh = omega * np.median(S)
        rank = max(1, int(np.count_nonzero(S > thresh)))
    rank = min(rank, len(S))
    S_trunc = S.copy()
    S_trunc[rank:] = 0
    matrix_hat = (U * S_trunc) @ Vh
    return matrix_hat, S, rank

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



def _extract_fid(acq: mrd.Acquisition) -> np.ndarray:
    """Extract a single 1-D complex FID from an acquisition, trimming discard points."""
    fid = np.squeeze(np.asarray(acq.data)).astype(np.complex128).ravel()
    pre = acq.head.discard_pre or 0
    post = acq.head.discard_post or 0
    if pre or post:
        fid = fid[pre:len(fid) - post]
    return fid

def reconstruct_spectral(head: mrd.Header, input: Iterable[mrd.Acquisition],
                         line_broadening: float, rank: int = None) -> Iterable[mrd.StreamItem]:
    """
    Reconstruct a single-voxel FID spectral time series with truncated-SVD denoising
    (Francischello et al., NMR Biomed. 2021;34:e4285).

    Each acquisition is one complex FID acquired at a time point of the metabolic-flux
    experiment. The FIDs are stacked as the columns of an (nsamples, nspectra) matrix,
    denoised via low-rank approximation (see `denoise_svd`), then line-broadened and
    Fourier transformed to spectra. Denoising is done on the raw complex signal, so no
    phase correction is required beforehand.

    Yields (as mrd.StreamItem):
        - each raw acquisition (preserved),
        - one "singular_values" NdArray (real) for knee inspection / rank selection,
        - one "global_spect" complex spectrum NdArray per time point (denoised).
    Args:
        - head: mrd.Header (used for the ppm axis via h1resonance_frequency_hz)
        - input: iterable of mrd.Acquisition (one FID per time point)
        - line_broadening: Lorentzian line-broadening factor in Hz
        - rank: singular values to retain (None -> Gavish-Donoho auto estimate)
    """
    acqs = list(input)
    if not acqs:
        return
    fids = [_extract_fid(acq) for acq in acqs]
    nsamples = min(f.shape[0] for f in fids)
    # stack as (nsamples, nspectra); truncate to common length for safety
    matrix = np.stack([f[:nsamples] for f in fids], axis=1)

    matrix_hat, singular_values, used_rank = denoise_svd(matrix, rank)
    print(f"SVD denoising: matrix={matrix.shape}, rank={used_rank}, "
          f"top singular values={np.round(singular_values[:min(10, len(singular_values))], 4)}",
          file=sys.stderr)

    # frequency / ppm axis from the dwell time
    dwell_s = acqs[0].head.sample_time_ns / 1e9
    freq_hz = np.fft.fftshift(np.fft.fftfreq(nsamples, d=dwell_s))
    carrier_mhz = (head.experimental_conditions.h1resonance_frequency_hz or 0) / 1e6
    if carrier_mhz > 0:
        xaxis = freq_hz / carrier_mhz          # ppm
    else:
        xaxis = freq_hz                        # Hz (no carrier available)
    xscale = [mrd.ArrayMetaValue.Float64(float(x)) for x in xaxis]

    # line-broadening apodization applied to the denoised FIDs before FFT
    t = np.arange(nsamples) * dwell_s
    apod = np.exp(-np.pi * line_broadening * t)

    # preserve the raw acquisitions
    for acq in acqs:
        yield mrd.StreamItem.Acquisition(acq)

    # singular values (for inspecting the knee and re-running with an explicit --rank)
    yield mrd.StreamItem.NdArrayDouble(mrd.NDArrayDouble(
        head=mrd.NDArrayHeader(
            dimension_labels=[mrd.ArrayDimension.N],
            array_type=mrd.ArrayType.USER_MAP,
            meta=mrd.ArrayMeta({
                "description": [mrd.ArrayMetaValue.String("singular_values")],
                "svd_rank": [mrd.ArrayMetaValue.Int64(int(used_rank))],
            })),
        data=singular_values.astype(np.float64)))

    # one denoised spectrum per time point
    for j, acq in enumerate(acqs):
        spectrum = np.fft.fftshift(np.fft.fft(matrix_hat[:, j] * apod))
        yield mrd.StreamItem.NdArrayComplexDouble(mrd.NDArrayComplexDouble(
            head=mrd.NDArrayHeader(
                dimension_labels=[mrd.ArrayDimension.FREQUENCY],
                array_type=mrd.ArrayType.USER_MAP,
                meta=mrd.ArrayMeta({
                    "description": [mrd.ArrayMetaValue.String("global_spect")],
                    "xscale": xscale,
                    "repetition": [mrd.ArrayMetaValue.Int64(int(j))],
                    "line broadening(Hz)": [mrd.ArrayMetaValue.Float64(float(line_broadening))],
                    "acquisition_time_stamp_ns": [mrd.ArrayMetaValue.Int64(
                        int(acq.head.acquisition_time_stamp_ns or 0))],
                })),
            data=spectrum.astype(np.complex128)))

def append_spectral_header(header: mrd.Header, line_broadening: float, rank: int):
    """Record the denoising parameters on the header's user parameters."""
    if header.user_parameters is None:
        header.user_parameters = mrd.UserParametersType()
    header.user_parameters.user_parameter_double.append(
        mrd.UserParameterDoubleType(name="line_broadening_factor", value=line_broadening))
    if rank is not None:
        header.user_parameters.user_parameter_long.append(
            mrd.UserParameterLongType(name="svd_rank", value=int(rank)))

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

def reconstruct_mrs(input: BinaryIO, output: BinaryIO, line_broadening: float,
                    denoise: bool = False, rank: int = None):
    with mrd.BinaryMrdReader(input) as reader:
        with mrd.BinaryMrdWriter(output) as writer:
            header = reader.read_header()
            if denoise:
                # single-voxel FID spectral time series -> truncated-SVD denoising
                append_spectral_header(header, line_broadening, rank)
                writer.write_header(header)
                writer.write_data(
                    reconstruct_spectral(header,
                        acquisition_reader(reader.read_data()),
                        line_broadening, rank
                    )
                )
            else:
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
    parser.add_argument("-d", "--denoise", action="store_true", help="Apply truncated-SVD denoising to a FID spectral time series")
    parser.add_argument("-r", "--rank", type=int, default=None, required=False, help="Number of singular values to retain (default: auto via Gavish-Donoho)")
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
                    reconstruct_mrs(input, output, args.line_broadening, args.denoise, args.rank)
            else:
                raise ValueError(f"No mrd2 files found in {args.folder}")
    elif args.input:
        if not args.input.is_file():
            raise ValueError(f"{args.input} is not a file")
        if args.output is None:
            raise ValueError("--output must be specified with --input")
        if not args.output.parent.is_dir():
            raise ValueError(f"Output directory {args.output.parent} does not exist")
        input = open(args.input, "rb")
        output = open(args.output, "wb")
        reconstruct_mrs(input, output, args.line_broadening, args.denoise, args.rank)
    else:
        raise ValueError("Either --folder or --input must be specified")