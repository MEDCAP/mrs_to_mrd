"""
Reconstruct a converted .mrd2 file into spectra, Lorentzian peak fits and metabolite maps.

tyger args:
    - python mrd2recon.py
    - -i
    - $(INPUT_PIPE)
    - -o
    - $(OUTPUT_PIPE)

Across multiple converted .mrd2 files to reconstruct, have multiple jobs.

With --folder, every .mrd2 that is not itself a _recon.mrd2 is reconstructed to <name>_recon.mrd2
beside it. MRStomrd2 names what it converts <meas_id>_<sequence>.mrd2, so a directory can hold more
than one scan and there is no single filename to look for.

Which reconstruction runs is decided by the sequence, via mrs_organize.sequence_family, not by a
command line flag. Peaks are specified the way the legacy recon specified them, as a name with
optional modifiers followed by a chemical shift in ppm:

    python mrd2recon.py -i raw.mrd2 -o recon.mrd2 \
        -bic_tm 0.0 -urea 2.3 -pyr_s 9.7 -ala_tm 15.2 -hyd_tm 18.1 -lac_m 21.8

    _s  source peak, the injected substrate
    _t  tiny peak, not a candidate for "which peak is the tallest one"
    _m  a derived metabolite

Everything the reconstruction produces is written to the output stream as mrd NdArrays, and the
raw acquisitions are passed through unchanged. Nothing is plotted and nothing is written outside
the output file; use mrdplot.py to look at the result.

Requires the mrd python package on the `dev` branch of mrd-fork, which is where the NdArray
types live.
"""

import argparse
import math
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, BinaryIO, Iterable, List, Optional, Sequence, Tuple

import numpy as np

import mrd

from lorentzian_fitter import LorentzianFitter, candidate_centers, estimate_width_fwhm
from mrs_organize import sequence_family
from svd_denoise import denoise_svd

# a voxel spectrum is fitted only if it rises this far above the estimated noise floor
NOISE_THRESHOLD_MULTIPLIER = 3.0
# how far a voxel spectrum may be rolled when aligning it against the reference spectrum
PHASE_SEARCH_RANGE = 15
# spectral zero fill factor; 1 means the spectral axis is exactly one point per echo
FIDPAD = 1

# Where each echo of an EPSI readout starts. A gradient switch holds more points than the sequence
# keeps, the extras being the ramps either side of the flat top, and only the flat top carries the
# uniform kx spacing the transform assumes. discard_pre drops points as though the ramp were split
# evenly across both ends, which it is not, so the window is addressed as a pad back from it: echo i
# starts at i * totalppswitch + discard_pre - pad. MRStomrd2 records the ramp time as the 'tramp'
# user parameter in us, and the ramp is what the pad is derived from; a file converted before that
# was recorded falls back to this constant, i.e. to discard_pre exactly as the sequence set it
EPSIGRE_DEFAULT_PAD = 0
# leading samples of a readout to replace with zero, for a sequence whose first points are unusable
# for a reason the ramp time does not account for. Only the first echo can reach samples this early
EPSIGRE_ZERO_LEAD = 13


# ---------- peak specification -------------------------------------------


@dataclass(frozen=True)
class PeakSpec:
    """The metabolite peaks to fit, as named on the command line or recorded in a header."""
    names: List[str]
    offsets: np.ndarray     # chemical shifts in ppm, same order as names
    modifiers: List[str]    # the letters after the underscore, same order as names
    wigglefactor: float = 1.0

    def __len__(self) -> int:
        return len(self.names)

    @property
    def biggest_idx(self) -> List[int]:
        """Peaks eligible to be the tallest one in the spectrum, i.e. those without _t."""
        return [i for i, m in enumerate(self.modifiers) if "t" not in m]

    @property
    def metabolite_idx(self) -> List[int]:
        """Peaks marked _m, the ones a kinetic model would treat as products."""
        return [i for i, m in enumerate(self.modifiers) if "m" in m]

    @property
    def source_idx(self) -> Optional[int]:
        """The peak marked _s, the injected substrate, if one was named."""
        for i, m in enumerate(self.modifiers):
            if "s" in m:
                return i
        return None


def split_peak_args(argv: Sequence[str], reserved) -> Tuple[PeakSpec, List[str]]:
    """
    Pull the peak specifications out of argv, leaving the rest for argparse.

    This runs before argparse rather than using parse_known_args, because a peak's value may
    itself be negative ('-bic_tm -0.4'), which argparse would read as another option.
    A token is a peak if it starts with '-', is not one of argparse's own options, and is
    followed by something that parses as a float.
    Args:
        - argv: the arguments, without the program name
        - reserved: the option strings argparse owns, e.g. {'-i', '--input', ...}, which is
          also what keeps '-w 0.5' from being read as a peak named 'w'
    Returns:
        - (PeakSpec, the arguments argparse should see)
    """
    names: List[str] = []
    offsets: List[float] = []
    modifiers: List[str] = []
    remaining: List[str] = []

    i = 0
    while i < len(argv):
        token = argv[i]
        if token.startswith("-") and token not in reserved and i + 1 < len(argv):
            try:
                offset = float(argv[i + 1])
            except ValueError:
                remaining.append(token)
                i += 1
                continue
            body = token[1:]
            name, _, mods = body.partition("_")
            names.append(name)
            offsets.append(offset)
            modifiers.append(mods)
            i += 2
            continue
        remaining.append(token)
        i += 1

    # wigglefactor is one of argparse's own options, so it is filled in by the caller
    spec = PeakSpec(names=names,
                    offsets=np.asarray(offsets, dtype=float),
                    modifiers=modifiers)
    return spec, remaining


def peak_spec_to_header(header: mrd.Header, spec: PeakSpec) -> None:
    """Record the peak list on the header so a downstream consumer can recover it."""
    if header.user_parameters is None:
        header.user_parameters = mrd.UserParametersType()
    for i, name in enumerate(spec.names):
        full = f"{name}_{spec.modifiers[i]}" if spec.modifiers[i] else name
        header.user_parameters.user_parameter_double.append(
            mrd.UserParameterDoubleType(name=full, value=float(spec.offsets[i])))
    header.user_parameters.user_parameter_double.append(
        mrd.UserParameterDoubleType(name="wigglefactor", value=float(spec.wigglefactor)))


def peak_spec_from_header(header: mrd.Header) -> PeakSpec:
    """Recover the peak list that peak_spec_to_header wrote into header.user_parameters."""
    names: List[str] = []
    offsets: List[float] = []
    modifiers: List[str] = []
    wigglefactor = 1.0

    user = getattr(header, "user_parameters", None)
    params = getattr(user, "user_parameter_double", None) or [] if user is not None else []
    for p in params:
        if p.name == "wigglefactor":
            wigglefactor = float(p.value)
            continue
        name, _, mods = p.name.partition("_")
        names.append(name)
        offsets.append(float(p.value))
        modifiers.append(mods)

    return PeakSpec(names=names,
                    offsets=np.asarray(offsets, dtype=float),
                    modifiers=modifiers,
                    wigglefactor=wigglefactor)


def append_recon_header(header: mrd.Header, *,
                        line_broadening: float,
                        spec: PeakSpec,
                        denoise: bool,
                        rank: Optional[int],
                        phantom_scaling: float) -> mrd.Header:
    """Record what this reconstruction was asked to do on the header it passes through."""
    if header.user_parameters is None:
        header.user_parameters = mrd.UserParametersType()
    header.user_parameters.user_parameter_double.append(
        mrd.UserParameterDoubleType(name="line_broadening_factor", value=float(line_broadening)))
    if phantom_scaling != 1.0:
        header.user_parameters.user_parameter_double.append(
            mrd.UserParameterDoubleType(name="phantom_scaling", value=float(phantom_scaling)))
    if denoise and rank is not None:
        header.user_parameters.user_parameter_long.append(
            mrd.UserParameterLongType(name="svd_rank", value=int(rank)))
    if len(spec) > 0:
        peak_spec_to_header(header, spec)
    return header


# ---------- NdArray emission ---------------------------------------------


def meta_values(value: Any) -> List[mrd.ArrayMetaValue]:
    """Encode a python value as the list of ArrayMetaValue that mrd meta entries hold."""
    if isinstance(value, str):
        return [mrd.ArrayMetaValue.String(value)]
    if isinstance(value, (bool, int, np.integer)):
        return [mrd.ArrayMetaValue.Int64(int(value))]
    if isinstance(value, (float, np.floating)):
        return [mrd.ArrayMetaValue.Float64(float(value))]
    if isinstance(value, (list, tuple, np.ndarray)):
        return [v for item in value for v in meta_values(item)]
    raise TypeError(f"cannot encode {type(value)} as array meta")


def make_meta(**fields: Any) -> mrd.ArrayMeta:
    """Build an ArrayMeta from keyword values, dropping the ones that are None."""
    return {key: meta_values(value) for key, value in fields.items() if value is not None}


def ndarray_stream_item(arr: mrd.NdArray) -> mrd.StreamItem:
    """Map an NdArray onto its StreamItem union arm; the dtype picks the arm."""
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


def emit(data: np.ndarray, *,
         head: mrd.NdArrayHeader = None,
         dimension_labels: List = None,
         array_type: mrd.ArrayType = mrd.ArrayType.USER_MAP,
         **meta: Any) -> mrd.StreamItem:
    """
    Wrap an array as a stream-ready NdArray.

    Real data is normalised to float64 and complex data to complex128 so the union arm is
    always one of the two the readers on the other side expect. Meta goes on the NdArray,
    which is where the dev schema puts it, not on the header.
    """
    array = np.asarray(data)
    array = array.astype(np.complex128) if np.iscomplexobj(array) else array.astype(np.float64)
    if head is None:
        head = mrd.NdArrayHeader(dimension_labels=dimension_labels, array_type=array_type)
    return ndarray_stream_item(mrd.NdArray(head=head, data=array, meta=make_meta(**meta)))


# ---------- k-space assembly ---------------------------------------------


def header_user_long(header: mrd.Header, name: str) -> Optional[int]:
    """The named long user parameter, or None when the header does not carry it."""
    user = getattr(header, "user_parameters", None)
    params = getattr(user, "user_parameter_long", None) or [] if user is not None else []
    for param in params:
        if param.name == name:
            return int(param.value)
    return None


def epsi_leading_pad(header: mrd.Header,
                     acq: mrd.Acquisition,
                     override: Optional[int] = None) -> Tuple[int, str]:
    """
    How far back from discard_pre each echo's window starts, and why.

    The ramp at the head of a switch is not on the flat top, so the usable window opens
    ceil(tramp / dwell) samples into the switch. discard_pre assumes the discarded points split
    evenly across both ends instead, and the pad is the difference between the two, which is what
    apply_line_broadening subtracts from every echo's start.
    Args:
        - acq: any acquisition of the readout, read for its dwell time and discard count
        - override: a pad given on the command line, which wins over the header
    Returns:
        - (pad, a one line account of where it came from, for the log)
    """
    if override is not None:
        return override, f"pad {override} from the command line"
    tramp_us = header_user_long(header, "tramp")
    if tramp_us is None:
        return EPSIGRE_DEFAULT_PAD, (f"pad {EPSIGRE_DEFAULT_PAD}, the header records no ramp time; "
                                     f"convert again to record one")
    if not acq.head.sample_time_ns:
        return EPSIGRE_DEFAULT_PAD, f"pad {EPSIGRE_DEFAULT_PAD}, the acquisition records no dwell time"
    # ceil, because a window opening part way through the last ramp sample still includes it
    ramp_samples = int(math.ceil(tramp_us * 1000.0 / acq.head.sample_time_ns))
    discard_pre = acq.head.discard_pre or 0
    pad = discard_pre - ramp_samples
    return pad, (f"pad {pad} from tramp={tramp_us}us over a {acq.head.sample_time_ns}ns dwell, "
                 f"i.e. a {ramp_samples} sample ramp against discard_pre={discard_pre}")


def switch_layout(acq: mrd.Acquisition) -> Tuple[int, int, int]:
    """
    How one EPSI readout divides into gradient switches.

    The switch count comes from user_int, where the converter records it, because idx.contrast
    carries the echo index the sequence acquired rather than the switches packed into one readout.
    user_int holds the whole width of a switch, ramps included; the kept width is what is left of
    it once the discard points come off both ends.
    Args:
        - acq: one acquisition of the readout
    Returns:
        - (switches, total points in one switch, points kept per switch)
    Raises:
        - ValueError on an acquisition converted before the switch layout was recorded
    """
    user_int = list(acq.head.user_int or [])
    if len(user_int) < 2 or int(user_int[0]) <= 0:
        raise ValueError("this acquisition records no switch layout in user_int; convert the scan "
                         "again with the current MRStomrd2")
    nswitch, totalppswitch = int(user_int[0]), int(user_int[1])
    kept = totalppswitch - (acq.head.discard_pre or 0) - (acq.head.discard_post or 0)
    return nswitch, totalppswitch, kept


def apply_line_broadening(acq: mrd.Acquisition,
                          line_broadening: float,
                          *,
                          leading_pad: int = 0,
                          zero_lead: int = 0) -> np.ndarray:
    """
    Split one EPSI readout into its switches and apply the line broadening apodization.

    An EPSI readout packs every switch into a single acquisition; each switch carries one point of
    the spectral dimension, so the apodization decays over switches, not over the points within
    one of them.

    leading_pad and zero_lead express a leading zero fill without copying the readout: the switch
    boundaries are placed as if leading_pad zeros had been prepended, and any sample whose
    position in the original readout is below zero_lead reads as zero. Sampling past the end of
    the readout, which the shift can cause for the last switch, also reads as zero.
    Args:
        - acq: one acquisition, one phase encode line
        - line_broadening: line broadening factor in Hz
        - leading_pad: zeros notionally prepended to the readout
        - zero_lead: leading samples of the original readout to replace with zero
    Returns:
        - (kept points, switches) complex array, discard points trimmed off each switch
    """
    # example data: 64 switches of 28 points over 1792 samples, 12 of each switch kept
    nswitch, totalppswitch, kept = switch_layout(acq)
    discard_pre = acq.head.discard_pre or 0
    result = np.zeros((kept, nswitch), dtype='complex')
    offsets = np.arange(kept)
    for iswitch in range(nswitch):
        tk = iswitch * acq.head.sample_time_ns * totalppswitch / 1.0e+9
        start = iswitch * totalppswitch + discard_pre - leading_pad
        source = start + offsets
        valid = (source >= zero_lead) & (source < acq.samples())
        result[valid, iswitch] = acq.data[0, source[valid]] * np.exp(-tk * line_broadening)
    return result


def spectral_axis(header: mrd.Header, acq: mrd.Acquisition) -> Tuple[np.ndarray, float, float]:
    """
    The spectral axis of an EPSI readout, in ppm.

    sample_time_ns is the dwell time of a single point, and one spectral point is acquired per
    readout switch, so the spectral sampling interval is a whole switch. The axis deliberately
    runs from 0 rather than being centred on zero: the peak centers are placed modulo the
    spectral width and the Lorentzian model carries explicit +-BW wraparound terms.
    Returns:
        - (xscale in ppm, spectral bandwidth in Hz, spectral bandwidth in ppm)
    """
    nswitch, totalppswitch, _ = switch_layout(acq)
    spectral_bw_hz = 1.0e+9 / (acq.head.sample_time_ns * totalppswitch)
    # the converter writes the 13C frequency here despite the field name, since that is the
    # frequency these spectra were actually acquired at
    center_freq_hz = header.experimental_conditions.h1resonance_frequency_hz
    bw_ppm = spectral_bw_hz / center_freq_hz * 1.0e+6
    nfreq = nswitch * FIDPAD
    xscale = np.arange(nfreq) / nfreq * bw_ppm
    return xscale, spectral_bw_hz, bw_ppm


def iter_repetitions(header: mrd.Header,
                     input: Iterable[mrd.Acquisition],
                     line_broadening: float,
                     *,
                     pad: Optional[int] = None,
                     zero_lead: Optional[int] = None) -> Iterable[Tuple[mrd.Acquisition, np.ndarray, list]]:
    """
    Assemble EPSI k-space one repetition at a time.

    Yields (reference acquisition, kspace, the acquisitions of that repetition), where kspace
    has shape (nviews, nsamples, nechoes * FIDPAD). Holding one repetition's acquisitions lets
    the caller pass them through to its own output without a second read.
    Args:
        - pad: sampling window pad, or None to derive it from the recorded ramp time
        - zero_lead: leading samples to read as zero, or None for the sequence's own default
    """
    enc = header.encoding[0]
    nviews = 1
    if enc.encoding_limits.phase is not None:
        nviews = enc.encoding_limits.phase.maximum + 1

    if zero_lead is None:
        # the 13 unusable points are a quirk of this one sequence, unlike the ramp, which every
        # EPSI readout has and which the pad already accounts for
        epsigre = header.measurement_information.sequence_name == "epsigre"
        zero_lead = EPSIGRE_ZERO_LEAD if epsigre else 0
    leading_pad = None                  # derived from the first acquisition, which carries the dwell

    current_rep = None
    kspace = None
    reference_acq = None
    acqs: list = []

    for acq in input:
        nechoes, totalppswitch, nsamples = switch_layout(acq)
        if leading_pad is None:
            leading_pad, why = epsi_leading_pad(header, acq, pad)
            window = acq.head.discard_pre - leading_pad
            print(f"EPSI sampling window: {why}, so each echo reads positions "
                  f"{window}..{window + nsamples - 1} of its {totalppswitch} point switch"
                  + (f", with the first {zero_lead} samples of the readout read as zero"
                     if zero_lead else ""), file=sys.stderr)
        if acq.head.idx.repetition != current_rep:
            if kspace is not None:
                yield reference_acq, kspace, acqs
            kspace = np.zeros((nviews, nsamples, nechoes * FIDPAD), dtype='complex')
            reference_acq = acq
            acqs = []
            current_rep = acq.head.idx.repetition
        view = acq.head.idx.kspace_encode_step_1 if acq.head.idx.kspace_encode_step_1 is not None else 0
        kspace[view, :, :nechoes] = apply_line_broadening(acq, line_broadening,
                                                          leading_pad=leading_pad,
                                                          zero_lead=zero_lead)
        acqs.append(acq)

    if kspace is not None:
        yield reference_acq, kspace, acqs


def denoise_kspace(kspace: np.ndarray, rank: Optional[int]):
    """
    Truncated-SVD denoise an EPSI k-space cube along its spectral dimension.

    The cube is flattened to (nechoes, nviews * nsamples) so every column is the spectral FID
    of one k-space location, which is the matrix layout the low-rank argument applies to.
    Returns:
        - (denoised cube of the original shape, the DenoiseResult)
    """
    nviews, nsamples, nechoes = kspace.shape
    matrix = kspace.reshape(nviews * nsamples, nechoes).T
    result = denoise_svd(matrix, rank)
    return result.matrix.T.reshape(nviews, nsamples, nechoes), result


def reconstruct_volume(kspace: np.ndarray) -> np.ndarray:
    """FFT an EPSI k-space cube over all three axes into (views, readout, frequency)."""
    axes = tuple(range(kspace.ndim))
    return np.fft.fftshift(np.fft.fftn(kspace, axes=axes), axes=axes)


# ---------- EPSI analysis ------------------------------------------------


def phase_align(volumes: np.ndarray,
                reference_spect: np.ndarray,
                *,
                noise_threshold: Optional[float] = None,
                search_range: int = PHASE_SEARCH_RANGE) -> Tuple[np.ndarray, np.ndarray]:
    """
    Roll and rotate every voxel spectrum to maximise its overlap with a reference spectrum.

    Each voxel is free to be shifted by a few spectral points and rotated by a global phase;
    picking the combination that maximises overlap with the brightest spectrum in the series
    puts every voxel into a common frame, so the sum over voxels adds coherently instead of
    cancelling. Voxels below the noise threshold are left alone and excluded from the sum.
    Args:
        - volumes: (nreps, nviews, nro, nfreq) complex
        - reference_spect: (nfreq,) complex, the brightest voxel spectrum in the series
    Returns:
        - (aligned copy of volumes, the global spectrum summed over aligned voxels)
    """
    out = np.array(volumes, dtype=complex, copy=True)
    global_spect = np.zeros(out.shape[-1], dtype=complex)
    reference_conj = np.conj(reference_spect)

    for ide in range(out.shape[0]):
        for j in range(out.shape[1]):
            for k in range(out.shape[2]):
                spect = out[ide, j, k, :]
                if noise_threshold is not None and np.max(np.abs(spect)) < noise_threshold:
                    continue
                best_overlap = -np.inf
                best_r = 0
                best_th = 0.0
                for r in range(-search_range, search_range + 1):
                    rolled = np.roll(spect, r)
                    S0 = float(np.sum(np.real(rolled * reference_conj)))
                    Spi2 = float(np.sum(np.real(rolled * 1j * reference_conj)))
                    overlap = S0 * S0 + Spi2 * Spi2
                    if overlap > best_overlap:
                        best_overlap = overlap
                        best_r = r
                        best_th = np.pi / 2 - np.arctan2(S0, Spi2)
                corrected = np.roll(spect, best_r) * np.exp(1j * best_th)
                out[ide, j, k, :] = corrected
                global_spect += corrected

    return out, global_spect


def fit_global_multipeak(global_spect: np.ndarray,
                         xscale: np.ndarray,
                         spec: PeakSpec) -> Tuple[LorentzianFitter, int, float]:
    """
    Fit the summed spectrum, trying each candidate for which peak is the tallest one.

    The peak offsets are a rigid pattern whose absolute position is unknown, so each peak that
    is not marked _t is tried as the one sitting under the tallest point and the hypothesis
    with the smallest residual wins. The resulting centers, widths and phases then lock the
    line shape for the per-voxel fits.
    Returns:
        - (fitter holding the winning fit, index of the winning peak, the width guess in ppm)
    """
    scaling = float(np.max(np.abs(global_spect)))
    if scaling == 0.0:
        raise ValueError("global spectrum is zero; cannot fit peaks")
    norm = global_spect / scaling
    bw_ppm = float(xscale[-1] - xscale[0] + (xscale[1] - xscale[0]))
    width_guess = estimate_width_fwhm(xscale, norm)
    widths_init = np.full(len(spec), width_guess)
    width_bounds = (width_guess / 2, width_guess * 1.5)

    candidates = spec.biggest_idx or list(range(len(spec)))
    best_fitter = None
    best_idx = candidates[0]
    for icg in candidates:
        fitter = LorentzianFitter(xscale, spec.wigglefactor)
        centers = candidate_centers(xscale, norm, spec.offsets, icg, bw_ppm)
        params = fitter.fit_global(norm, centers, widths_init, width_bounds=width_bounds)
        if best_fitter is None or params.loss < best_fitter.params.loss:
            best_fitter = fitter
            best_idx = icg

    return best_fitter, best_idx, width_guess


def fit_voxel_amplitudes(volumes: np.ndarray,
                         fitter: LorentzianFitter,
                         *,
                         noise_threshold: float = 0.0,
                         skip_initial_reps: int = 0) -> np.ndarray:
    """
    Fit peak amplitudes in every voxel with the line shape locked by the global fit.

    Only the amplitudes and a small baseline are free, so a voxel too noisy to support a full
    fit still yields a usable amplitude.
    Args:
        - volumes: (nreps, nviews, nro, nfreq) complex, already phase aligned
    Returns:
        - (npeaks, nreps, nviews, nro) float, peak heights in the units of the input
    """
    npeaks = len(fitter.params.centers)
    nreps, ny, nx, _ = volumes.shape
    amplitudes = np.zeros((npeaks, nreps, ny, nx))

    for ide in range(skip_initial_reps, nreps):
        print(f"Fitting voxels for repetition={ide}", file=sys.stderr)
        for j in range(ny):
            for k in range(nx):
                spect = volumes[ide, j, k, :]
                if np.max(np.abs(spect)) < noise_threshold:
                    continue
                amplitudes[:, ide, j, k] = np.abs(fitter.fit_amplitudes(spect))

    return amplitudes


def fit_phantom(volumes: np.ndarray,
                xscale: np.ndarray,
                wigglefactor: float) -> Tuple[np.ndarray, float]:
    """
    Single-peak fit of a phantom acquisition, for the scaling it provides.

    A phantom is one tube of a known concentration, so a single Lorentzian per voxel is the
    whole model and the summed spectrum needs no phase alignment. The scaling is the amplitude
    times width, i.e. the peak area, of the summed spectrum, averaged over repetitions.
    Args:
        - volumes: (nreps, nviews, nro, nfreq) complex
    Returns:
        - (per-voxel peak area map of shape (nreps, nviews, nro), the averaged scaling)
    """
    nreps, ny, nx, _ = volumes.shape
    phantom_map = np.zeros((nreps, ny, nx))
    scaling_total = 0.0

    for i in range(nreps):
        summed = volumes[i].reshape(-1, volumes.shape[-1]).sum(axis=0)
        scale_g = float(np.max(np.abs(summed)))
        if scale_g == 0.0:
            continue
        norm_g = summed / scale_g
        width = estimate_width_fwhm(xscale, norm_g)
        width_bounds = (width / 2, width * 1.5)
        center = np.array([xscale[int(np.argmax(np.abs(norm_g)))]])
        params = LorentzianFitter(xscale, wigglefactor).fit_global(
            norm_g, center, np.array([width]), width_bounds=width_bounds)
        scaling_total += abs(params.amplitudes[0] * params.widths[0]) * scale_g / nreps

        for j in range(ny):
            for k in range(nx):
                spect = volumes[i, j, k, :]
                scale_v = float(np.max(np.abs(spect)))
                if scale_v == 0.0:
                    continue
                voxel = LorentzianFitter(xscale, wigglefactor).fit_global(
                    spect / scale_v, center, np.array([width]), width_bounds=width_bounds)
                phantom_map[i, j, k] = abs(voxel.amplitudes[0] * voxel.widths[0]) * scale_v

    return phantom_map, (scaling_total or 1.0)


def reconstruct_phantom(path: Path, line_broadening: float, wigglefactor: float,
                        *, pad: Optional[int] = None,
                        zero_lead: Optional[int] = None) -> Tuple[np.ndarray, float]:
    """
    Reconstruct a separately converted phantom scan and fit it.

    The phantom is named explicitly on the command line rather than detected in the stream: the
    converter deliberately does not mark acquisitions as phantom data, since whether a scan is
    calibration is a question about the scan and not about its acquisitions.
    """
    with open(path, "rb") as stream:
        with mrd.BinaryMrdReader(stream) as reader:
            header = reader.read_header()
            volumes = []
            xscale = None
            for reference_acq, kspace, _ in iter_repetitions(
                    header, acquisition_reader(reader.read_data()), line_broadening,
                    pad=pad, zero_lead=zero_lead):
                if xscale is None:
                    xscale, _, _ = spectral_axis(header, reference_acq)
                volumes.append(reconstruct_volume(kspace))

    if not volumes:
        raise ValueError(f"phantom file {path} contained no acquisitions")
    phantom_map, phantom_scaling = fit_phantom(np.stack(volumes), xscale, wigglefactor)
    print(f"Phantom scaling = {phantom_scaling}", file=sys.stderr)
    return phantom_map, phantom_scaling


# ---------- EPSI reconstruction ------------------------------------------


def reconstruct_epsi(header: mrd.Header,
                     input: Iterable[mrd.Acquisition],
                     *,
                     line_broadening: float,
                     spec: PeakSpec,
                     denoise: bool = False,
                     rank: int = None,
                     skip_initial_reps: int = 0,
                     phantom_scaling: float = 1.0,
                     phantom_map: np.ndarray = None,
                     pad: Optional[int] = None,
                     zero_lead: Optional[int] = None) -> Iterable[mrd.StreamItem]:
    """
    Reconstruct an EPSI acquisition into spectra, a global peak fit and metabolite maps.

    The reconstruction is in two parts. Each repetition is assembled, optionally denoised,
    transformed and emitted as it arrives, so a consumer sees repetition n before n+1 is read.
    The analysis then runs once over the whole series, because aligning voxels against the
    brightest spectrum and fitting a summed spectrum both need every repetition in hand.
    """
    volumes: List[np.ndarray] = []
    xscale = None
    spectral_bw_hz = None
    current_max = -np.inf
    max_spect = None
    max_location = (0, 0, 0)
    last_repetition = 0

    for irep, (reference_acq, kspace, acqs) in enumerate(
            iter_repetitions(header, input, line_broadening, pad=pad, zero_lead=zero_lead)):
        for acq in acqs:
            yield mrd.StreamItem.Acquisition(acq)

        if xscale is None:
            xscale, spectral_bw_hz, _ = spectral_axis(header, reference_acq)

        if denoise:
            kspace, denoise_result = denoise_kspace(kspace, rank)
            yield emit(denoise_result.singular_values,
                       dimension_labels=[mrd.ArrayDimension.SAMPLES],
                       description="singular_values",
                       svd_rank=denoise_result.rank,
                       repetition=irep)

        print(f"Reconstructing repetition={reference_acq.head.idx.repetition}", file=sys.stderr)
        img = reconstruct_volume(kspace)
        volumes.append(img)
        last_repetition = reference_acq.head.idx.repetition

        # the brightest voxel of the series is the reference every other voxel is aligned to
        peak_per_voxel = np.abs(img).max(axis=-1)
        j, k = np.unravel_index(int(np.argmax(peak_per_voxel)), peak_per_voxel.shape)
        if peak_per_voxel[j, k] > current_max:
            current_max = float(peak_per_voxel[j, k])
            max_spect = np.copy(img[j, k, :])
            max_location = (irep, int(j), int(k))

        yield emit(img,
                   head=mrd.NdArrayHeader(
                       dimension_labels=[mrd.ArrayDimension.Y,
                                         mrd.ArrayDimension.X,
                                         mrd.ArrayDimension.CONTRAST],
                       # no ArrayType means "a reconstructed spectral image", so it goes in the
                       # generic bucket and is identified by its image type and description
                       array_type=mrd.ArrayType.USER_MAP,
                       image_type=mrd.ArrayImageType.COMPLEX,
                       measurement_uid=reference_acq.head.measurement_uid,
                       average=reference_acq.head.idx.average,
                       repetition=reference_acq.head.idx.repetition,
                       acquisition_time_stamp_ns=reference_acq.head.acquisition_time_stamp_ns),
                   description="epsi_image",
                   xscale_ppm=xscale,
                   **{"receiver bandwidth(Hz)": spectral_bw_hz,
                      "line broadening(Hz)": line_broadening})

    if not volumes:
        return

    volumes = np.stack(volumes)

    # the last repetition of a hyperpolarized series has decayed away, so it measures noise
    noise = float(np.mean(np.abs(volumes[-1])))
    noise_threshold = noise * NOISE_THRESHOLD_MULTIPLIER
    yield emit(np.array([noise]),
               dimension_labels=[mrd.ArrayDimension.SAMPLES],
               array_type=mrd.ArrayType.NOISE,
               description="noise",
               noise_threshold_multiplier=NOISE_THRESHOLD_MULTIPLIER,
               estimated_from_repetition=last_repetition)

    yield emit(max_spect,
               dimension_labels=[mrd.ArrayDimension.CONTRAST],
               description="max spectral",
               xscale_ppm=xscale,
               max_repetition=max_location[0],
               max_y=max_location[1],
               max_x=max_location[2])

    if phantom_map is not None:
        yield emit(phantom_map,
                   dimension_labels=[mrd.ArrayDimension.REPETITION,
                                     mrd.ArrayDimension.Y,
                                     mrd.ArrayDimension.X],
                   array_type=mrd.ArrayType.PHANTOM,
                   description="phantom peak area",
                   phantom_scaling=phantom_scaling)

    print("Aligning voxel spectra", file=sys.stderr)
    aligned, global_spect = phase_align(volumes, max_spect, noise_threshold=noise_threshold)

    if len(spec) == 0:
        print("No peaks specified, skipping the Lorentzian fits", file=sys.stderr)
        yield emit(global_spect,
                   dimension_labels=[mrd.ArrayDimension.CONTRAST],
                   description="global_spect",
                   xscale_ppm=xscale)
        return

    print(f"Fitting {len(spec)} peaks to the global spectrum", file=sys.stderr)
    fitter, biggest_idx, width_guess = fit_global_multipeak(global_spect, xscale, spec)
    params = fitter.params
    global_scaling = float(np.max(np.abs(global_spect)))

    yield emit(global_spect,
               dimension_labels=[mrd.ArrayDimension.CONTRAST],
               description="global_spect",
               xscale_ppm=xscale,
               width_guess_ppm=width_guess,
               biggest_peak_index=biggest_idx,
               biggest_peak_name=spec.names[biggest_idx])
    yield emit(fitter.eval() * global_scaling,
               dimension_labels=[mrd.ArrayDimension.CONTRAST],
               description="global_spect_fit",
               xscale_ppm=xscale,
               fit_loss=params.loss)

    for description, values in (("lorentzian_centers_ppm", params.centers),
                                ("lorentzian_widths_ppm", params.widths),
                                ("lorentzian_phases_rad", params.phases),
                                ("lorentzian_amplitudes", params.amplitudes * global_scaling)):
        yield emit(values,
                   dimension_labels=[mrd.ArrayDimension.BASIS],
                   description=description,
                   peak_names=spec.names)
    yield emit(np.array([params.baseline * global_scaling]),
               dimension_labels=[mrd.ArrayDimension.SAMPLES],
               description="lorentzian_baseline")

    metabolites = fit_voxel_amplitudes(aligned, fitter,
                                       noise_threshold=noise_threshold,
                                       skip_initial_reps=skip_initial_reps)
    # peak area rather than peak height. The voxel fit locks the line shape, so the width is
    # the globally fitted one and this is a per-peak rescaling of the amplitudes
    areas = metabolites * np.abs(params.widths)[:, None, None, None]

    map_meta = dict(peak_names=spec.names,
                    peak_offsets_ppm=spec.offsets,
                    source_peak_index=spec.source_idx,
                    metabolite_indices=spec.metabolite_idx or None,
                    phantom_scaling=phantom_scaling)
    for description, values in (("metabolite_amplitude", metabolites),
                                ("metabolite_area", areas)):
        yield emit(values,
                   dimension_labels=[mrd.ArrayDimension.BASIS,
                                     mrd.ArrayDimension.REPETITION,
                                     mrd.ArrayDimension.Y,
                                     mrd.ArrayDimension.X],
                   description=description,
                   **map_meta)


# ---------- spectral (single voxel FID) reconstruction --------------------


def extract_fid(acq: mrd.Acquisition) -> np.ndarray:
    """Extract a single 1-D complex FID from an acquisition, trimming discard points."""
    fid = np.squeeze(np.asarray(acq.data)).astype(np.complex128).ravel()
    pre = acq.head.discard_pre or 0
    post = acq.head.discard_post or 0
    if pre or post:
        fid = fid[pre:len(fid) - post]
    return fid


def reconstruct_spectral(header: mrd.Header,
                         input: Iterable[mrd.Acquisition],
                         *,
                         line_broadening: float,
                         denoise: bool = False,
                         rank: int = None) -> Iterable[mrd.StreamItem]:
    """
    Reconstruct a single-voxel FID spectral time series.

    Each acquisition is one complex FID acquired at a time point of the metabolic flux
    experiment. The FIDs are stacked as the columns of an (nsamples, nspectra) matrix,
    optionally denoised by low rank approximation, then line broadened and Fourier
    transformed. Denoising acts on the raw complex signal, so no phase correction is needed
    beforehand.

    Yields the raw acquisitions, one complex spectrum per time point, and, when denoising, the
    singular values so the knee can be inspected and an explicit --rank chosen.
    """
    acqs = list(input)
    if not acqs:
        return
    fids = [extract_fid(acq) for acq in acqs]
    nsamples = min(f.shape[0] for f in fids)
    # stack as (nsamples, nspectra); truncate to common length for safety
    matrix = np.stack([f[:nsamples] for f in fids], axis=1)

    singular_values = None
    used_rank = None
    if denoise:
        result = denoise_svd(matrix, rank)
        matrix, singular_values, used_rank = result.matrix, result.singular_values, result.rank
        print(f"SVD denoising: matrix={matrix.shape}, rank={used_rank}, "
              f"top singular values={np.round(singular_values[:min(10, len(singular_values))], 4)}",
              file=sys.stderr)

    # frequency / ppm axis from the dwell time
    dwell_s = acqs[0].head.sample_time_ns / 1e9
    freq_hz = np.fft.fftshift(np.fft.fftfreq(nsamples, d=dwell_s))
    carrier_mhz = (header.experimental_conditions.h1resonance_frequency_hz or 0) / 1e6
    xaxis = freq_hz / carrier_mhz if carrier_mhz > 0 else freq_hz   # ppm, else Hz

    # line-broadening apodization applied to the FIDs before FFT
    t = np.arange(nsamples) * dwell_s
    apod = np.exp(-np.pi * line_broadening * t)

    for acq in acqs:
        yield mrd.StreamItem.Acquisition(acq)

    if singular_values is not None:
        yield emit(singular_values,
                   dimension_labels=[mrd.ArrayDimension.SAMPLES],
                   description="singular_values",
                   svd_rank=used_rank)

    for j, acq in enumerate(acqs):
        spectrum = np.fft.fftshift(np.fft.fft(matrix[:, j] * apod))
        yield emit(spectrum,
                   dimension_labels=[mrd.ArrayDimension.FREQUENCY],
                   description="global_spect",
                   xscale_ppm=xaxis,
                   repetition=j,
                   acquisition_time_stamp_ns=int(acq.head.acquisition_time_stamp_ns or 0),
                   **{"line broadening(Hz)": line_broadening})


# ---------- driver -------------------------------------------------------


def acquisition_reader(input: Iterable[mrd.StreamItem]) -> Iterable[mrd.Acquisition]:
    """Yield just the acquisitions out of an mrd stream."""
    for item in input:
        if isinstance(item, mrd.StreamItem.Acquisition):
            yield item.value


def reconstruct_mrs(input: BinaryIO,
                    output: BinaryIO,
                    *,
                    line_broadening: float,
                    spec: PeakSpec,
                    denoise: bool = False,
                    rank: int = None,
                    skip_initial_reps: int = 0,
                    phantom_scaling: float = 1.0,
                    phantom_map: np.ndarray = None,
                    pad: Optional[int] = None,
                    zero_lead: Optional[int] = None) -> None:
    """Reconstruct one converted file, choosing the reconstruction from its sequence."""
    with mrd.BinaryMrdReader(input) as reader:
        with mrd.BinaryMrdWriter(output) as writer:
            header = reader.read_header()
            family = sequence_family(header.measurement_information.sequence_name)
            append_recon_header(header,
                                line_broadening=line_broadening,
                                spec=spec,
                                denoise=denoise,
                                rank=rank,
                                phantom_scaling=phantom_scaling)
            writer.write_header(header)
            acquisitions = acquisition_reader(reader.read_data())
            if family == "epsi":
                writer.write_data(
                    reconstruct_epsi(header, acquisitions,
                                     line_broadening=line_broadening,
                                     spec=spec,
                                     denoise=denoise,
                                     rank=rank,
                                     skip_initial_reps=skip_initial_reps,
                                     phantom_scaling=phantom_scaling,
                                     phantom_map=phantom_map,
                                     pad=pad,
                                     zero_lead=zero_lead))
            else:
                if family != "spectral":
                    print(f"Unrecognized sequence "
                          f"'{header.measurement_information.sequence_name}', "
                          f"reconstructing it as a single voxel FID series", file=sys.stderr)
                writer.write_data(
                    reconstruct_spectral(header, acquisitions,
                                         line_broadening=line_broadening,
                                         denoise=denoise,
                                         rank=rank))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Reconstruct MRS data from an mrd2 file. Peaks are named as "
                    "-<name>[_smt] <ppm>, e.g. -pyr_s 9.7 -lac_m 21.8")
    parser.add_argument("-f", "--folder", type=Path, required=False, help="Folder to search for converted .mrd2 files, ignoring _recon.mrd2 output")
    parser.add_argument("-i", "--input", type=Path, required=False, help="Input mrd2 file")
    parser.add_argument("-o", "--output", type=Path, required=False, help="Output mrd2 file")
    parser.add_argument("-lb", "--line-broadening", type=float, default=42, required=False, help="Line broadening factor in Hz")
    parser.add_argument("-d", "--denoise", action="store_true", help="Apply truncated-SVD denoising before the transform")
    parser.add_argument("-r", "--rank", type=int, default=None, required=False, help="Number of singular values to retain (default: auto via Gavish-Donoho)")
    parser.add_argument("-w", "--wigglefactor", type=float, default=1.0, required=False, help="How far peak centers may move from their specified offsets, in ppm")
    parser.add_argument("--phantom", type=Path, default=None, required=False, help="A separately converted phantom .mrd2 to fit for the scaling it provides")
    parser.add_argument("--skip-initial-reps", type=int, default=0, required=False, help="Repetitions to leave out of the per-voxel fits")
    parser.add_argument("--pad", type=int, default=None, required=False, help="EPSI sampling window, as samples back from discard_pre. Default: derived from the recorded ramp time, use MRSreader.py -w to check it against the data")
    parser.add_argument("--zero-lead", type=int, default=None, required=False, help="Leading samples of each readout to read as zero, overriding the sequence's own default")
    return parser


if __name__ == "__main__":
    parser = build_parser()
    # the peak arguments have to come out before argparse sees them, since a peak's value may be
    # negative and argparse would read that as another option
    reserved = {option for action in parser._actions for option in action.option_strings}
    spec, remaining = split_peak_args(sys.argv[1:], reserved)
    args = parser.parse_args(remaining)
    spec = PeakSpec(names=spec.names, offsets=spec.offsets, modifiers=spec.modifiers,
                    wigglefactor=args.wigglefactor)

    phantom_map = None
    phantom_scaling = 1.0
    if args.phantom is not None:
        if not args.phantom.is_file():
            raise ValueError(f"{args.phantom} is not a file")
        phantom_map, phantom_scaling = reconstruct_phantom(
            args.phantom, args.line_broadening, spec.wigglefactor,
            pad=args.pad, zero_lead=args.zero_lead)

    recon_kwargs = dict(line_broadening=args.line_broadening,
                        spec=spec,
                        denoise=args.denoise,
                        rank=args.rank,
                        skip_initial_reps=args.skip_initial_reps,
                        pad=args.pad,
                        zero_lead=args.zero_lead,
                        phantom_scaling=phantom_scaling,
                        phantom_map=phantom_map)

    if args.folder and args.input:
        raise ValueError("Cannot specify both --folder and --input")
    elif args.folder:
        if not args.folder.is_dir():
            raise ValueError(f"{args.folder} is not a directory")
        else:
            # MRStomrd2 names a converted scan <meas_id>_<sequence>.mrd2, so there is no fixed
            # filename to look for and one directory can hold several. Everything this run writes
            # ends in _recon.mrd2, and so does anything an earlier run left behind, which is what
            # separates the inputs from the outputs
            mrd2_filepaths: List[Path] = []
            for root, dirnames, filenames in os.walk(args.folder):
                for filename in sorted(filenames):
                    if filename.endswith(".mrd2") and not filename.endswith("_recon.mrd2"):
                        mrd2_filepaths.append(Path(os.path.join(root, filename)))
            if len(mrd2_filepaths) > 0:
                # for raw filepath, run reconstruct with parameters from cmd line arguments
                for i, input_filepath in enumerate(mrd2_filepaths):
                    # named after the input rather than its directory, so that two scans converted
                    # into one directory do not reconstruct over each other
                    recon_filepath = input_filepath.with_name(input_filepath.stem + "_recon.mrd2")
                    print(f"Reconstructing {i+1}/{len(mrd2_filepaths)} at: {recon_filepath}", file=sys.stderr)
                    with open(input_filepath, "rb") as input, open(recon_filepath, "wb") as output:
                        reconstruct_mrs(input, output, **recon_kwargs)
            else:
                raise ValueError(f"No mrd2 files found in {args.folder}")
    elif args.input:
        if not args.input.is_file():
            raise ValueError(f"{args.input} is not a file")
        if args.output is None:
            raise ValueError("--output must be specified with --input")
        if not args.output.parent.is_dir():
            raise ValueError(f"Output directory {args.output.parent} does not exist")
        with open(args.input, "rb") as input, open(args.output, "wb") as output:
            reconstruct_mrs(input, output, **recon_kwargs)
    else:
        raise ValueError("Either --folder or --input must be specified")
