"""
Reconstruct a converted EPSI .mrd2 file into lorentzian peak fits and metabolite maps.

EPSI only: anything else raises.

tyger args:
- python mrd2recon.py
    - -i
    - $(INPUT_PIPE)
    - -o
    - $(OUTPUT_PIPE)

With --folder, every .mrd2 that is not itself a _recon.mrd2 is reconstructed to
<name>_recon.mrd2 beside it. With --input, one file is reconstructed to --output.

    python mrd2recon.py -f {directory} {metabolite shift parameters}
    python mrd2recon.py -i raw.mrd2 -o recon.mrd2 \
        -bic_tm 0.0 -urea 2.3 -pyr_s 9.7 -ala_tm 15.2 -hyd_tm 18.1 -lac_m 21.8

    _s  source peak, the injected substrate
    _t  tiny peak, not a candidate for "which peak is the tallest one"
    _m  a derived metabolite

The fit runs twice: once on the summed spectrum, to settle the relative ppm shift of each
metabolite, then once per voxel with that line shape. How far a voxel may depart from it is
given by three windows, all zero by default, which is to say the line shape is pinned and
only the amplitudes vary:

    -df   how far a peak center may move from the global fit, in ppm
    -dw   how far a peak width may move from the global fit, in ppm
    -dph  how far a peak phase may move from the global fit, in radians

The global fit has its own two constraints, GLOBAL_CENTER_WINDOW_PPM and WIDTH_BOUND_SCALE,
set at the top of this module to what mrd2_recon_to_incorporate.py enforced.

An averaged phantom prescan, flagged IS_NOISE_MEASUREMENT by the converter, is a separate
encoding with its own matrix and its own spectral axis, so it is assembled into its own buffer
and fitted on its own. It runs through the same alignment and the same fitter as the series -
only the noise floor, the peak list and the output names differ - and yields phantom_* arrays
alongside the metabolite_* ones. Its phantom_global_peak_areas is the single number the
prescan exists to provide: it is recorded, never applied, because what makes two experiments
comparable is the consumer's decision.

Everything the reconstruction produces goes to the output stream as mrd NdArrays, not as the
mrd Image field, and the raw acquisitions are kept in the recon file unchanged.

Where this departs from mrd2_recon_to_incorporate.py, deliberately:

    - the readout geometry is read from each acquisition, not hardcoded at necho=64, nro=12
    - the first spectral point is kept; the legacy's `if iecho > 0` left it zero, discarding
      a real FID point
    - every repetition is fitted; the legacy started at 4, a debug shortcut it kept
    - phantom voxels are fitted with the line shape pinned, like the series. The legacy fitted
      them unbounded. phantom_global_peak_areas is unaffected, since it comes from the global
      fit; only the per-voxel phantom map differs
    - voxel amplitudes start from the spectrum rather than from zero, which is a starting
      point rather than a result, but this is not a convex fit
"""

import argparse
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, BinaryIO, Iterable, List, Optional, Sequence, Tuple

import numpy as np

import mrd

from lorentzian_fitter import LorentzianFitter, candidate_centers, estimate_width_fwhm

# a voxel spectrum is fitted only if above the estimated noise floor
NOISE_THRESHOLD_MULTIPLIER = 3.0
# how far a voxel spectrum may be rolled when aligning it against the reference spectrum
PHASE_SEARCH_RANGE = 15
# Spectral zero fill factor; 1 means the spectral axis is exactly one point per switch.
#
# 2 rather than 1 because at 1 the lines are narrower than the grid can represent. On
# cirrhrat_0_1 the spectral step is 0.261 ppm and four of six fitted widths came back between
# 0.155 and 0.217 - sub-pixel. The model is A/(1 + i*delta/w), whose magnitude falls off as
# A*w/delta, so an unresolvably narrow peak becomes a tall spike with a heavy 1/delta tail:
# urea, 20x the height of anything near it, put 19488 of tail at 5.21 ppm where the whole
# spectrum reads 37842. Alanine then fitted the residue at 5.77 instead of the 5.2-5.4 its
# signal actually sits at.
#
# Zero filling to 2 halves the step to 0.130 ppm, which those same widths resolve. Alanine
# moves to 5.37, no width is sub-pixel, nothing clamps, and the fit improves on both channels -
# on cirrhrat_0_1 real 0.073 -> 0.048 and imag 0.076 -> 0.036. Flooring the width at one
# spectral step fixes the placement too, but by clamping four of six widths and at a much worse
# residual, so this is the better of the two.
FIDPAD = 1

# The global fit's two constraints. Both exist because lorn was linearised: main's lorn.py
# parameterized through arctan, so `c = centers + arctan(x)/pi * wigglefactor` capped a center
# at +-0.5 ppm and `w = widths * (1 + arctan(x) * 1.8/pi)` held a width to 0.1 to 1.9 times the
# guess, with no bounds argument passed to minimize at all. lorn_to_incorporate.py replaced both
# with plain offsets, `c = centers + x0` and `w = x0`, which confines nothing - so
# mrd2_recon_to_incorporate.py had to supply bounds by hand. It supplied one for the width, at
# the narrower (0.5, 1.5), and none for the center, and that omission is why a tiny peak can
# slide onto a strong neighbour and be fitted as a second component of its line: the residual
# falls while both peaks' maps are ruined. On cirrhrat_43_1, unbounded, centers walked 2.6 to
# 3.0 ppm off their placement and two peaks ended 0.01 ppm apart.
#
# The center window is a (lo, hi) clip on each peak's own width rather than one fixed number,
# because a center is determined about as well as its line is narrow: on cirrhrat_43_1 the
# fitted widths run 0.133 to 0.871 ppm, and a single window is either too loose for the sharp
# peaks or too tight for the broad ones. The clip floors it so a very narrow peak is not pinned
# to nothing, and caps it so a very broad one cannot wander onto its neighbour.
#
# The width bound is absolute, as multiples of the width guess, and (0.1, 1.9) restores what
# the arctan allowed. At (0.5, 1.5) five of six of this data's true widths fall outside the
# bound, so the peaks sit on it instead of on their own linewidth.
#
# append_recon_header records whichever values a run used, so a recon file says which regime
# produced it.
CenterWindow = Optional[Tuple[float, float]]
GLOBAL_CENTER_WINDOW_PPM: CenterWindow = (0.1, 1.0)
WIDTH_BOUND_SCALE = (0.1, 1.9)


# ---------- peak specification -------------------------------------------


@dataclass(frozen=True)
class PeakSpec:
    """The metabolite peaks to fit, as named on the command line or recorded in a header."""
    names: List[str]
    offsets: np.ndarray     # chemical shifts in ppm, same order as names
    modifiers: List[str]    # the letters after the underscore, same order as names

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


# A phantom holds one known metabolite, so its spectrum is fitted with a single peak rather
# than a pattern of offsets. One peak at offset 0 reduces candidate_centers to the argmax of
# the spectrum, which is where the legacy placed it, and leaves the hypothesis loop exactly one
# candidate to try - so the phantom fit runs through the same code as the metabolite fit.
PHANTOM_SPEC = PeakSpec(names=["phantom"], offsets=np.zeros(1), modifiers=[""])


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
          also what keeps '-dw 0.1' from being read as a peak named 'dw'
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

    spec = PeakSpec(names=names,
                    offsets=np.asarray(offsets, dtype=float),
                    modifiers=modifiers)
    return spec, remaining



def append_recon_header(header: mrd.Header, *,
                        line_broadening: float,
                        spec: PeakSpec,
                        fit_df: float = 0.0,
                        fit_dw: float = 0.0,
                        fit_dph: float = 0.0,
                        global_df: CenterWindow = GLOBAL_CENTER_WINDOW_PPM,
                        width_scale: Tuple[float, float] = WIDTH_BOUND_SCALE) -> mrd.Header:
    """
    Record the additional parameters required to run on recon on existing header

    The peak list goes on as one double per peak, named '<name>_<modifiers>', so a downstream
    consumer can recover which peak each map belongs to.
    """
    if header.user_parameters is None:
        header.user_parameters = mrd.UserParametersType()
    header.user_parameters.user_parameter_double.append(
        mrd.UserParameterDoubleType(name="line_broadening_factor", value=float(line_broadening)))
    header.user_parameters.user_parameter_double.append(
        mrd.UserParameterDoubleType(name="frequency_window_ppm", value=float(fit_df)))
    header.user_parameters.user_parameter_double.append(
        mrd.UserParameterDoubleType(name="linewidth_window_ppm", value=float(fit_dw)))
    header.user_parameters.user_parameter_double.append(
        mrd.UserParameterDoubleType(name="phase_window_rad", value=float(fit_dph)))
    # the center window is a clip on each peak's own width, so it takes two numbers; inf for
    # both says the centers were left unbounded
    gdf_lo, gdf_hi = ((float("inf"), float("inf")) if global_df is None
                      else (global_df, global_df) if np.isscalar(global_df)
                      else (global_df[0], global_df[1]))
    header.user_parameters.user_parameter_double.append(
        mrd.UserParameterDoubleType(name="global_frequency_window_lo_ppm", value=float(gdf_lo)))
    header.user_parameters.user_parameter_double.append(
        mrd.UserParameterDoubleType(name="global_frequency_window_hi_ppm", value=float(gdf_hi)))
    header.user_parameters.user_parameter_double.append(
        mrd.UserParameterDoubleType(name="width_bound_lo", value=float(width_scale[0])))
    header.user_parameters.user_parameter_double.append(
        mrd.UserParameterDoubleType(name="width_bound_hi", value=float(width_scale[1])))
    for i, name in enumerate(spec.names):
        full = f"{name}_{spec.modifiers[i]}" if spec.modifiers[i] else name
        header.user_parameters.user_parameter_double.append(
            mrd.UserParameterDoubleType(name=full, value=float(spec.offsets[i])))
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


def emit(data: np.ndarray, *,
         head: mrd.NdArrayHeader = None,
         dimension_labels: List = None,
         array_type: mrd.ArrayType = mrd.ArrayType.USER_MAP,
         **meta: Any) -> mrd.StreamItem:
    """
    Wrap an array as a stream-ready NdArray.

    Real data is normalised to float64 and complex data to complex128, which is also what
    picks the StreamItem union arm: those are the only two the readers on the other side
    expect, so the cast and the choice of arm are the same decision. Meta goes on the
    NdArray, which is where the dev schema puts it, not on the header, and a None value is
    dropped rather than encoded.
    """
    array = np.asarray(data)
    complex_data = np.iscomplexobj(array)
    array = array.astype(np.complex128) if complex_data else array.astype(np.float64)
    if head is None:
        head = mrd.NdArrayHeader(dimension_labels=dimension_labels, array_type=array_type)
    arr = mrd.NdArray(head=head, data=array,
                      meta={key: meta_values(value)
                            for key, value in meta.items() if value is not None})
    return (mrd.StreamItem.NdArrayComplexDouble(arr) if complex_data
            else mrd.StreamItem.NdArrayDouble(arr))


# ---------- k-space assembly ---------------------------------------------


def apply_line_broadening(acq: mrd.Acquisition,
                          line_broadening: float,
                          *,
                          nswitches: int,
                          npoints_per_switch: int) -> np.ndarray:
    """
    Reshape channel 0 of one readout into (kept points, switches) and apodize it.

    The readout is laid out as (discard_pre, npoints_per_switch, discard_post) repeated once per
    switch, so cirrhrat_43_1's 1792 samples are 64 switches of 28 points with 12 kept. Only
    discard_pre and npoints_per_switch are needed to find that flat top; what follows it is
    whatever is left over, and EncodingBuffer.check has already confirmed the geometry is the
    same for every acquisition of this encoding.

    The apodization decays along the switch axis, because one spectral point is acquired per
    switch and so tk advances by a whole switch: that axis is the FID time axis. Which is why
    the exponent is applied after the transpose, where the switch axis is last. Moving the
    transpose would silently apodize along the readout instead.
    Returns a (npoints_per_switch, nswitches) complex array
    """
    totalppswitch = acq.samples() // nswitches
    # version1: discard pre and post
    kspace_per_switch = acq.data[0,:].reshape(nswitches,totalppswitch)[:,acq.head.discard_pre:acq.head.discard_pre+npoints_per_switch]
    # version2: do not discard and take full total points per switch to line broadening
    # kspace_per_switch = acq.data[0,:].reshape(nswitches,totalppswitch)
    tk = np.arange(nswitches) * acq.head.sample_time_ns * totalppswitch / 1.0e+9
    return kspace_per_switch.T * np.exp(-tk * line_broadening)


def header_nswitches(header: mrd.Header) -> int:
    """
    The switch count, which is shared by every matrix in a file and recorded once on the header.

    Separate from group_layout because the EPSI test in reconstruct_mrs runs against the whole
    stream, series and prescan together, where the other three numbers are meaningless and the
    homogeneity check group_layout applies would rightly reject the mix.
    Returns 0 when the header records none, which is how a non-EPSI file reads.
    """
    if header.user_parameters is None:
        return 0
    for item in header.user_parameters.user_parameter_long:
        if item.name == "nswitches":
            return int(item.value)
    return 0


@dataclass
class EncodingBuffer:
    """
    One encoding's k-space, with the geometry and the spectral axis that belong to it.

    The series and the averaged prescan are separate matrices - on cirrhrat_43_1 the series
    keeps 12 points of a switch across 8 views and the prescan 18 across 12 - so every number
    describing one of them is held here rather than in a variable the two take turns writing.
    That is what the legacy could not do. It kept xscale, the width guess and the lorn module
    globals in one place for both, which is why it needed an explicit restore after the phantom
    voxel loop clobbered them, and why its answers depended on processing the phantom first.
    """
    img: np.ndarray             # (nreps, ky, kx, t); k-space until fft_kspace_to_image runs
    xscale: np.ndarray          # this encoding's own ppm axis
    kept: int                   # flat top points of one switch, which is the readout axis
    samples: int                # the three the homogeneity check compares
    discard_pre: int
    discard_post: int
    is_prescan: bool

    def check(self, acq: mrd.Acquisition) -> None:
        """
        Raise unless this acquisition was acquired at the geometry the buffer was built from.

        Every acquisition of one encoding shares a readout geometry, so a minority layout is a
        corrupted header rather than a variant. Reading the geometry once and never looking
        again let a single mutated discard_post rewrite kept from 12 to 20 for all 216
        acquisitions of a series, silently.
        """
        layout = (acq.samples(), acq.head.discard_pre or 0, acq.head.discard_post or 0)
        if layout != (self.samples, self.discard_pre, self.discard_post):
            raise ValueError(
                f"encoding {acq.head.encoding_space_ref} was built at (samples, discard_pre, "
                f"discard_post)={(self.samples, self.discard_pre, self.discard_post)} but "
                f"carries an acquisition at {layout}")


def make_buffer(header: mrd.Header,
                acq: mrd.Acquisition,
                nswitches: int,
                center_freq_hz: float) -> EncodingBuffer:
    """
    Build the buffer one encoding fills, from the first acquisition that names it.

    The readout geometry comes from the acquisition rather than the header user parameters,
    because those describe the series only - generate_header is built from a rawdata file - and
    the prescan beside it keeps a different number of points out of a different switch period.
    Reading npoints_per_switch off the header handed the prescan the series' 12 point window
    out of its own 18, and the series' spectral bandwidth along with it.

    The view and repetition counts are the two a readout cannot supply, so they come from the
    encoding the converter wrote for this matrix and pointed the acquisition at.
    Raises:
        - ValueError if the encoding declares no matrix, or if the discards leave no readout
    """
    ref = acq.head.encoding_space_ref or 0
    if ref >= len(header.encoding) or header.encoding[ref].encoding_limits is None:
        raise ValueError(f"an acquisition points at encoding {ref}, which this header "
                         f"does not describe")
    limits = header.encoding[ref].encoding_limits
    if limits.phase is None or limits.repetition is None:
        raise ValueError(f"encoding {ref} declares no view or no repetition count, so there is "
                         f"no matrix to fill")

    totalppswitch = acq.samples() // nswitches
    discard_pre = acq.head.discard_pre or 0
    discard_post = acq.head.discard_post or 0
    kept = totalppswitch - discard_pre - discard_post
    if kept <= 0:
        raise ValueError(f"discard_pre={discard_pre} and discard_post={discard_post} leave "
                         f"{kept} points of a {totalppswitch} point switch")

    # one spectral point is acquired per switch, so the spectral dwell is a whole switch period
    spectral_bw_hz = 1.0e+9 / (acq.head.sample_time_ns * totalppswitch)
    bw_ppm = spectral_bw_hz / center_freq_hz * 1.0e+6
    nfreq = nswitches * FIDPAD
    # a prescan is not a repetition of the series and carries no repetition axis of its own; in
    # tar mode every prescan file is handed rep_idx=0, so two of them would collide here
    return EncodingBuffer(
        img=np.zeros((limits.repetition.maximum + 1,    # repetition
                      limits.phase.maximum + 1,         # ky
                      kept,                             # kx
                      nfreq),                           # t
                     dtype='complex'),
        xscale=np.arange(nfreq) / nfreq * bw_ppm,
        kept=kept,
        samples=acq.samples(),
        discard_pre=discard_pre,
        discard_post=discard_post,
        is_prescan=bool(acq.head.flags & mrd.AcquisitionFlags.IS_NOISE_MEASUREMENT))


def fft_kspace_to_image(kspace: np.ndarray) -> np.ndarray:
    """FFT one repetition's k-space cube over all three axes into (views, readout, frequency)."""
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
    Takes (nreps, nviews, nro, nfreq) against the brightest voxel spectrum in the series, and
    returns (an aligned copy, the global spectrum summed over aligned voxels).
    """
    out = np.array(volumes, dtype=complex, copy=True)
    global_spect = np.zeros(out.shape[-1], dtype=complex)
    reference_conj = np.conj(reference_spect)

    for ide in range(out.shape[0]):
        for j in range(out.shape[1]):
            for k in range(out.shape[2]):
                spect = out[ide, j, k, :]
                # a repetition the stream never delivered is a zero sla0, and stays one
                if not np.any(spect):
                    continue
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
                         spec: PeakSpec,
                         center_window: CenterWindow = GLOBAL_CENTER_WINDOW_PPM,
                         width_scale: Tuple[float, float] = WIDTH_BOUND_SCALE
                         ) -> Tuple[LorentzianFitter, int, float]:
    """
    Fit the summed spectrum, trying each candidate for which peak is the tallest one.

    The peak offsets are a rigid pattern whose absolute position is unknown, so each peak that
    is not marked _t is tried as the one sitting under the tallest point and the hypothesis
    with the smallest residual wins. The resulting centers, widths and phases then lock the
    line shape for the per-voxel fits.

    center_window is how far the fit may move a center from that rigid placement, in ppm. The
    offsets are known chemistry, so a center is nearly determined before the fit starts; left
    unbounded, a tiny peak slides onto a strong neighbour and is fitted as a second component
    of its line. On a 6 peak kidney series hyd_tm walked 1.2 ppm onto urea, took a third of its
    amplitude, and cut the residual doing it. At 0 the placement is final. The default 0.5 is
    what lorn.py enforced implicitly through `c = centers + arctan(x0)/pi * wigglefactor` at
    the wigglefactor of 1.0 the EPSI path always ran at.

    width_scale is the (lo, hi) width bound as multiples of the width guess. It is what stops
    the optimizer walking a width through zero, and it is a real constraint on the fit: too
    narrow a range and a peak sits on a bound instead of on its own linewidth, which makes the
    model line the wrong shape and its height wrong with it. The default (0.1, 1.9) is the
    range `w = w0 * (1 + arctan(x0) * 1.8 / pi)` allowed. The (0.5, 1.5) this replaces came
    from mrd2_recon_to_incorporate.py, which linearised the transforms and re-supplied a
    narrower bound; on cirrhrat_43_1 it clamped four of six peaks and cost 0.7 of residual.
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
    width_bounds = (width_guess * width_scale[0], width_guess * width_scale[1])

    candidates = spec.biggest_idx or list(range(len(spec)))
    best_fitter = None
    best_idx = candidates[0]
    for icg in candidates:
        fitter = LorentzianFitter(xscale)
        centers = candidate_centers(xscale, norm, spec.offsets, icg, bw_ppm)
        params = fitter.fit_global(norm, centers, widths_init, width_bounds=width_bounds,
                                   center_window=center_window)
        if best_fitter is None or params.loss < best_fitter.params.loss:
            best_fitter = fitter
            best_idx = icg

    return best_fitter, best_idx, width_guess


def fit_voxel_peaks(volumes: np.ndarray,
                    fitter: LorentzianFitter,
                    *,
                    noise_threshold: float = 0.0,
                    fit_df: float = 0.0,
                    fit_dw: float = 0.0,
                    fit_dph: float = 0.0) -> Tuple[np.ndarray, np.ndarray]:
    """
    Fit every voxel with its line shape held at, or near, the one the global fit settled on.

    The windows say how far a voxel may depart from that shape. At their default of zero the
    line shape is pinned and only the amplitudes and a baseline are free, so a voxel too noisy
    to support a full fit still yields a usable amplitude. Opening them lets a voxel whose
    shim, and so whose line, differs from the average of the slice fit its own. The windows are
    in ppm, ppm and radians. Returns ((npeaks, nreps, nviews, nro) peak heights, the matching
    peak areas, each from its own voxel's fitted width).
    """
    npeaks = len(fitter.params.centers)
    nreps, ny, nx, _ = volumes.shape
    amplitudes = np.zeros((npeaks, nreps, ny, nx))
    areas = np.zeros((npeaks, nreps, ny, nx))

    for ide in range(nreps):
        print(f"Fitting voxels for repetition={ide}", file=sys.stderr)
        for j in range(ny):
            for k in range(nx):
                spect = volumes[ide, j, k, :]
                if not np.any(spect):
                    continue
                if np.max(np.abs(spect)) < noise_threshold:
                    continue
                params = fitter.fit_windowed(spect, df=fit_df, dw=fit_dw, dph=fit_dph)
                amplitudes[:, ide, j, k] = np.abs(params.amplitudes)
                areas[:, ide, j, k] = np.abs(params.amplitudes * params.widths)

    return amplitudes, areas


def fit_and_emit_peaks(aligned: np.ndarray,
                       global_spect: np.ndarray,
                       xscale: np.ndarray,
                       spec: PeakSpec,
                       *,
                       noise_threshold: float,
                       fit_df: float = 0.0,
                       fit_dw: float = 0.0,
                       fit_dph: float = 0.0,
                       label: str = "metabolite",
                       global_df: CenterWindow = GLOBAL_CENTER_WINDOW_PPM,
                       width_scale: Tuple[float, float] = WIDTH_BOUND_SCALE
                       ) -> Iterable[mrd.StreamItem]:
    """
    Fit the peaks on an aligned series and emit everything that describes the fit.

    Everything downstream of the alignment: the global fit, the line shape it settled on, and
    the per-voxel maps. Takes the aligned (nreps, nviews, nro, nfreq) series and the sum over
    its voxels, which the line shape is fitted to. Nothing flows back to the caller, so a run
    with no peaks named stops after the summed spectrum.

    `label` names every array this emits, since the series and the prescan are fitted
    separately and share one output stream; without it a reader would see two `global_spect`
    arrays and no way to tell which encoding each belongs to.
    """
    if len(spec) == 0:
        print("No peaks specified, skipping the Lorentzian fits", file=sys.stderr)
        yield emit(global_spect,
                   dimension_labels=[mrd.ArrayDimension.CONTRAST],
                   description=f"{label}_global_spect",
                   encoding=label,
                   xscale_ppm=xscale)
        return

    print(f"Fitting {len(spec)} peaks to the {label} global spectrum", file=sys.stderr)
    fitter, biggest_idx, width_guess = fit_global_multipeak(global_spect, xscale, spec,
                                                            center_window=global_df,
                                                            width_scale=width_scale)
    params = fitter.params
    global_scaling = float(np.max(np.abs(global_spect)))

    yield emit(global_spect,
               dimension_labels=[mrd.ArrayDimension.CONTRAST],
               description=f"{label}_global_spect",
               encoding=label,
               xscale_ppm=xscale,
               width_guess_ppm=width_guess,
               biggest_peak_index=biggest_idx,
               biggest_peak_name=spec.names[biggest_idx])
    yield emit(fitter.eval() * global_scaling,
               dimension_labels=[mrd.ArrayDimension.CONTRAST],
               description=f"{label}_global_spect_fit",
               encoding=label,
               xscale_ppm=xscale,
               fit_loss=params.loss)

    # the area of each peak in the summed spectrum, amplitude times width. For a phantom, whose
    # spec is one peak, this single number is the scaling the prescan exists to provide: it is
    # what makes two experiments comparable. Recorded, never applied - the legacy carried it
    # alongside the maps the same way, and what to do with it is the consumer's decision
    for description, values in ((f"{label}_lorentzian_centers_ppm", params.centers),
                                (f"{label}_lorentzian_widths_ppm", params.widths),
                                (f"{label}_lorentzian_phases_rad", params.phases),
                                (f"{label}_lorentzian_amplitudes",
                                 params.amplitudes * global_scaling),
                                (f"{label}_global_peak_areas",
                                 np.abs(params.amplitudes * params.widths) * global_scaling)):
        yield emit(values,
                   dimension_labels=[mrd.ArrayDimension.BASIS],
                   description=description,
                   encoding=label,
                   peak_names=spec.names)
    yield emit(np.array([params.baseline * global_scaling]),
               dimension_labels=[mrd.ArrayDimension.SAMPLES],
               description=f"{label}_lorentzian_baseline",
               encoding=label)

    # peak height and peak area. The area is each voxel's own amplitude times its own fitted
    # width, so with the width window closed it is a per-peak rescaling of the amplitudes and
    # with it open it is a genuinely per-voxel integral
    metabolites, areas = fit_voxel_peaks(aligned, fitter,
                                         noise_threshold=noise_threshold,
                                         fit_df=fit_df, fit_dw=fit_dw, fit_dph=fit_dph)

    map_meta = dict(encoding=label,
                    peak_names=spec.names,
                    peak_offsets_ppm=spec.offsets,
                    source_peak_index=spec.source_idx,
                    metabolite_indices=spec.metabolite_idx or None,
                    fit_df_ppm=fit_df,
                    fit_dw_ppm=fit_dw,
                    fit_dph_rad=fit_dph)
    for description, values in ((f"{label}_amplitude", metabolites),
                                (f"{label}_area", areas)):
        yield emit(values,
                   dimension_labels=[mrd.ArrayDimension.BASIS,
                                     mrd.ArrayDimension.REPETITION,
                                     mrd.ArrayDimension.Y,
                                     mrd.ArrayDimension.X],
                   description=description,
                   **map_meta)


def stream_epsi_image(img: np.ndarray,
                      xscale: np.ndarray,
                      spec: PeakSpec,
                      *,
                      is_prescan: bool,
                      label: str,
                      fit_df: float = 0.0,
                      fit_dw: float = 0.0,
                      fit_dph: float = 0.0,
                      global_df: CenterWindow = GLOBAL_CENTER_WINDOW_PPM,
                      width_scale: Tuple[float, float] = WIDTH_BOUND_SCALE
                      ) -> Iterable[mrd.StreamItem]:
    """
    Align one encoding's images, fit them, and emit everything that describes the fit.

    The series and the phantom prescan run through here alike, which is what the legacy did
    too: one loop over `allimgsets` did the brightest voxel search, the roll and rotate
    alignment, the global spectrum accumulation and the FWHM width guess for both. Only three
    things branch on which one this is - the noise floor, the peak list, and what the outputs
    are called - and all three arrive as arguments, so neither can reach into the other's
    state the way the legacy's module level globals let them.
    Args:
        - img: this encoding's (nreps, ny, nx, nfreq) images, already transformed
        - xscale: this encoding's own ppm axis; a prescan's switch period differs from the
          series', so its spectral bandwidth does too
        - spec: the metabolite pattern, or PHANTOM_SPEC for a prescan
        - is_prescan: selects the noise floor, the one thing that cannot be read off the data
        - label: names the arrays this emits, so both encodings can share one stream
    """
    if is_prescan:
        # every voxel of a phantom is signal, which is why the legacy set a phantom's noise to
        # zero outright rather than gating it
        noise_threshold = 0.0
    else:
        # the last repetition of a hyperpolarized series has decayed away, so it measures noise
        noise = float(np.mean(np.abs(img[-1])))
        noise_threshold = noise * NOISE_THRESHOLD_MULTIPLIER
        if noise == 0.0:
            print(f"warning: repetition {img.shape[0] - 1} is empty, so the noise floor is 0 "
                  f"and every voxel above zero will be fitted", file=sys.stderr)

    # the brightest voxel is the reference every other voxel is aligned to. One argmax over the
    # whole series picks the same voxel the legacy's triple loop did: both take the first
    # maximum, and C order over this array is repetition order
    peak_per_voxel = np.abs(img).max(axis=-1)   # across t
    max_rep, max_y, max_x = (int(i) for i in
                             np.unravel_index(int(np.argmax(peak_per_voxel)),
                                              peak_per_voxel.shape))
    print(f"Aligning phase for {label}", file=sys.stderr)
    aligned, global_spect = phase_align(img, img[max_rep, max_y, max_x, :],
                                        noise_threshold=noise_threshold)

    yield emit(aligned,
               dimension_labels=[mrd.ArrayDimension.REPETITION,
                                 mrd.ArrayDimension.Y,
                                 mrd.ArrayDimension.X,
                                 mrd.ArrayDimension.CONTRAST],
               description=f"{label}_image_aligned",
               encoding=label,
               noise_threshold=noise_threshold,
               phase_search_range=PHASE_SEARCH_RANGE,
               reference_repetition=max_rep,
               reference_y=max_y,
               reference_x=max_x)

    yield from fit_and_emit_peaks(aligned, global_spect, xscale, spec,
                                  label=label,
                                  noise_threshold=noise_threshold,
                                  fit_df=fit_df, fit_dw=fit_dw, fit_dph=fit_dph,
                                  global_df=global_df, width_scale=width_scale)


def reconstruct_mrs(input: BinaryIO,
                    output: BinaryIO,
                    *,
                    line_broadening: float,
                    spec: PeakSpec,
                    fit_df: float = 0.0,
                    fit_dw: float = 0.0,
                    fit_dph: float = 0.0,
                    global_df: CenterWindow = GLOBAL_CENTER_WINDOW_PPM,
                    width_scale: Tuple[float, float] = WIDTH_BOUND_SCALE) -> None:
    """
    Reconstruct one converted EPSI file.

    Acquisitions are routed by encoding_space_ref, which the converter sets per matrix, so the
    series and the averaged prescan each fill their own buffer and neither can be assembled
    into the other's. Nothing reads the FIRST_IN_PHASE / LAST_IN_REPETITION flags: a cube is
    complete when the stream ends, not when a flag says so, which is what makes the result
    independent of the order the acquisitions arrive in and leaves a series missing a view with
    zero rows rather than everything shifted by one.

    The transform runs once per repetition on the whole (ky, kx, t) cube. Transforming each
    readout as it arrived would leave the phase encode axis in k-space, because a DFT along ky
    is a sum over every view and no single acquisition holds more than one of them.
    Raises:
        - ValueError if the header describes no EPSI matrix, records no resonance frequency,
          or if the stream holds no acquisitions
    """
    with mrd.BinaryMrdReader(input) as reader:
        header = reader.read_header()

        nswitches = header_nswitches(header)
        if not nswitches:
            raise ValueError("this header records no switch count, so it describes no EPSI "
                             "matrix")
        # the converter writes the 13C frequency here despite the field name, since that is the
        # frequency these spectra were actually acquired at
        center_freq_hz = header.experimental_conditions.h1resonance_frequency_hz
        if not center_freq_hz:
            raise ValueError("the header records no resonance frequency, so there is no ppm "
                             "axis")

        with mrd.BinaryMrdWriter(output) as writer:
            # record what this reconstruction was run with, so a recon file says how it was made
            header = append_recon_header(header,
                                         line_broadening=line_broadening,
                                         spec=spec,
                                         fit_df=fit_df,
                                         fit_dw=fit_dw,
                                         fit_dph=fit_dph,
                                         global_df=global_df,
                                         width_scale=width_scale)
            writer.write_header(header)

            buffers: dict = {}
            acquisitions: List[mrd.Acquisition] = []
            for item in reader.read_data():
                if not isinstance(item, mrd.StreamItem.Acquisition):
                    continue
                acq = item.value
                acquisitions.append(acq)
                ref = acq.head.encoding_space_ref or 0
                if ref not in buffers:
                    buffers[ref] = make_buffer(header, acq, nswitches, center_freq_hz)
                buffer = buffers[ref]
                buffer.check(acq)
                buffer.img[acq.head.idx.repetition or 0,
                           acq.head.idx.kspace_encode_step_1 or 0,
                           :,
                           :nswitches] = apply_line_broadening(
                               acq, line_broadening,
                               nswitches=nswitches,
                               npoints_per_switch=buffer.kept)
            if not acquisitions:
                raise ValueError("this stream holds no acquisitions, so there is nothing to "
                                 "reconstruct")

            # one transform per repetition, over all three axes of its cube
            print(f"Transforming {len(buffers)} encoding(s) from "
                  f"{len(acquisitions)} acquisitions", file=sys.stderr)
            for buffer in buffers.values():
                for rep in range(buffer.img.shape[0]):
                    buffer.img[rep] = fft_kspace_to_image(buffer.img[rep])

            for ref, buffer in sorted(buffers.items()):
                label = "phantom" if buffer.is_prescan else "metabolite"
                writer.write_data(stream_epsi_image(
                    buffer.img, buffer.xscale,
                    PHANTOM_SPEC if buffer.is_prescan else spec,
                    is_prescan=buffer.is_prescan,
                    label=label,
                    fit_df=fit_df, fit_dw=fit_dw, fit_dph=fit_dph,
                    global_df=global_df, width_scale=width_scale))

            # the raw acquisitions are kept in the recon file unchanged, as the legacy kept them
            writer.write_data(mrd.StreamItem.Acquisition(a) for a in acquisitions)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Reconstruct MRS data from an mrd2 file. Peaks are named as "
                    "-<name>[_smt] <ppm>, e.g. -pyr_s 9.7 -lac_m 21.8")
    parser.add_argument("-f", "--folder", type=Path, required=False, help="Folder to search for converted .mrd2 files, ignoring _recon.mrd2 output")
    parser.add_argument("-i", "--input", type=Path, required=False, help="Input mrd2 file")
    parser.add_argument("-o", "--output", type=Path, required=False, help="Output mrd2 file")
    parser.add_argument("-lb", "--line-broadening", type=float, default=42, required=False, help="Line broadening factor in Hz")
    parser.add_argument("-df", "--fit-df", type=float, default=0.0, required=False, 
                        help="How far a peak center may move from the global fit during the per-voxel fit, in ppm. Default 0, i.e. held at the global fit")
    parser.add_argument("-dw", "--fit-dw", type=float, default=0.0, required=False, 
                        help="How far a peak width may move from the global fit during the per-voxel fit, in ppm. Default 0, i.e. held at the global fit")
    parser.add_argument("-dph", "--fit-dph", type=float, default=0.0, required=False, 
                        help="How far a peak phase may move from the global fit during the per-voxel fit, in radians. Default 0, i.e. held at the global fit")

    parser.add_argument("--fidpad", type=int, default=FIDPAD, required=False,
                        help=f"Spectral zero fill factor. Default {FIDPAD}. At 1 the spectral "
                             f"step is one point per switch, which on these scans is wider than "
                             f"the lines, so a peak becomes a sub-pixel spike whose tail "
                             f"displaces its weaker neighbours")

    # the peak arguments have to come out before argparse sees them, since a peak's value may be
    # negative and argparse would read that as another option
    reserved = {option for action in parser._actions for option in action.option_strings}
    spec, remaining = split_peak_args(sys.argv[1:], reserved)
    args = parser.parse_args(remaining)

    # read as a module global by spectral_axis and the buffers, so set it before anything runs
    FIDPAD = args.fidpad

    recon_kwargs = dict(line_broadening=args.line_broadening,
                        spec=spec,
                        fit_df=args.fit_df,
                        fit_dw=args.fit_dw,
                        fit_dph=args.fit_dph)

    if args.folder and args.input:
        raise ValueError("Cannot specify both --folder (local only) and --input")
    elif args.folder:
        if not args.folder.is_dir():
            raise ValueError(f"{args.folder} is not a directory")
        else:
            mrd2_filepaths: List[Path] = []
            for root, dirnames, filenames in os.walk(args.folder):
                for filename in sorted(filenames):
                    if filename.endswith(".mrd2") and not filename.endswith("_recon.mrd2"):
                        mrd2_filepaths.append(Path(os.path.join(root, filename)))
            if len(mrd2_filepaths) > 0:
                for i, input_filepath in enumerate(mrd2_filepaths):
                    # generate output file with suffix '_recon.mrd2' of raw filename
                    recon_filepath = input_filepath.with_name(input_filepath.stem + "_recon.mrd2")
                    print(f"Reconstructing {i+1}/{len(mrd2_filepaths)} at: {recon_filepath}", file=sys.stderr)
                    with open(input_filepath, "rb") as input, open(recon_filepath, "wb") as output:
                        reconstruct_mrs(input, output, **recon_kwargs)
                    # except ValueError as err:
                    #     # a folder holds whatever was converted into it, and this reconstruction is
                    #     # EPSI only. One spectral scan among the inputs should not abandon the rest,
                    #     # so drop the empty output it left behind and carry on
                    #     recon_filepath.unlink(missing_ok=True)
                    #     print(f"  skipping {input_filepath.name}: {err}", file=sys.stderr)
            else:
                raise ValueError(f"No raw mrd2 files found in {args.folder}")
    elif args.input:
        # no is_file() check: under tyger --input is a named pipe, for which it is False
        if args.output is None:
            raise ValueError("--input needs --output")
        with open(args.input, "rb") as input, open(args.output, "wb") as output:
            reconstruct_mrs(input, output, **recon_kwargs)
    else:
        raise ValueError("Either --folder or --input must be specified")
