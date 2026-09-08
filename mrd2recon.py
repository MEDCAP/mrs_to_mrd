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

The global fit has its own two knobs, -gdf for how far it may move a center from the rigid
placement and --width-bounds for the linewidth range, both as multiples of the width guess.

An averaged phantom prescan, flagged IS_NOISE_MEASUREMENT by the converter, is reconstructed
on its own rather than as a repetition of the series. It contributes a phantom_image and a
phantom_scaling, the fitted area of its single peak, which is recorded on the metabolite maps
rather than applied to them.

Everything the reconstruction produces goes to the output stream as mrd NdArrays, not as the
mrd Image field, and the raw acquisitions are kept in the recon file unchanged.
"""

import argparse
import os
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any, BinaryIO, Iterable, Iterator, List, Optional, Sequence, Tuple

import numpy as np

import mrd

from lorentzian_fitter import LorentzianFitter, candidate_centers, estimate_width_fwhm

# a voxel spectrum is fitted only if above the estimated noise floor
NOISE_THRESHOLD_MULTIPLIER = 3.0
# how far a voxel spectrum may be rolled when aligning it against the reference spectrum
PHASE_SEARCH_RANGE = 15
# spectral zero fill factor; 1 means the spectral axis is exactly one point per echo
FIDPAD = 1


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
                        global_df: float = 0.5,
                        width_scale: Tuple[float, float] = (0.1, 1.9)) -> mrd.Header:
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
    header.user_parameters.user_parameter_double.append(
        mrd.UserParameterDoubleType(name="global_frequency_window_ppm", value=float(global_df)))
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
    whatever is left over, and group_layout has already checked the geometry is the same for
    every acquisition in the group.

    The apodization decays along the switch axis, because one spectral point is acquired per
    switch and so tk advances by a whole switch: that axis is the FID time axis. Which is why
    the exponent is applied after the transpose, where the switch axis is last. Moving the
    transpose would silently apodize along the readout instead.
    Returns a (npoints_per_switch, nswitches) complex array
    Raises:
        - ValueError if the readout is too short to hold the switches it should
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


@dataclass(frozen=True)
class Layout:
    """The acquisition matrix one group of acquisitions was acquired at."""
    nswitches: int
    totalppswitch: int      # the whole switch, flat top and discards together
    kept: int               # points of the flat top, which is the readout axis
    nviews: int
    nreps: int


def group_layout(header: mrd.Header, acqs: List[mrd.Acquisition]) -> Layout:
    """
    The acquisition matrix one group of acquisitions was acquired at.

    Three of the five numbers come from the readout, and are right for any matrix: nswitches is
    shared by every matrix in a file and is on the header as a user parameter, the whole switch
    is nsamples / nswitches, and the kept flat top is that less the discard counts. On
    cirrhrat_43_1 that gives 12 points a switch for the series and 18 for the averaged prescan
    beside it, which is correct for both.

    The view and repetition counts are the two a readout cannot supply, and the views are what
    differs between the matrices: 8 for that series, 12 for its prescan. The converter records a
    matrix per encoding and points each acquisition at its own through encoding_space_ref, so
    read them from there. Files converted before that leave the ref null and carry a single
    encoding describing the series only, so fall back to the highest index actually present -
    right whenever a group is complete, and an undercount on a scan missing a view.

    Everything here is a loop bound or an array dimension downstream, and none of it is
    checkable once the numbers have been handed on, so this is where the readout is validated.
    Raises rather than picking a majority, because a group that disagrees with itself means
    something upstream is wrong and a quietly reasonable answer would hide it.
    Raises:
        - ValueError if the header records no switch count, if the group mixes geometries, if it
          keeps no points, or if it carries a view or repetition index the encoding does not
          declare
    """
    nswitches = header_nswitches(header)
    if not nswitches:
        raise ValueError("this header records no switch count, so it describes no EPSI matrix")

    # every acquisition in one group shares a readout geometry, so a minority layout is a
    # corrupted header, not a variant. Taking it from acqs[0] let a single mutated
    # discard_post rewrite kept from 12 to 20 for all 216 acquisitions of a series, silently
    layouts = Counter((a.samples(), a.head.discard_pre or 0, a.head.discard_post or 0,
                       a.head.encoding_space_ref)
                      for a in acqs)
    if len(layouts) > 1:
        detail = "; ".join(f"{n} at (samples, discard_pre, discard_post, encoding_ref)={lay}"
                           for lay, n in layouts.most_common())
        raise ValueError(f"this group mixes {len(layouts)} readout geometries: {detail}")
    (nsamples, discard_pre, discard_post, ref), _ = layouts.most_common(1)[0]

    totalppswitch = nsamples // nswitches
    kept = totalppswitch - discard_pre - discard_post
    if kept <= 0:
        raise ValueError(f"discard_pre={discard_pre} and discard_post={discard_post} leave "
                         f"{kept} points of a {totalppswitch} point switch")

    views_present = {a.head.idx.kspace_encode_step_1 or 0 for a in acqs}
    reps_present = {a.head.idx.repetition or 0 for a in acqs}
    # every encoding limit is optional, so an encoding that declares neither reads like a file
    # converted before the per-encoding matrices and falls back the same way
    limits = (header.encoding[ref].encoding_limits
              if ref is not None and ref < len(header.encoding) else None)
    nviews = (limits.phase.maximum + 1
              if limits is not None and limits.phase is not None else max(views_present) + 1)
    nreps = (limits.repetition.maximum + 1
             if limits is not None and limits.repetition is not None
             else max(reps_present) + 1)
    if max(views_present) >= nviews:
        raise ValueError(f"encoding {ref} declares {nviews} views but view index "
                         f"{max(views_present)} is present")
    if max(reps_present) >= nreps:
        raise ValueError(f"encoding {ref} declares {nreps} repetitions but repetition index "
                         f"{max(reps_present)} is present")
    return Layout(nswitches=nswitches, totalppswitch=totalppswitch, kept=kept,
                  nviews=nviews, nreps=nreps)


def group_by_repetition(acqs: List[mrd.Acquisition]) -> dict:
    """
    The acquisitions of each repetition, keyed by repetition index and in stream order.

    The repetition is the unit the reconstruction is built from, so the grouping is explicit
    rather than a side effect of accumulating cubes. The first acquisition of each group is
    the one whose header describes that repetition.
    """
    by_rep: dict = {}
    for acq in acqs:
        by_rep.setdefault(acq.head.idx.repetition or 0, []).append(acq)
    return by_rep


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
                         center_window: float = 0.5,
                         width_scale: Tuple[float, float] = (0.1, 1.9)
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
                       global_df: float = 0.5,
                       width_scale: Tuple[float, float] = (0.1, 1.9)) -> Iterable[mrd.StreamItem]:
    """
    Fit the peaks on an aligned series and emit everything that describes the fit.

    Everything downstream of the alignment: the global fit, the line shape it settled on, and
    the per-voxel maps. Takes the aligned (nreps, nviews, nro, nfreq) series and the sum over
    its voxels, which the line shape is fitted to. Nothing flows back to the caller, so a run
    with no peaks named stops after the summed spectrum.
    """
    if len(spec) == 0:
        print("No peaks specified, skipping the Lorentzian fits", file=sys.stderr)
        yield emit(global_spect,
                   dimension_labels=[mrd.ArrayDimension.CONTRAST],
                   description="global_spect",
                   xscale_ppm=xscale)
        return

    print(f"Fitting {len(spec)} peaks to the global spectrum", file=sys.stderr)
    fitter, biggest_idx, width_guess = fit_global_multipeak(global_spect, xscale, spec,
                                                            center_window=global_df,
                                                            width_scale=width_scale)
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

    # peak height and peak area. The area is each voxel's own amplitude times its own fitted
    # width, so with the width window closed it is a per-peak rescaling of the amplitudes and
    # with it open it is a genuinely per-voxel integral
    metabolites, areas = fit_voxel_peaks(aligned, fitter,
                                         noise_threshold=noise_threshold,
                                         fit_df=fit_df, fit_dw=fit_dw, fit_dph=fit_dph)

    map_meta = dict(peak_names=spec.names,
                    peak_offsets_ppm=spec.offsets,
                    source_peak_index=spec.source_idx,
                    metabolite_indices=spec.metabolite_idx or None,
                    fit_df_ppm=fit_df,
                    fit_dw_ppm=fit_dw,
                    fit_dph_rad=fit_dph)
                    # derived from the phantom prescan and recorded, not applied: it is what
                    # makes two experiments comparable, and that is the consumer's decision
    for description, values in (("metabolite_amplitude", metabolites),
                                ("metabolite_area", areas)):
        yield emit(values,
                   dimension_labels=[mrd.ArrayDimension.BASIS,
                                     mrd.ArrayDimension.REPETITION,
                                     mrd.ArrayDimension.Y,
                                     mrd.ArrayDimension.X],
                   description=description,
                   **map_meta)


def reconstruct_phantom(header: mrd.Header,
                        phantom_acqs: List[mrd.Acquisition],
                        line_broadening: float,
                        width_scale: Tuple[float, float] = (0.1, 1.9)
                        ) -> Tuple[float, Optional[np.ndarray]]:
    """
    Reconstruct the urea phantom prescan and derive the scaling factor it exists to provide.

    The phantom holds one known metabolite, so its spectrum is fitted with a single peak rather
    than the metabolite pattern, and the quantity of interest is that peak's area, amplitude
    times width: that is what makes two experiments comparable. No noise gate applies here -
    every voxel of a phantom is signal, which is why the legacy sets a phantom's noise to zero
    outright. Several phantom sets average into one number.

    The prescan is assembled here rather than by the caller. Its matrix need not match the
    series', and on cirrhrat_43_1 it does not, so it carries its own layout and its own
    spectral axis and there is nothing for the caller to hold on its behalf.

    width_scale is the same (lo, hi) bound the metabolite fit uses, as multiples of the width
    guess. It matters more here than there. The scaling is the peak's area, amplitude times
    width, so a width held off its true value scales every experiment this number normalises.
    Returns:
        - (the scaling, 1.0 when there is no phantom to derive one from,
           the per-voxel area map of the last phantom set, or None)
    """
    layout = group_layout(header, phantom_acqs)
    nswitches, totalppswitch, kept, nviews = (layout.nswitches, layout.totalppswitch,
                                              layout.kept, layout.nviews)
    xscale, _ = spectral_axis(header, phantom_acqs[0], nswitches, totalppswitch)

    scalings: List[float] = []
    voxel_map = None

    for rep, rep_acqs in sorted(group_by_repetition(phantom_acqs).items()):
        kspace = np.zeros((nviews, kept, nswitches * FIDPAD), dtype='complex')
        for acq in rep_acqs:
            kspace[acq.head.idx.kspace_encode_step_1 or 0, :, :nswitches] = \
                apply_line_broadening(acq, line_broadening, 
                                      nswitches=nswitches,
                                      npoints_per_switch=kept)
        missing = sorted(set(range(nviews))
                         - {a.head.idx.kspace_encode_step_1 or 0 for a in rep_acqs})
        if missing:
            print(f"warning: prescan {rep} is missing views {missing}, which keep zero rows "
                  f"through the transform", file=sys.stderr)

        img = fft_kspace_to_image(kspace)
        # every voxel contributes; a phantom has no noise floor to gate against
        global_spect = img.sum(axis=(0, 1))
        scaling = float(np.max(np.abs(global_spect)))
        if scaling == 0.0:
            continue
        norm = global_spect / scaling

        width_guess = estimate_width_fwhm(xscale, norm)
        fitter = LorentzianFitter(xscale)
        center = np.array([xscale[int(np.argmax(np.abs(norm)))]])
        params = fitter.fit_global(norm, center, np.array([width_guess]),
                                   width_bounds=(width_guess * width_scale[0],
                                                 width_guess * width_scale[1]))
        scalings.append(float(np.abs(params.amplitudes[0] * params.widths[0]) * scaling))

        # the same area, per voxel, with the line shape held at the one the set settled on
        voxel_map = np.zeros(img.shape[:2])
        for j in range(img.shape[0]):
            for k in range(img.shape[1]):
                voxel = fitter.fit_windowed(img[j, k, :])
                voxel_map[j, k] = float(np.abs(voxel.amplitudes[0] * voxel.widths[0]))

    if not scalings:
        return 1.0, None
    return float(np.mean(scalings)), voxel_map


def epsi_images(acqs: List[mrd.Acquisition],
                line_broadening: float,
                layout: Layout) -> Iterator[Tuple[int, mrd.Acquisition, np.ndarray]]:
    """
    Transform one repetition's k-space cube into an image, a repetition at a time.

    The repetition is the smallest unit that can be transformed. One acquisition is one view, so
    the view axis is complete only once every acquisition of that repetition is in hand, and
    transforming all three axes of the cube is what turns it into (views, readout, frequency).
    The cube is local to the iteration; only the image outlives it.
    Yields (repetition index, that repetition's first acquisition, its image), in repetition
    order. The acquisition is the one whose header describes that repetition.
    """
    for rep, rep_acqs in sorted(group_by_repetition(acqs).items()):
        kspace = np.zeros((layout.nviews, layout.kept, layout.nswitches * FIDPAD),
                          dtype='complex')
        for acq in rep_acqs:
            kspace[acq.head.idx.kspace_encode_step_1 or 0, :, :layout.nswitches] = \
                apply_line_broadening(acq, line_broadening,
                                      nswitches=layout.nswitches,
                                      totalppswitch=layout.totalppswitch,
                                      npoints_per_switch=layout.kept)
        missing = sorted(set(range(layout.nviews))
                         - {a.head.idx.kspace_encode_step_1 or 0 for a in rep_acqs})
        if missing:
            print(f"warning: repetition {rep} is missing views {missing}, which keep zero rows "
                  f"through the transform", file=sys.stderr)
        yield rep, rep_acqs[0], fft_kspace_to_image(kspace)


def stream_epsi_image(header: mrd.Header,
                      img,
                      prescan_img,
                      *,
                      line_broadening: float,
                      spec: PeakSpec,
                      fit_df: float = 0.0,
                      fit_dw: float = 0.0,
                      fit_dph: float = 0.0,
                      global_df: float = 0.5,
                      width_scale: Tuple[float, float] = (0.1, 1.9),
                      xscale: np.ndarray
                      ) -> Iterable[mrd.StreamItem]:
    """
    FFT -> fitting
    """
    # the last repetition of a hyperpolarized series has decayed away, so it measures noise.
    noise = float(np.mean(np.abs(img[-1])))
    noise_threshold = noise * NOISE_THRESHOLD_MULTIPLIER
    if noise == 0.0:
        print(f"warning: repetition {noise_repetition} is empty, so the noise floor is 0 and "
              f"every voxel above zero will be fitted", file=sys.stderr)
    # the brightest voxel of the series is the reference every other voxel is aligned to. One
    # argmax over the series picks the same voxel the per-repetition scan it replaces did: both
    # take the first maximum, and C order over this array is repetition order
    peak_per_voxel = np.abs(img).max(axis=-1)   # across t
    max_rep, max_y, max_x = (int(i) for i in
                             np.unravel_index(int(np.argmax(peak_per_voxel)),
                                              peak_per_voxel.shape))
    print("Aligning phase", file=sys.stderr)
    aligned, global_spect = phase_align(img, img[max_rep, max_y, max_x, :], noise_threshold=noise_threshold)

    yield emit(aligned,
               dimension_labels=[mrd.ArrayDimension.REPETITION,
                                 mrd.ArrayDimension.Y,
                                 mrd.ArrayDimension.X,
                                 mrd.ArrayDimension.CONTRAST],
               description="epsi_image_aligned",
               noise_threshold=noise_threshold,
               phase_search_range=PHASE_SEARCH_RANGE,
               reference_repetition=max_rep,
               reference_y=max_y,
               reference_x=max_x)



    # # the phantom prescan, reconstructed on its own and on its own spectral axis: its matrix
    # # need not match the series', and on cirrhrat_43_1 it does not
    # phantom_scaling, phantom_map = 1.0, None
    # if phantom_acqs:
    #     phantom_scaling, phantom_map = reconstruct_phantom(
    #         header, phantom_acqs, line_broadening, width_scale)
    # print(f"Phantom scaling {phantom_scaling:.6g} "
    #       f"from {len(phantom_acqs)} prescan acquisition(s)", file=sys.stderr)
    # if phantom_map is not None:
    #     yield emit(phantom_map,
    #                dimension_labels=[mrd.ArrayDimension.Y, mrd.ArrayDimension.X],
    #                description="phantom_image",
    #                phantom_scaling=phantom_scaling)



    yield from fit_and_emit_peaks(aligned, global_spect, xscale, spec,
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
                    global_df: float = 0.5,
                    width_scale: Tuple[float, float] = (0.1, 1.9)) -> None:
    """
    Reconstruct one converted EPSI file

    The test runs before the writer is opened, so a non-EPSI file fails with this message
    rather than unwinding the writer mid protocol behind a ProtocolError about unwritten data.
    That is what lets the --folder loop drop the empty output and carry on.
    Raises:
        - ValueError if the file holds no acquisitions, or is not an EPSI acquisition
    """
    with mrd.BinaryMrdReader(input) as reader:
        header = reader.read_header()
        with mrd.BinaryMrdWriter(output) as writer:
            # add parameters used in recon to header
            header = append_recon_header(header,
                                         line_broadening=line_broadening,
                                         spec=spec,
                                         fit_df=fit_df,
                                         fit_dw=fit_dw,
                                         fit_dph=fit_dph,
                                         global_df=global_df,
                                         width_scale=width_scale)
            writer.write_header(header)
            # accumulate kspace across views, fft for each repetition
            user_long = {}
            for item in header.user_parameters.user_parameter_long:
                user_long[item.name] = item.value
            for item in reader.read_data():
                if isinstance(item, mrd.StreamItem.Acquisition):
                    # write back acq to recon file
                    acq = item.value
                    # sample period * total_points_per_switch including ramps and rephasing
                    spectral_bw_hz = 1.0e+9 / (acq.head.sample_time_ns * acq.samples() // user_long['nswitches'])
                    # the converter writes the 13C frequency here despite the field name, since that is the
                    # frequency these spectra were actually acquired at
                    center_freq_hz = header.experimental_conditions.h1resonance_frequency_hz
                    if not center_freq_hz:
                        raise ValueError("the header records no resonance frequency, so there is no ppm axis")
                    bw_ppm = spectral_bw_hz / center_freq_hz * 1.0e+6
                    nfreq = user_long['nswitches'] * FIDPAD
                    xscale = np.arange(nfreq) / nfreq * bw_ppm
                    if acq.head.flags & mrd.AcquisitionFlags.FIRST_IN_REPETITION:
                        # reference_acq to pass the acq header
                        enc = header.encoding[acq.head.encoding_space_ref].encoding_limits
                        img = np.zeros((enc.repetition.maximum + 1,         # repetition
                                        enc.phase.maximum + 1,              # ky
                                        user_long['npoints_per_switch'],    # kx
                                        user_long['nswitches']),            # t
                                        dtype='complex')
                    # accumulate each view into kspace
                    if acq.head.flags & mrd.AcquisitionFlags.FIRST_IN_PHASE:
                        kspace = np.zeros((enc.phase.maximum + 1,              # ky
                                           user_long['npoints_per_switch'],    # kx
                                           user_long['nswitches']),            # t
                                           dtype='complex')
                    # apply line broadening-> shape(points, switches)
                    kspace[acq.head.idx.kspace_encode_step_1] = apply_line_broadening(
                                                                    acq,
                                                                    line_broadening=line_broadening,
                                                                    nswitches=user_long['nswitches'],
                                                                    npoints_per_switch=user_long['npoints_per_switch'])
                    # once kspace is filled up with each phase, fft and store to img as single repetition
                    if acq.head.flags & mrd.AcquisitionFlags.LAST_IN_PHASE:
                        img[acq.head.idx.repetition] = fft_kspace_to_image(kspace)
                    if acq.head.flags & mrd.AcquisitionFlags.LAST_IN_REPETITION:
                        # initialize kspace for rawdata and prescan
                        writer.write_data(stream_epsi_image(header,
                                                            img,
                                                            reference_acq=acq,
                                                            line_broadening=line_broadening,
                                                            spec=spec,
                                                            fit_df=fit_df,
                                                            fit_dw=fit_dw,
                                                            fit_dph=fit_dph,
                                                            global_df=global_df,
                                                            width_scale=width_scale,
                                                            xscale=xscale))

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

    # the peak arguments have to come out before argparse sees them, since a peak's value may be
    # negative and argparse would read that as another option
    reserved = {option for action in parser._actions for option in action.option_strings}
    spec, remaining = split_peak_args(sys.argv[1:], reserved)
    args = parser.parse_args(remaining)

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
