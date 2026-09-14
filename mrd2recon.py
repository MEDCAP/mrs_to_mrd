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
                        width_scale: Tuple[float, float] = (0.1, 1.9),
                        drift_slope: float = 0.0) -> mrd.Header:
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
    # so a recon file says which window it was reconstructed at, not just which peaks were fitted
    header.user_parameters.user_parameter_double.append(
        mrd.UserParameterDoubleType(name="readout_drift_slope", value=float(drift_slope)))
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
                          totalppswitch: int,
                          npoints_per_switch: int,
                          drift_slope: float = 0.0) -> np.ndarray:
    """
    Read the flat top out of each switch of one readout and apodize it.

    The readout is laid out as (discard_pre, npoints_per_switch, discard_post) repeated once per
    switch, so cirrhrat_43_1's 1792 samples are 64 switches of 28 points with 12 kept. Only
    discard_pre and npoints_per_switch are needed to find that flat top; what follows it is
    whatever is left over, and the layout the caller passes has already been checked against
    every acquisition of the group.

    drift_slope lets the window follow the echo along the switch train, in samples per switch.
    The switch period is not exactly nsamples/nswitches - cirrhrat_43_1 runs about 28.19 against
    the 28 it records - so a window pinned at one offset drifts off the echo by the end of the
    train. Switch i is read from `discard_pre + slope * i` instead of from discard_pre.

    Every switch still yields exactly npoints_per_switch samples, whatever the slope, so the
    matrix this feeds keeps the same shape for every view: what varies is which samples are kept,
    never how many. That is what makes this the right place for the correction and a roll of the
    samples the wrong one - a roll is cyclic within the switch, so once the shift exceeds the
    discarded ramp it brings the rephasing ramp into the window. Indexing the flat sample axis
    rather than the reshaped grid reads across the switch boundary instead, which is where the
    drifting plateau's samples actually are.

    The window is clamped to the readout at both ends, which costs at most the last switch or
    two: at 0.1937 on cirrhrat_43_1 the final window wants sample 1792.2 of 1792.

    The apodization decays along the switch axis, because one spectral point is acquired per
    switch and so tk advances by a whole switch: that axis is the FID time axis. Which is why
    the exponent is applied after the transpose, where the switch axis is last. Moving the
    transpose would silently apodize along the readout instead.
    Returns a (npoints_per_switch, nswitches) complex array
    Raises:
        - ValueError if the readout is too short to hold the switches it should
    """
    samples = acq.data[0, :]
    if samples.shape[0] < nswitches * totalppswitch:
        raise ValueError(f"this readout holds {samples.shape[0]} samples, too few for "
                         f"{nswitches} switches of {totalppswitch} points")
    discard_pre = acq.head.discard_pre or 0
    # where each switch's flat top starts on the flat sample axis, so a window that runs past the
    # nominal switch boundary reads the next switch's samples rather than wrapping to its own start
    starts = np.rint(np.arange(nswitches) * (totalppswitch + drift_slope)
                     + discard_pre).astype(int)
    starts = np.clip(starts, 0, samples.shape[0] - npoints_per_switch)
    flat_top = samples[starts[:, None] + np.arange(npoints_per_switch)[None, :]]
    tk = np.arange(nswitches) * acq.head.sample_time_ns * totalppswitch / 1.0e+9
    return flat_top.T * np.exp(-tk * line_broadening)


def header_nswitches(header: mrd.Header) -> int:
    """
    The switch count, which is shared by every matrix in a file and recorded once on the header.

    Separate from layout_from_acq because the EPSI test in reconstruct_mrs has to run on the
    header alone, before the writer is opened and so before any acquisition has been read.
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


def acq_geometry(acq: mrd.Acquisition) -> Tuple:
    """The readout geometry every acquisition of one group has to agree on."""
    return (acq.samples(), acq.head.discard_pre or 0, acq.head.discard_post or 0,
            acq.head.encoding_space_ref)


def layout_from_acq(header: mrd.Header, acq: mrd.Acquisition) -> Layout:
    """
    The acquisition matrix a group was acquired at, read from one of its acquisitions.

    Three of the five numbers come from the readout, and are right for any matrix: nswitches is
    shared by every matrix in a file and is on the header as a user parameter, the whole switch
    is nsamples / nswitches, and the kept flat top is that less the discard counts. On
    cirrhrat_43_1 that gives 12 points a switch for the series and 18 for the averaged prescan
    beside it, which is correct for both - where the header's single npoints_per_switch, being
    one number for the whole file, can only be right for one of them.

    The view and repetition counts are the two a readout cannot supply, and the views are what
    differs between the matrices: 8 for that series, 12 for its prescan. The converter records a
    matrix per encoding and points each acquisition at its own through encoding_space_ref, so
    read them from there. Files converted before that leave the ref null and carry a single
    encoding describing the series only, so fall back to that one.

    The views are counted from kspace_encoding_step_1 rather than from phase, because
    kspace_encode_step_1 is the index the acquisitions are actually filed under. MRStomrd2
    writes the same count to both limits but never sets idx.phase at all, so phase is the field
    that would go stale first; it stays as the fallback, for a converter that set only it.

    One acquisition is enough because the group is homogeneous, which here is a checked claim
    rather than an assumption: RepetitionImages compares every later acquisition against this
    one, so a single mutated discard_post raises instead of quietly rewriting kept for the
    whole group.
    Raises:
        - ValueError if the header records no switch count, if the readout keeps no points, or
          if the encoding declares no view or repetition count to size the matrix from
    """
    nswitches = header_nswitches(header)
    if not nswitches:
        raise ValueError("this header records no switch count, so it describes no EPSI matrix")

    nsamples, discard_pre, discard_post, ref = acq_geometry(acq)
    totalppswitch = nsamples // nswitches
    kept = totalppswitch - discard_pre - discard_post
    if kept <= 0:
        raise ValueError(f"discard_pre={discard_pre} and discard_post={discard_post} leave "
                         f"{kept} points of a {totalppswitch} point switch")

    limits = None
    if ref is not None and ref < len(header.encoding):
        limits = header.encoding[ref].encoding_limits
    elif header.encoding:
        limits = header.encoding[0].encoding_limits
    views = None if limits is None else (limits.kspace_encoding_step_1 or limits.phase)
    reps = None if limits is None else limits.repetition
    if views is None or reps is None:
        raise ValueError(f"encoding {ref} declares no view or repetition count, so there is no "
                         f"matrix to reconstruct into")
    return Layout(nswitches=nswitches, totalppswitch=totalppswitch, kept=kept,
                  nviews=views.maximum + 1, nreps=reps.maximum + 1)


def fft_kspace_to_image(kspace: np.ndarray) -> np.ndarray:
    """FFT one repetition's k-space cube over all three axes into (views, readout, frequency)."""
    axes = tuple(range(kspace.ndim))
    return np.fft.fftshift(np.fft.fftn(kspace, axes=axes), axes=axes)


class RepetitionImages:
    """
    Accumulate one group's acquisitions into per-repetition images, a repetition at a time.

    A repetition is the smallest unit that can be transformed. One acquisition is one view, so
    the view axis is complete only once every acquisition of that repetition is in hand, and
    transforming all three axes of the cube is what turns it into (views, readout, frequency).
    Acquisitions arrive one at a time off the stream, so a cube is held open until the
    repetition index changes and is transformed then. Only the image outlives it, which keeps a
    single repetition of k-space in memory rather than the whole series, while still retaining
    every image: the alignment reference is the brightest voxel across the whole series and the
    global spectrum is the sum over it, so everything after the transform needs the series in
    hand.

    The repetition axis of `images` is the recorded repetition index rather than a running count
    of the repetitions delivered, so it lines up with acq.head.idx.repetition and with the
    repetition axis of the metabolite maps, which makes it a time axis a kinetic model can use.
    A repetition the stream never delivers is left as zeros and is skipped by both fits.
    """

    def __init__(self, header: mrd.Header, first_acq: mrd.Acquisition, line_broadening: float,
                 drift_slope: float = 0.0):
        self.layout = layout_from_acq(header, first_acq)
        self.line_broadening = line_broadening
        self.drift_slope = drift_slope
        self._geometry = acq_geometry(first_acq)
        self.images = np.zeros((self.layout.nreps, self.layout.nviews, self.layout.kept,
                                self.layout.nswitches * FIDPAD), dtype='complex')
        # the first acquisition of each repetition, whose header describes that repetition
        self.reference_acq: dict = {}
        self.reps_present: List[int] = []
        self._rep: Optional[int] = None
        self._cube: Optional[np.ndarray] = None
        self._views_seen: set = set()

    def add(self, acq: mrd.Acquisition) -> None:
        """
        Take one acquisition, transforming the repetition before it once this one starts.

        The geometry check is per acquisition rather than over an assembled group, because a
        stream is only ever seen one acquisition at a time. It is the same guarantee the batch
        check gave: a minority layout is a corrupted header, not a variant, and a group that
        disagrees with itself means something upstream is wrong that a quietly reasonable answer
        would hide. Taking the layout from the first acquisition alone once let a single mutated
        discard_post rewrite kept from 12 to 20 for all 216 acquisitions of a series, silently.
        Raises:
            - ValueError if this acquisition's readout geometry differs from the group's, if it
              carries a view or repetition index the encoding does not declare, or if its
              repetition has already been transformed
        """
        geometry = acq_geometry(acq)
        if geometry != self._geometry:
            raise ValueError(f"this group mixes readout geometries: (samples, discard_pre, "
                             f"discard_post, encoding_ref)={geometry} against {self._geometry} "
                             f"on the first acquisition of the group")

        rep = acq.head.idx.repetition or 0
        if rep >= self.layout.nreps:
            raise ValueError(f"the encoding declares {self.layout.nreps} repetitions but "
                             f"repetition index {rep} is present")
        if rep != self._rep:
            self._transform_open_cube()
            # the cube for a repetition is opened once and transformed when the next one starts,
            # so a repetition that comes back later would be rebuilt from its tail alone
            if rep in self.reps_present:
                raise ValueError(f"repetition {rep} arrives in more than one run of the stream, "
                                 f"so its k-space cannot be assembled in a single pass")
            self._rep = rep
            self._cube = np.zeros((self.layout.nviews, self.layout.kept,
                                   self.layout.nswitches * FIDPAD), dtype='complex')
            self._views_seen = set()
            self.reference_acq[rep] = acq

        view = acq.head.idx.kspace_encode_step_1 or 0
        if view >= self.layout.nviews:
            raise ValueError(f"the encoding declares {self.layout.nviews} views but view index "
                             f"{view} is present")
        self._cube[view, :, :self.layout.nswitches] = apply_line_broadening(
            acq, self.line_broadening,
            nswitches=self.layout.nswitches,
            totalppswitch=self.layout.totalppswitch,
            npoints_per_switch=self.layout.kept,
            drift_slope=self.drift_slope)
        self._views_seen.add(view)

    def close(self) -> None:
        """Transform the repetition still open, once the stream has run out."""
        self._transform_open_cube()
        absent = sorted(set(range(self.layout.nreps)) - set(self.reps_present))
        if absent:
            print(f"warning: the stream carries none of repetitions {absent} that the encoding "
                  f"declares, which stay zero through the reconstruction", file=sys.stderr)

    def _transform_open_cube(self) -> None:
        """Transform the open cube into its image and drop it. Idempotent, so close() is safe."""
        if self._cube is None:
            return
        missing = sorted(set(range(self.layout.nviews)) - self._views_seen)
        if missing:
            print(f"warning: repetition {self._rep} is missing views {missing}, which keep zero "
                  f"rows through the transform", file=sys.stderr)
        self.images[self._rep] = fft_kspace_to_image(self._cube)
        self.reps_present.append(self._rep)
        self._cube = None


def spectral_axis(header: mrd.Header,
                  acq: mrd.Acquisition,
                  nswitches: int,
                  totalppswitch: int) -> Tuple[np.ndarray, float]:
    """
    The spectral axis of an EPSI readout, in ppm.

    sample_time_ns is the dwell time of a single point, and one spectral point is acquired per
    readout switch, so the spectral sampling interval is a whole switch. The axis deliberately
    runs from 0 rather than being centred on zero: the peak centers are placed modulo the
    spectral width and the Lorentzian model carries explicit +-BW wraparound terms.
    Returns:
        - (xscale in ppm, spectral bandwidth in Hz)
    """
    spectral_bw_hz = 1.0e+9 / (acq.head.sample_time_ns * totalppswitch)
    # the converter writes the 13C frequency here despite the field name, since that is the
    # frequency these spectra were actually acquired at
    center_freq_hz = header.experimental_conditions.h1resonance_frequency_hz
    if not center_freq_hz:
        raise ValueError("the header records no resonance frequency, so there is no ppm axis")
    bw_ppm = spectral_bw_hz / center_freq_hz * 1.0e+6
    nfreq = nswitches * FIDPAD
    xscale = np.arange(nfreq) / nfreq * bw_ppm
    return xscale, spectral_bw_hz


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
                # a repetition the stream never delivered is a zero slab, and stays one
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
                       width_scale: Tuple[float, float] = (0.1, 1.9),
                       phantom_scaling: float = 1.0) -> Iterable[mrd.StreamItem]:
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
                    fit_dph_rad=fit_dph,
                    # derived from the phantom prescan and recorded, not applied: it is what
                    # makes two experiments comparable, and that is the consumer's decision
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


def reconstruct_phantom(header: mrd.Header,
                        prescan: RepetitionImages,
                        width_scale: Tuple[float, float] = (0.1, 1.9)
                        ) -> Tuple[float, Optional[np.ndarray]]:
    """
    Fit the urea phantom prescan and derive the scaling factor it exists to provide.

    The phantom holds one known metabolite, so its spectrum is fitted with a single peak rather
    than the metabolite pattern, and the quantity of interest is that peak's area, amplitude
    times width: that is what makes two experiments comparable. No noise gate applies here -
    every voxel of a phantom is signal, which is why the legacy sets a phantom's noise to zero
    outright. Several phantom sets average into one number.

    The prescan arrives already transformed, in its own RepetitionImages, because its matrix
    need not match the series' and on cirrhrat_43_1 it does not. That accumulator carries the
    prescan's own layout, which is what this reads its spectral axis off.

    width_scale is the same (lo, hi) bound the metabolite fit uses, as multiples of the width
    guess. It matters more here than there. The scaling is the peak's area, amplitude times
    width, so a width held off its true value scales every experiment this number normalises.
    Returns:
        - (the scaling, 1.0 when there is no phantom to derive one from,
           the per-voxel area map of the last phantom set, or None)
    """
    if not prescan.reps_present:
        return 1.0, None
    xscale, _ = spectral_axis(header, prescan.reference_acq[prescan.reps_present[0]],
                              prescan.layout.nswitches, prescan.layout.totalppswitch)

    scalings: List[float] = []
    voxel_map = None

    for rep in prescan.reps_present:
        img = prescan.images[rep]
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


def reconstruct_epsi(header: mrd.Header,
                     series: RepetitionImages,
                     prescan: Optional[RepetitionImages],
                     *,
                     line_broadening: float,
                     spec: PeakSpec,
                     fit_df: float = 0.0,
                     fit_dw: float = 0.0,
                     fit_dph: float = 0.0,
                     global_df: float = 0.5,
                     width_scale: Tuple[float, float] = (0.1, 1.9)
                     ) -> Iterable[mrd.StreamItem]:
    """
    Reconstruct an EPSI acquisition into spectra, a global peak fit and metabolite maps.

    Everything downstream of the transform. The series arrives already in the image domain, as
    the (repetitions, views, readout, frequency) array its RepetitionImages accumulated one
    repetition at a time off the stream, so nothing here touches k-space or an acquisition's
    data. What it still needs from the acquisitions is their headers, which the accumulator
    keeps one per repetition, and which carry the timestamps that make the repetition axis a
    time axis.

    The averaged phantom prescan arrives as its own RepetitionImages, on its own matrix, and
    contributes only its scaling factor and its own image.

    Both sides of the alignment are emitted. The pre-alignment image goes out per repetition,
    carrying that repetition's own acquisition header, and is the record of what the alignment
    was given; the aligned series goes out once as a 4-D array, and is what the fits actually
    see and what a consumer needs to re-fit without realigning.

    Yields NdArrays only. The raw acquisitions are passed through by the caller, in a stream
    call of their own.
    """
    if not series.reps_present:
        return

    layout = series.layout
    recon_array = series.images
    reps_present = series.reps_present
    xscale, spectral_bw_hz = spectral_axis(header, series.reference_acq[reps_present[0]],
                                           layout.nswitches, layout.totalppswitch)

    for rep in reps_present:
        reference_acq = series.reference_acq[rep]
        yield emit(recon_array[rep],
                   head=mrd.NdArrayHeader(
                       dimension_labels=[mrd.ArrayDimension.Y,
                                         mrd.ArrayDimension.X,
                                         mrd.ArrayDimension.CONTRAST],
                       # the schema has no arm for "a reconstructed spectral image", so it
                       # goes in the generic bucket and the description identifies it
                       array_type=mrd.ArrayType.USER_MAP,
                       image_type=mrd.ArrayImageType.COMPLEX,
                       measurement_uid=reference_acq.head.measurement_uid,
                       average=reference_acq.head.idx.average,
                       repetition=rep,
                       acquisition_time_stamp_ns=reference_acq.head.acquisition_time_stamp_ns),
                   description="epsi_image",
                   xscale_ppm=xscale,
                   **{"receiver bandwidth(Hz)": spectral_bw_hz,
                      "line broadening(Hz)": line_broadening})

    # the phantom prescan, fitted on its own spectral axis: its matrix need not match the
    # series', and on cirrhrat_43_1 it does not
    phantom_scaling, phantom_map = 1.0, None
    nprescan = 0 if prescan is None else len(prescan.reps_present)
    if prescan is not None:
        phantom_scaling, phantom_map = reconstruct_phantom(header, prescan, width_scale)
    print(f"Phantom scaling {phantom_scaling:.6g} from {nprescan} prescan repetition(s)",
          file=sys.stderr)
    if phantom_map is not None:
        yield emit(phantom_map,
                   dimension_labels=[mrd.ArrayDimension.Y, mrd.ArrayDimension.X],
                   description="phantom_image",
                   phantom_scaling=phantom_scaling)

    # the last repetition of a hyperpolarized series has decayed away, so it measures noise.
    # It is the last one delivered rather than recon_array[-1], which is zeros when the stream
    # stopped short of the declared count. A repetition missing views dilutes this downward,
    # since its zero rows count toward the mean
    noise_repetition = max(reps_present)
    noise = float(np.mean(np.abs(recon_array[noise_repetition])))
    noise_threshold = noise * NOISE_THRESHOLD_MULTIPLIER
    if noise == 0.0:
        print(f"warning: repetition {noise_repetition} is empty, so the noise floor is 0 and "
              f"every voxel above zero will be fitted", file=sys.stderr)
    yield emit(np.array([noise]),
               dimension_labels=[mrd.ArrayDimension.SAMPLES],
               array_type=mrd.ArrayType.NOISE,
               description="noise",
               noise_threshold_multiplier=NOISE_THRESHOLD_MULTIPLIER,
               estimated_from_repetition=noise_repetition)

    # the brightest voxel of the series is the reference every other voxel is aligned to. One
    # argmax over the series picks the same voxel the per-repetition scan it replaces did: both
    # take the first maximum, and C order over this array is repetition order
    peak_per_voxel = np.abs(recon_array).max(axis=-1)
    max_rep, max_y, max_x = (int(i) for i in
                             np.unravel_index(int(np.argmax(peak_per_voxel)),
                                              peak_per_voxel.shape))
    max_spect = np.copy(recon_array[max_rep, max_y, max_x, :])

    yield emit(max_spect,
               dimension_labels=[mrd.ArrayDimension.CONTRAST],
               description="max spectral",
               xscale_ppm=xscale,
               max_repetition=max_rep,
               max_y=max_y,
               max_x=max_x)

    print("Aligning voxel spectra", file=sys.stderr)
    aligned, global_spect = phase_align(recon_array, max_spect, noise_threshold=noise_threshold)

    yield emit(aligned,
               dimension_labels=[mrd.ArrayDimension.REPETITION,
                                 mrd.ArrayDimension.Y,
                                 mrd.ArrayDimension.X,
                                 mrd.ArrayDimension.CONTRAST],
               description="epsi_image_aligned",
               xscale_ppm=xscale,
               noise_threshold=noise_threshold,
               phase_search_range=PHASE_SEARCH_RANGE,
               reference_repetition=max_rep,
               reference_y=max_y,
               reference_x=max_x)

    yield from fit_and_emit_peaks(aligned, global_spect, xscale, spec,
                                  noise_threshold=noise_threshold,
                                  fit_df=fit_df, fit_dw=fit_dw, fit_dph=fit_dph,
                                  global_df=global_df, width_scale=width_scale,
                                  phantom_scaling=phantom_scaling)


def reconstruct_mrs(input: BinaryIO,
                    output: BinaryIO,
                    *,
                    line_broadening: float,
                    spec: PeakSpec,
                    fit_df: float = 0.0,
                    fit_dw: float = 0.0,
                    fit_dph: float = 0.0,
                    global_df: float = 0.5,
                    width_scale: Tuple[float, float] = (0.1, 1.9),
                    drift_slope: float = 0.0) -> None:
    """
    Reconstruct one converted EPSI file

    The stream is read once. Each acquisition is written straight back out and, on the way past,
    added to its group's accumulator, which transforms a repetition as soon as the next one
    starts. So the raw acquisitions reach the output unchanged, only one repetition of k-space
    is ever held, and by the time the passthrough is exhausted every image the fits need is in
    hand. Acquisitions and NdArrays go out in two stream calls, acquisitions first, which is the
    order a consumer reads them in.

    Acquisitions flagged IS_NOISE_MEASUREMENT are the averaged phantom prescan, not repetitions
    of the series, so they accumulate into a second group with its own matrix. Which group an
    acquisition belongs to is the only thing this decides; everything else about a group is the
    accumulator's.

    nswitches > 1 is the EPSI test, and it is the converter's own: MRStomrd2 sets discard_pre
    and discard_post only under that condition, so it is exactly the set of files whose readouts
    are switch structured. The header carries nswitches for every conversion, so its presence
    alone says nothing - a spectral file records nswitches=1. The test runs before the writer is
    opened, so a non-EPSI file fails with this message rather than unwinding the writer mid
    protocol behind a ProtocolError about unwritten data. A file that turns out to hold no
    acquisitions cannot be caught that early, so it lets the writer close on a header alone and
    raises after. Either way the --folder loop drops the output and carries on.
    Raises:
        - ValueError if the file is not an EPSI acquisition, or holds no acquisitions
    """
    with mrd.BinaryMrdReader(input) as reader:
        header = reader.read_header()
        nswitches = header_nswitches(header)
        if nswitches <= 1:
            raise ValueError(f"this header records nswitches={nswitches}; mrd2recon "
                             f"reconstructs EPSI only")

        # keyed by whether the acquisition is prescan, so a group's accumulator is built from
        # the first acquisition that belongs to it and sizes itself off that acquisition's own
        # encoding
        groups: dict = {}

        def passthrough() -> Iterator[mrd.StreamItem]:
            """Write every acquisition back out unchanged, accumulating it on the way past."""
            for item in reader.read_data():
                if not isinstance(item, mrd.StreamItem.Acquisition):
                    continue
                acq = item.value
                is_prescan = bool(acq.head.flags & mrd.AcquisitionFlags.IS_NOISE_MEASUREMENT)
                if is_prescan not in groups:
                    groups[is_prescan] = RepetitionImages(header, acq, line_broadening,
                                                          drift_slope)
                groups[is_prescan].add(acq)
                yield item

        with mrd.BinaryMrdWriter(output) as writer:
            header = append_recon_header(header,
                                         line_broadening=line_broadening,
                                         spec=spec,
                                         fit_df=fit_df,
                                         fit_dw=fit_dw,
                                         fit_dph=fit_dph,
                                         global_df=global_df,
                                         width_scale=width_scale,
                                         drift_slope=drift_slope)
            writer.write_header(header)
            # the raw acquisitions pass through unchanged, and are accumulated as they go
            writer.write_data(passthrough())
            for group in groups.values():
                group.close()

            series = groups.get(False)
            if series is not None:
                writer.write_data(reconstruct_epsi(header,
                                                   series,
                                                   groups.get(True),
                                                   line_broadening=line_broadening,
                                                   spec=spec,
                                                   fit_df=fit_df,
                                                   fit_dw=fit_dw,
                                                   fit_dph=fit_dph,
                                                   global_df=global_df,
                                                   width_scale=width_scale))

    if series is None:
        raise ValueError("this file holds no acquisitions")


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
    parser.add_argument("-gdf", "--global-df", type=float, default=0.5, required=False,
                        help="How far the global fit may move a peak center from its rigid-pattern placement, in ppm. Default 0.5; 0 pins the placement")
    parser.add_argument("--width-bounds", type=float, nargs=2, default=(0.1, 1.9), metavar=("LO", "HI"), required=False,
                        help="Linewidth bounds for both fits, as multiples of the width guess. Default 0.1 1.9")
    parser.add_argument("--drift-slope", type=float, default=0.0, required=False,
                        help="How far the readout window follows the echo along the switch train, "
                             "in samples per switch. Default 0, i.e. the same window in every "
                             "switch. The switch period is not exactly nsamples/nswitches, so a "
                             "pinned window drifts off the echo by the end of a long train")

    # the peak arguments have to come out before argparse sees them, since a peak's value may be
    # negative and argparse would read that as another option
    reserved = {option for action in parser._actions for option in action.option_strings}
    spec, remaining = split_peak_args(sys.argv[1:], reserved)
    args = parser.parse_args(remaining)

    recon_kwargs = dict(line_broadening=args.line_broadening,
                        spec=spec,
                        fit_df=args.fit_df,
                        fit_dw=args.fit_dw,
                        fit_dph=args.fit_dph,
                        global_df=args.global_df,
                        width_scale=tuple(args.width_bounds),
                        drift_slope=args.drift_slope)

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
                    try:
                        with open(input_filepath, "rb") as input, open(recon_filepath, "wb") as output:
                            reconstruct_mrs(input, output, **recon_kwargs)
                    except ValueError as err:
                        # a folder holds whatever was converted into it, and this reconstruction is
                        # EPSI only. One spectral scan among the inputs should not abandon the rest,
                        # so drop the empty output it left behind and carry on
                        recon_filepath.unlink(missing_ok=True)
                        print(f"  skipping {input_filepath.name}: {err}", file=sys.stderr)
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
