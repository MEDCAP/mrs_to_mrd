"""
Take the EPSI echo drift out of a converted MRD v2 stream.

The middle stage of the pipeline:  tar -> MRStomrd2.py -> mrd2shift.py -> mrd2recon.py

THE DEFECT.  The scanner records no_samples and no_switches separately, so everything downstream
divides: 1792/64 = 28 samples per switch on the cirrhrat data.  That division is wrong - the real
gradient period is 28.1937 samples.  Block i therefore starts 0.1937*i samples out of step with the
gradient cycle it describes, and the echo appears to walk +12.2 samples across a 64 switch train.
The reconstruction reads the same window out of every switch, so early switches are read correctly
and late ones collect gradient ramp instead of signal.

Two axes, and which is which matters: position within a switch is kx, the readout direction; the
switch index is time, and becomes the spectral axis after the transform.

THE MEASUREMENT.  fit_peak_lines fits the per-switch peaks as two parallel lines sharing one slope.
Two lines because a switch crosses k-space centre twice - the readout echo on the plateau and the
rephasing echo after it - so the brightest position in a switch is whichever of them won there, and
what looks like scatter is two orderly families.  Sharing the slope means every switch constrains
it rather than only the ones that happened to peak on the readout echo.

THE BASE.  Before any of that, the whole readout is moved onto the position the sequence asks for by
pushing zeros in at its front - measure_pad works out how many from where switch 0's readout echo
sits against the `discard_pre + npoints_per_switch/2` the sequence puts it at, and prepend_zeros
applies it.  Switch 0 is barely displaced by the drift, so once it is on target the first few
switches can be trusted and everything left to correct is the accumulation along the train.  The
number this measures is the `ramp + 3` an older converter hard-coded per tramp: +7 on the tramp 112
cirrhrat data, +5 on the tramp 100 kidney data.  --pad overrides it, and --pad 0 turns it off.

THE CORRECTION.  On top of that base, each switch i is displaced by delta_i = -slope * i, anchored on
switch 0 so the train is pulled back onto the position the base put switch 0 at.  Where that
displacement is applied is the whole question, and --method selects it:

    regrid      exact, resampling the whole readout onto the period it was acquired at (default)
    contiguous  exact, over the whole switch, along the contiguous readout rather than wrapping
    shift       exact, over the whole switch, cyclically
    roll        the fractional part rounded away, over the contiguous readout      (control)
    alloc       the fractional part linearly interpolated, over the contiguous readout (control)

A `phase` method, a phase ramp applied inside the reconstruction's window alone, was dropped: it
leaves the readout where it was by construction, so it cannot straighten the echo, and measured on
cirrhrat_43_1 it made the global fit worse than no correction at all (residual 1.648 against 1.284).

    python mrd2shift.py -i raw.mrd2 -o straight.mrd2
    cat raw.mrd2 | python mrd2shift.py -i - -o - | ...
    python mrd2shift.py -i raw.mrd2 -p methods.png          # every method drawn, nothing written

Both default to $INPUT_PIPE/$OUTPUT_PIPE, so a Tyger codespec needs no arguments.  The whole stream
is held in memory: the drift is measured from every acquisition before the first can be written, and
a buffer FIFO can only be read once.  Every diagnostic goes to stderr - anything on stdout lands in
the middle of the stream and the next stage dies on its magic bytes.
"""

from __future__ import annotations

import argparse
import contextlib
import os
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

import mrd

# How wide a drift the search covers, in samples per switch, and how finely.  Half a sample sweeps
# the echo one and a half switch widths over a 64 switch train; the kidney scans run -0.03 to +0.27
# and cirrhrat_43_1 fits +0.1937.
DRIFT_SLOPE_LIMIT = 0.5
DRIFT_FINE_STEP = 0.002
# How far off its own fitted line a family may sit before the peaks are called noise.  Real data
# reads 0.4-0.9 positions; a readout with no echo reads 2.1-2.6 with the switches split about evenly.
DRIFT_MAX_RMS = 1.5
# What was taken out, recorded so a second pass refuses rather than correcting twice.
SLOPE_PARAMETER = "echo_drift_slope_applied"
PAD_PARAMETER = "echo_drift_pad_applied"
METHOD_PARAMETER = "echo_drift_method"
METHODS = ("shift", "regrid", "contiguous", "roll", "alloc")
DEFAULT_METHOD = "regrid"


# ---------- the stream ------------------------------------------------------------------------


@contextlib.contextmanager
def open_binary(path: str, mode: str):
    """A binary stream for a path, or stdin/stdout when the path is -, which is left unclosed."""
    if path == "-":
        yield sys.stdin.buffer if mode == "rb" else sys.stdout.buffer
        return
    with open(path, mode) as stream:
        yield stream


def read_stream(path: str) -> Tuple[mrd.Header, List[mrd.StreamItem]]:
    """Read a whole stream into memory, keeping non-acquisition items so a rewrite passes them on."""
    # opening a FIFO for read blocks until the buffer sidecar opens the write end
    with open_binary(path, "rb") as stream:
        with mrd.BinaryMrdReader(stream) as reader:
            header = reader.read_header()
            items = list(reader.read_data())
    return header, items


def write_stream(path: str, header: mrd.Header, items: Sequence[mrd.StreamItem]) -> None:
    """Write the header and every item back out, opened even when nothing changed."""
    with open_binary(path, "wb") as output:
        # closing the writer is what emits the end-of-stream sentinel, hence the with block
        with mrd.BinaryMrdWriter(output) as writer:
            writer.write_header(header)
            writer.write_data(items)


def header_long(header: mrd.Header, name: str) -> Optional[int]:
    """The named long user parameter, or None when the header does not carry it."""
    user = getattr(header, "user_parameters", None)
    for param in (getattr(user, "user_parameter_long", None) or []):
        if param.name == name:
            return int(param.value)
    return None


def switch_layout(header: mrd.Header, acq: mrd.Acquisition) -> Tuple[int, int, int, int]:
    """
    (nswitch, total, kept, start) - the switch train, and the window the reconstruction reads.

    MRStomrd2 records nswitches and npoints_per_switch on the header and discard_pre per
    acquisition, and mrd2recon slices `discard_pre : discard_pre + npoints_per_switch`, so that
    slice is what a correction aimed at the reconstruction has to target.

    `kept` comes from this acquisition's own discards where they account for the switch, and only
    then from the header.  npoints_per_switch is written once for a whole file, and a file holds the
    series and a prescan acquired on a different switch period: on cirrhrat_43_1 the header says 12,
    which is the series' kept and closes its 28 point switch, while the prescan's 34 point switch is
    `4 + 18 + 3*4` and closes at 18.  Reading 12 off the header for both puts the prescan's window
    and its expected echo position six samples from where they are.
    """
    nswitch = header_long(header, "nswitches") or 1
    total = acq.samples() // max(nswitch, 1)
    start = acq.head.discard_pre or 0
    own = total - start - (acq.head.discard_post or 0)
    if own > 0 and 4 * start + own == total:
        return nswitch, total, own, start
    return nswitch, total, (header_long(header, "npoints_per_switch") or 0) or total, start


def acquisition_cube(acqs: Sequence[mrd.Acquisition], nswitch: int,
                     total: int) -> Tuple[np.ndarray, int]:
    """
    One layout's acquisitions as a (switch, position, view, repetition) cube.

    Views and repetitions are taken from the indices the converter set rather than their counts, so
    a group missing one still lands in the right column.
    """
    views = sorted({int(a.head.idx.kspace_encode_step_1 or 0) for a in acqs})
    reps = sorted({int(a.head.idx.repetition or 0) for a in acqs})
    view_at = {v: i for i, v in enumerate(views)}
    rep_at = {r: i for i, r in enumerate(reps)}

    used = nswitch * total
    cube = np.zeros((nswitch, total, len(views), len(reps)), dtype=np.complex128)
    filled = np.zeros((len(views), len(reps)), dtype=bool)
    for acq in acqs:
        samples = np.asarray(acq.data)
        # MRS acquires on one channel; summed rather than assumed, so a multi coil stream still reads
        line = samples[0] if samples.shape[0] == 1 else samples.sum(axis=0)
        v = view_at[int(acq.head.idx.kspace_encode_step_1 or 0)]
        r = rep_at[int(acq.head.idx.repetition or 0)]
        cube[:, :, v, r] = line[:used].reshape(nswitch, total)
        filled[v, r] = True
    return cube, int((~filled).sum())


# Both, because the converter has used each in turn: MRStomrd2 marks an averaged prescan
# IS_NOISE_MEASUREMENT today and mrd2recon reads that, while several docstrings and fid_recon still
# name IS_NAVIGATION_DATA. Accepting either is what lets this read a stream from before or after
# that change, and mrd2_recon_from_main.PRESCAN_FLAGS does the same for the same reason.
PRESCAN_FLAGS = int(mrd.AcquisitionFlags.IS_NOISE_MEASUREMENT) | int(
    mrd.AcquisitionFlags.IS_NAVIGATION_DATA)


def is_prescan(acq: mrd.Acquisition) -> bool:
    """Whether this acquisition is the averaged prescan that calibrates the data, not data."""
    return bool(int(acq.head.flags) & PRESCAN_FLAGS)


# ---------- measuring the drift ---------------------------------------------------------------


def fit_peak_lines(signal: np.ndarray) -> dict:
    """
    Fit the two parallel lines the per-switch peaks lie on.

    A switch crosses k-space centre twice, so its brightest position is whichever echo won there and
    a plot of those positions is two straight lines, not one line with outliers.  One slope is fitted
    through every switch of both families with an intercept each, so all of them constrain the number
    the correction uses, and the separation is measured from the same fit.
    Args:
        - signal: (switch, position) magnitude
    Returns:
        - the peaks, the shared slope, each line's intercept and rms, and which family each switch
          fell in
    """
    nswitch, total = signal.shape
    switches = np.arange(nswitch)
    peaks = np.argmax(signal, axis=1)

    # the one slope and separation leaving every switch nearest to one of two parallel lines,
    # searched rather than derived because the assignment depends on the slope and vice versa
    # The separation is not a free parameter and must not be searched.  The two crossings of k-space
    # centre sit ppsw/2 + 2*ramp apart, and since period = ppsw + 4*ramp that is identically
    # period/2 - so ppsw cancels and the spacing carries no information beyond the period.  Worse,
    # period/2 is exactly where a two-line model cannot resolve: line + period/2 and line - period/2
    # are the same line relabelled, so least squares is *repelled* from the truth.  Measured on
    # cirrhrat_43_1 the loss reads 22.2 at separation 12, 81.3 at the true 14, and 22.2 again at 16.
    # Pinning it leaves the shared slope, which is the only thing here that was ever identifiable.
    separation = total / 2
    best = None
    for slope in np.arange(-DRIFT_SLOPE_LIMIT, DRIFT_SLOPE_LIMIT + DRIFT_FINE_STEP / 2,
                           DRIFT_FINE_STEP):
        residual = (peaks - peaks[0] - slope * switches + total / 2) % total - total / 2
        nearer = np.minimum(np.abs(residual),
                            np.abs((residual - separation + total / 2) % total - total / 2))
        cost = float((nearer ** 2).sum())
        if best is None or cost < best[0]:
            best = (cost, float(slope))
    _, slope = best

    residual = (peaks - peaks[0] - slope * switches + total / 2) % total - total / 2
    off_second = (residual - separation + total / 2) % total - total / 2
    first = np.abs(residual) <= np.abs(off_second)
    # unwrapped about whichever line each switch fell on, so a family crossing the switch boundary is
    # a straight line to fit rather than one that jumps by `total` at the crossing
    unwrapped = peaks[0] + slope * switches + np.where(first, residual, separation + off_second)

    counts = (int(first.sum()), int((~first).sum()))
    if min(counts) >= 2:
        design = np.column_stack((switches, first, ~first)).astype(float)
        solution = np.linalg.lstsq(design, unwrapped, rcond=None)[0]
        slope, first_c, second_c = float(solution[0]), float(solution[1]), float(solution[2])
    else:
        fitted = np.polyfit(switches, unwrapped, 1)
        slope, first_c = float(fitted[0]), float(fitted[1])
        second_c = first_c + separation

    def spread(family, intercept):
        if not family.any():
            return float('nan')
        return float(np.sqrt(((unwrapped[family] - (intercept + slope * switches[family])) ** 2).mean()))

    return dict(peaks=peaks, slope=slope, separation=float(second_c - first_c),
                first_intercept=first_c, second_intercept=second_c, first_family=first,
                on_first=counts[0], on_second=counts[1],
                first_rms=spread(first, first_c), second_rms=spread(~first, second_c))


def measure_drift(cube: np.ndarray, nswitch: int, total: int) -> dict:
    """
    How fast the echo walks along the switch train, fitted from the per-switch peaks.

    This replaces an FFT-coherent search that scored candidate slopes by how sharp they left the
    readout profile.  On cirrhrat_43_1 that search refuses the scan outright and, forced, returns
    +0.2700 where the truth is +0.1937: de-drifting at +0.1937 and refitting leaves -0.001 residual
    and lands the readout echo on the position the sequence predicts, which +0.2700 does not.
    Returns:
        - usable, reason, the slope, the drift over the train, the implied period, and the fit
    """
    lines = fit_peak_lines(np.abs(cube).sum(axis=(2, 3)))
    slope = float(lines['slope'])
    drift = slope * (nswitch - 1)
    spread = [r for r in (lines['first_rms'], lines['second_rms']) if np.isfinite(r)]
    worst = max(spread) if spread else float('inf')

    report = dict(usable=False, reason="", slope=slope, drift=drift, period=total + slope,
                  rms=worst, lines=lines, nswitch=nswitch, total=total)
    if worst >= DRIFT_MAX_RMS:
        report['reason'] = (f"the per switch peaks sit {worst:.2f} positions off the lines fitted "
                            f"through them, past the {DRIFT_MAX_RMS} that separates a readout with "
                            f"an echo from one where the argmax is landing on noise")
    elif abs(drift) < 1.0:
        report['reason'] = (f"the fitted drift, {drift:+.2f} samples over {nswitch} switches, is "
                            f"below the one sample it would take to move anything")
    else:
        report['usable'] = True
        report['reason'] = (f"{drift:+.2f} samples over {nswitch} switches, {slope:+.4f} per switch, "
                            f"through {lines['on_first']}/{lines['on_second']} switches on the two "
                            f"echoes, rms {worst:.2f}")
    return report


def switch_offsets(nswitch: int, slope: float) -> np.ndarray:
    """
    The displacement per switch, anchored on switch 0.

    delta_i = -slope * i, so switch 0 does not move and every later switch is pulled back onto it.
    Anchored there and not on the middle of the train because the zero-pad base is what decides the
    absolute position, and it decides it from switch 0: switch 0 carries the least drift and the most
    signal, so it is the one reading of the echo position worth trusting.  A mid-train anchor would
    hold the *mean* position instead, and the two would compose into a readout sitting half the
    drift - six samples of a twelve point window on cirrhrat - past where the pad aimed it.

    The cost is that the largest displacement now falls at the end of the train, slope * (M-1)
    rather than half that.  `shift` is the only method that still wraps inside a switch, so it is the
    only one for which that means more of the switch's own ramp pulled into the window at the far
    end; contiguous, regrid, roll and alloc all read from the neighbouring switch instead and do not
    pay it.  Which is a difference -p is there to show.
    """
    return -slope * np.arange(nswitch)


def geometry_closes(total: int, kept: int, start: int) -> bool:
    """
    Whether a layout's own numbers account for its switch, `ramp + kept + 3 ramps == total`.

    The epsi readout loop runs rampup -> kept -> rampdown -> rephasing, so a switch holds
    `4 * ramp + npoints_per_switch` samples, and `discard_pre` is the ramp.  When that does not add
    up, the layout is not described by the numbers attached to it and nothing derived from them -
    above all `start + kept // 2`, where the echo is supposed to be - can be trusted.

    Both of cirrhrat_43_1's layouts do close, once switch_layout takes `kept` from each
    acquisition's own discards rather than from the file-wide header: the series at `4 + 12 + 3*4`
    = 28, the prescan at `4 + 18 + 3*4` = 34.  This stays as the guard for the case where they do
    not - the kidney data, where `4*2.5 + 12` is 22 against a 20 point switch, is one.
    """
    return 4 * start + kept == total


def readout_anchor(lines: dict, signal: np.ndarray, total: int) -> float:
    """
    Where switch 0's *readout* echo sits, read off the fitted lines rather than its own argmax.

    A switch crosses k-space centre twice and the two crossings sit exactly total/2 apart, so the two
    candidate anchors are equally far from any target and proximity cannot choose between them.
    Brightness can: the readout echo is on the gradient plateau where the signal is, and the
    rephasing one is not, so the family whose peaks average brighter is the readout family.

    Read off the fitted line rather than off switch 0 itself because switch 0 peaks on whichever of
    its two echoes won there.  On cirrhrat_43_1 that is the rephasing one - its argmax reads 17 where
    its neighbours read 3, 4, 5 - and anchoring on that would aim the whole correction at the wrong
    echo, putting the readout echo at position 1 of 28.
    Args:
        - lines: what fit_peak_lines returned for `signal`
        - signal: (switch, position) magnitude, the array `lines` was fitted from
    Returns:
        - switch 0's readout echo position, wrapped into one switch
    """
    peaks = lines['peaks']
    first = lines['first_family']
    brightness = signal[np.arange(len(peaks)), peaks]
    on_first = float(brightness[first].mean()) if first.any() else -np.inf
    on_second = float(brightness[~first].mean()) if (~first).any() else -np.inf
    intercept = lines['first_intercept'] if on_first >= on_second else lines['second_intercept']
    return float(intercept % total)


def measure_pad(anchor: float, total: int, start: int, kept: int) -> int:
    """
    How many zeros to push in at the front so switch 0's echo lands where the sequence puts it.

    `start + kept // 2` is ramp + npoints_per_switch/2, the middle of the gradient plateau.  It is
    taken from discard_pre as MRStomrd2 writes it - the ramp - and not from the centred
    (total - kept)//2 an older reading of this file assumed; on cirrhrat those differ by 4.

    Reduced the short way round the switch, since a position is cyclic within one: an anchor at 26 of
    28 reaching a target of 2 has moved +4, not -24.
    """
    expected = start + kept // 2
    return int(round((expected - anchor + total / 2) % total - total / 2))


# ---------- applying it -----------------------------------------------------------------------


def prepend_zeros(data: np.ndarray, pad: int) -> np.ndarray:
    """
    Move the whole readout `pad` samples later by pushing zeros in at its front.

    Not a roll.  The readout is one continuous time series, so what falls off the end of the train is
    gone rather than wrapped back to its start, and the positions switch 0 gains at its front were
    never measured - zero says that, where a wrapped sample would claim the tail of the train was
    acquired before its head.  A negative pad drops that many leading samples and zero-fills the tail
    instead, which is the same statement made at the other end.
    Args:
        - data: (coils, samples)
        - pad: samples to push in at the front, negative to drop them from it
    """
    out = np.zeros_like(data)
    if pad > 0:
        out[:, pad:] = data[:, :data.shape[1] - pad]
    elif pad < 0:
        out[:, :data.shape[1] + pad] = data[:, -pad:]
    else:
        out[:] = data
    return out


def _phase_shift(block: np.ndarray, offset: float, axis: int = 0) -> np.ndarray:
    """One exact cyclic displacement by the shift theorem, along `axis`."""
    n = block.shape[axis]
    kernel = np.fft.fftfreq(n) * n
    shape = [1] * block.ndim
    shape[axis] = n
    spectrum = np.fft.fft(block, axis=axis)
    return np.fft.ifft(spectrum * np.exp(-2j * np.pi * kernel.reshape(shape) * offset / n),
                       axis=axis)


def apply_switch(data: np.ndarray, offsets: np.ndarray, nswitch: int, total: int) -> np.ndarray:
    """
    Displace each switch over the whole switch, exactly, so the window ends up holding different
    samples.

    That is what re-centres the k-space line on the echo and fixes the asymmetric truncation which
    broadens the spatial point spread - it sharpens the metabolite map, where a phase confined to the
    window cannot.  Cyclic: what a switch's own displacement pulls in at the end furthest from the
    anchor is that switch's own ramp and rephasing points, not a real neighbour, and the weight this
    implicitly puts on them never reaches zero the way a truncated sinc does - see _phase_shift's
    impulse response.  `roll` and `alloc` both read the real neighbour instead; see apply_roll and
    apply_alloc.
    """
    used = nswitch * total
    out = np.array(data, copy=True).astype(np.complex128)
    body = out[:, :used].reshape(out.shape[0], nswitch, total)
    for i, offset in enumerate(offsets):
        body[:, i, :] = _phase_shift(body[:, i, :], float(offset), axis=1)
    out[:, :used] = body.reshape(out.shape[0], used)
    return out


def apply_roll(data: np.ndarray, offsets: np.ndarray, nswitch: int, total: int) -> np.ndarray:
    """
    Displace each switch by the nearest whole sample, read from the real contiguous readout.

    `roll` throws away the same fractional part `shift` corrects for - the question it exists to
    answer is what that fraction is worth - but there is no reason to also throw away where the
    whole-sample part of the shift points.  `np.roll` wraps within one switch, so what fills the far
    end of a large displacement is that switch's own ramp and rephasing samples: junk, standing in
    for a real neighbour that the readout actually has.  This reads that neighbour instead, at
    whatever integer position `switch_offsets` rounds to, over the same contiguous coordinate
    `apply_contiguous` uses - it differs from it only in reading one real sample per position rather
    than interpolating sixteen.
    Args:
        - offsets: per-switch displacement in samples, from switch_offsets; rounded here to the
          nearest whole sample before it addresses anything
    """
    used = nswitch * total
    whole = np.rint(offsets).astype(int)
    source = (np.arange(used).reshape(nswitch, total) - whole[:, None]).ravel()
    inside = (source >= 0) & (source < used)
    out = np.array(data, copy=True).astype(np.complex128)
    moved = np.zeros((data.shape[0], used), dtype=np.complex128)
    moved[:, inside] = out[:, source[inside]]
    out[:, :used] = moved
    return out


def apply_alloc(data: np.ndarray, offsets: np.ndarray, nswitch: int, total: int) -> np.ndarray:
    """
    Displace each switch by a linear blend of its two nearest real neighbours.

    The two-tap answer to what roll and shift answer more crudely and more precisely: split the
    fractional part of the displacement between the two samples straddling the source position, in
    proportion to how close each one is - a source of 1142.748 reads mostly the sample at 1143
    (weight 0.748, since it is only 0.252 away) and a little of 1142 (weight 0.252).  Read from the
    contiguous readout, over the same source coordinate apply_contiguous and apply_roll use, so a
    large displacement pulls in a real neighbouring switch rather than wrapping - unlike the cyclic
    two-tap blend this replaces, which shaded the field of view instead of blurring it for the same
    reason `shift` and `roll` used to wrap: the two taps came from within the same switch.
    """
    used = nswitch * total
    source = (np.arange(used, dtype=float).reshape(nswitch, total)
              - np.asarray(offsets)[:, None]).ravel()
    base = np.floor(source).astype(int)
    frac = source - base
    out = np.array(data, copy=True).astype(np.complex128)
    moved = np.zeros((data.shape[0], used), dtype=np.complex128)
    for index, weight in ((base, 1.0 - frac), (base + 1, frac)):
        inside = (index >= 0) & (index < used)
        moved[:, inside] += weight[inside] * out[:, index[inside]]
    out[:, :used] = moved
    return out


def resample(line: np.ndarray, source: np.ndarray, half_width: int = 8) -> np.ndarray:
    """
    Read a readout at arbitrary fractional sample positions, band-limited.

    The readout is one continuous time series, so a sample between two of its points is recoverable
    by interpolation rather than by rounding.  A Lanczos-windowed sinc, so truncating the kernel does
    not ring, and zero outside the readout - at the two physical ends there is no neighbour, and zero
    says "not measured" where clamping would repeat an edge.
    Args:
        - line: (coils, samples)
        - source: float position to read for each output sample, same length as the readout
    """
    used = source.size
    base = np.floor(source).astype(int)
    out = np.zeros((line.shape[0], used), dtype=np.complex128)
    for tap in range(-half_width + 1, half_width + 1):
        index = base + tap
        frac = source - index
        weight = np.sinc(frac) * np.sinc(frac / half_width)
        inside = (index >= 0) & (index < line.shape[1])
        if inside.any():
            out[:, inside] += weight[inside] * line[:, index[inside]]
    return out


def apply_regrid(data: np.ndarray, slope: float, nswitch: int, total: int,
                 lead: float = 0.0) -> np.ndarray:
    """
    Resample the whole readout onto the period it was really acquired at.

    The other methods displace each switch; this re-grids once.  The defect is that the stream is
    divided by `total` where the gradient period is `total + slope`, so switch i of the real sequence
    begins at sample `i*(total + slope)` rather than `i*total`.  Reading there and writing to `i*total`
    puts every switch back on the grid the reconstruction assumes, with no per-switch discontinuity
    at the boundaries and nothing wrapped.

    Anchored on switch 0, like switch_offsets and for the same reason: reading switch i straight from
    `i * (total + slope)` leaves every echo at the within-period offset switch 0 already had, which
    is the position the zero-pad base has put where the sequence asks for it.  `lead` moves that
    grid, and the base having already placed the readout there, it is left at zero.
    """
    used = nswitch * total
    out = np.array(data, copy=True).astype(np.complex128)
    position = np.arange(total, dtype=float)
    source = np.concatenate([i * (total + slope) + position + lead
                             for i in range(nswitch)])
    out[:, :used] = resample(out[:, :used], source)
    return out


def apply_contiguous(data: np.ndarray, offsets: np.ndarray, nswitch: int, total: int,
                     half_width: int = 8) -> np.ndarray:
    """
    Displace each switch along the contiguous readout instead of wrapping inside it.

    `shift` is cyclic within a switch, so the samples it pulls in at the ends of the train - where
    the displacement is largest - are that switch's own ramp and rephasing points.  Those are junk,
    and burying the weak metabolites under them is what costs the whole-switch family its hydrate.
    The readout is one continuous time series, so the samples that really sit beside the window are
    the neighbouring switch's, and this takes those: band-limited interpolation with a windowed sinc
    over the contiguous readout, zero only at the two physical ends where there is no neighbour.
    Args:
        - half_width: taps either side; 8 puts the sinc truncation well below the noise
    """
    used = nswitch * total
    out = np.array(data, copy=True).astype(np.complex128)
    source = (np.arange(used, dtype=float).reshape(nswitch, total)
              - np.asarray(offsets)[:, None]).ravel()
    out[:, :used] = resample(out[:, :used], source, half_width)
    return out


def apply_base_and_method(data: np.ndarray, pad: int, offsets: np.ndarray, slope: float,
                          nswitch: int, total: int, method: str) -> np.ndarray:
    """
    The zero-pad base and then one drift method, in that order, over one readout.

    The order is the whole point: the pad puts switch 0's echo where the sequence asks for it, and
    the method straightens the walk away from there.  Reversed, the method would be straightening a
    train about a mean the pad then moves, and the two would no longer compose.

    Returned rather than written back, so the same call serves the stream correction and the -p
    figure, which tries every method on a copy and must leave the stream alone.
    """
    padded = prepend_zeros(data, pad) if pad else data
    if method == "regrid":
        return apply_regrid(padded, slope, nswitch, total)
    if method == "roll":
        return apply_roll(padded, offsets, nswitch, total)
    if method == "alloc":
        return apply_alloc(padded, offsets, nswitch, total)
    if method == "contiguous":
        return apply_contiguous(padded, offsets, nswitch, total)
    return apply_switch(padded, offsets, nswitch, total)


def record(header: mrd.Header, slope: float, pad: int, method: str) -> None:
    """Note what was taken out and how, so a second pass refuses rather than correcting twice."""
    if header.user_parameters is None:
        header.user_parameters = mrd.UserParametersType()
    for param in header.user_parameters.user_parameter_double:
        if param.name == SLOPE_PARAMETER:
            param.value = float(param.value) + slope
            break
    else:
        header.user_parameters.user_parameter_double.append(
            mrd.UserParameterDoubleType(name=SLOPE_PARAMETER, value=float(slope)))
    for param in header.user_parameters.user_parameter_long:
        if param.name == PAD_PARAMETER:
            param.value = int(param.value) + pad
            break
    else:
        header.user_parameters.user_parameter_long.append(
            mrd.UserParameterLongType(name=PAD_PARAMETER, value=int(pad)))
    for param in header.user_parameters.user_parameter_string:
        if param.name == METHOD_PARAMETER:
            param.value = method
            return
    header.user_parameters.user_parameter_string.append(
        mrd.UserParameterStringType(name=METHOD_PARAMETER, value=method))


def header_double(header: mrd.Header, name: str) -> Optional[float]:
    user = getattr(header, "user_parameters", None)
    for param in (getattr(user, "user_parameter_double", None) or []):
        if param.name == name:
            return float(param.value)
    return None


def pooled_lines(acqs: Sequence[mrd.Acquisition]) -> List[np.ndarray]:
    """
    Each acquisition's readout as one coil-summed (1, samples) copy.

    Copied, and taken before anything is applied to the stream: the correction writes back into
    acq.data in place, so a figure built from the acquisitions afterwards would be drawing a second
    correction on top of the first.
    """
    lines = []
    for acq in acqs:
        samples = np.asarray(acq.data)
        # MRS acquires on one channel; summed rather than assumed, as acquisition_cube does
        line = samples[0:1] if samples.shape[0] == 1 else samples.sum(axis=0, keepdims=True)
        lines.append(np.array(line, copy=True))
    return lines


def pooled_signal(lines: Sequence[np.ndarray], nswitch: int, total: int,
                  transform=None) -> np.ndarray:
    """
    A (switch, position) magnitude summed over every pooled readout, after `transform`.

    The same aggregate acquisition_cube's callers take, built from the lines pooled_lines put aside
    rather than through the (view, repetition) cube - which is what lets a method be tried on a copy
    without writing anything back into the stream.
    Args:
        - lines: (1, samples) readouts, from pooled_lines
        - transform: (coils, samples) -> (coils, samples), applied per readout; identity if None
    """
    used = nswitch * total
    signal = np.zeros((nswitch, total))
    for line in lines:
        moved = transform(line) if transform is not None else line
        signal += np.abs(moved[0, :used]).reshape(nswitch, total)
    return signal


def plot_method_matrix(panels: Sequence[Tuple[str, np.ndarray]], lines: dict, nswitch: int,
                       total: int, start: int, kept: int, label: str, path: Path) -> None:
    """
    Draw where the echo sits in every switch: before the correction, and under every method.

    The picture the whole correction is about: position within a switch across, switch number up, so
    an echo that walks is a slanted stripe and a corrected one is vertical.  Every panel is drawn the
    same way so they can be read against each other - the window the reconstruction reads is shaded,
    and the position the sequence puts the echo at is a dotted line, which is what the pad aims
    switch 0 at and what a working method holds the rest of the train on.

    The fitted lines go on the first panel only, since that is the data they were fitted to.
    Args:
        - panels: ordered (title, signal) pairs, each signal a (switch, position) magnitude
        - lines: what fit_peak_lines returned for the first panel
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    columns = min(len(panels), 4)
    rows = -(-len(panels) // columns)
    figure, axes = plt.subplots(rows, columns, figsize=(4.7 * columns, 5.8 * rows), squeeze=False)
    switches = np.arange(nswitch)
    expected = start + kept // 2
    for index, (title, signal) in enumerate(panels):
        ax = axes[index // columns][index % columns]
        ax.imshow(signal, aspect='auto', origin='lower', interpolation='nearest',
                  extent=(-0.5, total - 0.5, -0.5, nswitch - 0.5))
        ax.axvspan(start - 0.5, start + kept - 0.5, color='w', alpha=0.16,
                   label=f"window {start}..{start + kept - 1}")
        ax.axvline(expected, color='C2', lw=1.1, ls=':', label=f"expected {expected}")
        peaks = np.argmax(signal, axis=1)
        ax.plot(peaks, switches, 'x', color='C3', ms=4, label='brightest position')
        if index == 0:
            for intercept, style in ((lines['first_intercept'], '-'),
                                     (lines['second_intercept'], '--')):
                ax.plot((intercept + lines['slope'] * switches) % total, switches, style,
                        color='w', lw=1.2, alpha=0.8)
            ax.plot([], [], '-', color='w', lw=1.2, label='fitted echo lines')
        ax.set_xlabel(f"position within the {total} point switch")
        ax.set_ylabel("switch")
        ax.set_title(title)
        ax.legend(fontsize=7, loc='upper right')
    for blank in range(len(panels), rows * columns):
        axes[blank // columns][blank % columns].axis('off')
    figure.suptitle(label)
    figure.tight_layout()
    figure.savefig(path, dpi=140, bbox_inches='tight')
    plt.close(figure)
    print(f"  wrote {path}", file=sys.stderr)


# ---------- the raw side, for MRStomrd2 -w ------------------------------------------------------


def check_peak_position(group) -> int:
    """
    Plot where the two echoes sit in each switch, one figure per raw .MRD of one scan.

    The same two-line fit the stream side uses, read off the files before conversion.  One figure per
    file and never pooled, because pooling is what hides a file whose echo sits somewhere the others'
    does not; the repetitions inside one file are aggregated for the opposite reason, since one
    repetition of one view is mostly noise.
    Args:
        - group: a ScanGroup from MRSorganize; only its real data is read
    Returns:
        - how many files were plotted
    """
    # imported here so conversion needs neither the reader nor a plotting stack
    import matplotlib.pyplot as plt
    from MRSreader import MRSdata

    plotted = 0
    for filepath in group.rawdata_file_list:
        mrs = MRSdata()
        mrs.read_from_file(filepath)
        nswitch = max(mrs.nswitches, 1)
        if nswitch <= 1 or mrs.rawdata is None or mrs.rawdata.size == 0:
            print(f"Skipping {filepath}: no EPSI readout to find echoes in", file=sys.stderr)
            continue
        total = mrs.nsamples // nswitch
        kept = mrs.npoints_per_switch or total
        ramp = int(mrs.tramp // (mrs.sample_period // 10)) if mrs.sample_period else 0
        # magnitude summed over every axis but the samples one: views, slices, echoes and
        # repetitions are all just repeats of the same readout for finding an echo
        aggregated = np.abs(mrs.rawdata).sum(axis=tuple(range(1, mrs.rawdata.ndim)))
        signal = aggregated[:nswitch * total].reshape(nswitch, total)

        lines = fit_peak_lines(signal)
        switches = np.arange(nswitch)
        print(f"\n{filepath}: {nswitch} switches of {total}, {mrs.nrepetitions} repetition(s)",
              file=sys.stderr)
        print(f"  slope {lines['slope']:+.4f} per switch, a period of {total + lines['slope']:.4f}; "
              f"lines at {lines['first_intercept'] % total:.2f} ({lines['on_first']} switches, rms "
              f"{lines['first_rms']:.2f}) and {lines['second_intercept'] % total:.2f} "
              f"({lines['on_second']}, rms {lines['second_rms']:.2f})", file=sys.stderr)

        figure, axes = plt.subplots(figsize=(7, 9))
        axes.imshow(signal, aspect='auto', origin='lower', interpolation='nearest',
                    extent=(-0.5, total - 0.5, -0.5, nswitch - 0.5))
        axes.axvspan(ramp - 0.5, ramp + kept - 0.5, color='w', alpha=0.15,
                     label=f"window {ramp}..{ramp + kept - 1}")
        axes.plot(lines['peaks'][lines['first_family']], switches[lines['first_family']],
                  'x', color='C3', ms=6, label=f"first echo ({lines['on_first']})")
        if not lines['first_family'].all():
            axes.plot(lines['peaks'][~lines['first_family']], switches[~lines['first_family']],
                      '+', color='C1', ms=7, label=f"second echo ({lines['on_second']})")
        for intercept in (lines['first_intercept'], lines['second_intercept']):
            axes.plot((intercept + lines['slope'] * switches) % total, switches, '.',
                      color='w', ms=2, alpha=0.7)
        axes.set_xlabel(f"position within the {total} point switch")
        axes.set_ylabel("switch")
        axes.set_title(f"{Path(filepath).name}: {lines['slope']:+.4f} per switch")
        axes.legend(fontsize=8, loc='upper right')
        figure.tight_layout()
        plt.show()
        plt.close(figure)
        plotted += 1
    return plotted


# ---------- driver ----------------------------------------------------------------------------


def correct_stream(input_path: str, output_path: Optional[str], method: str = DEFAULT_METHOD,
                   force: bool = False, drop_first: bool = False,
                   plot: Optional[str] = None, pad: Optional[int] = None) -> int:
    """
    Measure the zero-pad base and the drift of a converted stream, and write it back corrected.

    The acquisitions are grouped by switch layout, because a file holds the series and the averaged
    prescan beside it and their layouts differ.  The biggest group's real data is what the reported
    numbers come from - pooling an averaged prescan in with the repetitions would measure neither.

    The prescan is never shifted.  It is calibration rather than data, and the scaling the
    reconstruction takes from it has to mean the same thing before and after a correction.  Any other
    layout is measured on its own data rather than having the series' numbers scaled onto it, and one
    whose echo cannot be fitted, or whose own geometry does not account for its switch, is left alone
    too - being left alone is the better failure.

    With no output path this is a dry run: -p still draws the figure, and nothing is written.
    Returns:
        - 0 when the stream was handled, 1 when the input carried nothing to correct
    """
    header, items = read_stream(input_path)
    acqs = [item.value for item in items if isinstance(item, mrd.StreamItem.Acquisition)]
    if not acqs:
        print(f"no acquisitions in {input_path}", file=sys.stderr)
        return 1
    if header.measurement_information.sequence_name == "epsigre43_FB_13C":
        print(f"{input_path} EVO1 data so there is no drift to take out",
              file=sys.stderr)
        if output_path:
            write_stream(output_path, header, items)
        return 0
    if (header_long(header, "nswitches") or 1) <= 1:
        print(f"{input_path} records no switch train, so there is no drift to take out",
              file=sys.stderr)
        if output_path:
            write_stream(output_path, header, items)
        return 0

    groups: Dict[Tuple[int, int, int, int], List[mrd.Acquisition]] = {}
    for acq in acqs:
        groups.setdefault(switch_layout(header, acq), []).append(acq)
    biggest = max(groups, key=lambda key: len(groups[key]))
    nswitch, total, kept, start = biggest
    pool = [a for a in groups[biggest] if not is_prescan(a)] or groups[biggest]

    cube, missing = acquisition_cube(pool, nswitch, total)
    if missing:
        print(f"  {missing} view/repetition slot(s) no acquisition filled, read as zero",
              file=sys.stderr)
    raw_signal = np.abs(cube).sum(axis=(2, 3))
    # put the readouts aside before the correction writes back into them, so -p draws every method
    # over the raw data rather than over whichever one was already applied to the stream
    raw_lines = pooled_lines(pool) if plot else []
    drift = measure_drift(cube, nswitch, total)
    anchor = readout_anchor(drift['lines'], raw_signal, total)
    measured_pad = measure_pad(anchor, total, start, kept)
    applied_pad = measured_pad if pad is None else int(pad)
    expected = start + kept // 2

    print(f"{nswitch} switches of {total} samples, reading {start}..{start + kept - 1}; "
          f"{len(pool)} acquisitions pooled", file=sys.stderr)
    print(f"  {'drift' if drift['usable'] else 'no usable drift'}: {drift['reason']}",
          file=sys.stderr)
    print(f"  slope {drift['slope']:+.4f} per switch, a period of {drift['period']:.4f} "
          f"where nsamples/{nswitch} records {total}", file=sys.stderr)
    # the measurement beside the `ramp + 3` an older converter hard-coded per tramp, which is what it
    # should reproduce: +7 on the tramp 112 cirrhrat data, +5 on the tramp 100 kidney data
    print(f"  switch 0's readout echo sits at {anchor:.2f} and the sequence puts it at {expected} "
          f"(ramp {start} + {kept}/2), so the base is {measured_pad:+d} zeros against the "
          f"{start + 3:+d} of ramp+3"
          + (f"; --pad {applied_pad:+d} overrides it" if pad is not None else ""), file=sys.stderr)

    applied = drift['slope'] if drift['usable'] else 0.0
    already = header_double(header, SLOPE_PARAMETER)
    corrected_signal = raw_signal
    if not drift['usable'] and not applied_pad:
        print("  nothing to take out, copying the stream through unchanged", file=sys.stderr)
    elif already is not None and not force:
        print(f"  this stream already records {already:+.4f} taken out of it, so it is copied "
              f"through unchanged; pass --force to correct it again", file=sys.stderr)
    else:
        print(f"  padding {applied_pad:+d} zeros and taking {applied:+.4f} per switch out of the "
              f"{total} point layout, --method {method}", file=sys.stderr)
        corrected = 0
        skipped_prescan = 0
        for key, group in groups.items():
            n, t, k, s = key
            # The prescan is calibration, not data: it is a separate acquisition that happens to
            # travel in the same file, the drift measured off the series says nothing about it, and
            # the scaling the reconstruction takes from it has to mean the same thing before and
            # after a correction.  Left alone on cirrhrat_43_1 it reads 625.678 whatever --method
            # ran; corrected along with the series it read 284.540 to 599.589, a calibration moving
            # with the thing it is supposed to calibrate.
            real = [a for a in group if not is_prescan(a)]
            skipped_prescan += len(group) - len(real)
            if not real:
                continue
            group = real
            if key == biggest:
                scaled, pad_here = applied, applied_pad
            else:
                # Measured on this layout's own data rather than scaled off the biggest group's.
                # The drift is a fractional timebase disagreement, so scaling by switch length is
                # the right shape of answer, but only where the layout agrees: on cirrhrat_43_1 the
                # scaled +0.2352 fits the prescan's own peaks worse than its own +0.0937 does, and
                # both fit it badly.  A layout that cannot measure its own drift is left alone
                # rather than corrected by extrapolation - being left alone is the better failure.
                own_cube, _ = acquisition_cube(group, n, t)
                own_signal = np.abs(own_cube).sum(axis=(2, 3))
                own = measure_drift(own_cube, n, t)
                scaled = own['slope'] if own['usable'] else 0.0
                # the pad needs `start + kept // 2`, which means nothing when the layout's own
                # numbers do not account for its switch
                pad_here = (measure_pad(readout_anchor(own['lines'], own_signal, t), t, s, k)
                            if geometry_closes(t, k, s) else 0)
                if not own['usable'] or not pad_here:
                    print(f"    the {t} point layout is left as it is: "
                          + ("its geometry does not close, "
                             f"{s} + {k} + 3*{s} is {4 * s + k} of {t}, so where its echo belongs "
                             "cannot be worked out" if not geometry_closes(t, k, s)
                             else f"no usable drift, {own['reason']}"), file=sys.stderr)
                if not own['usable'] and not pad_here:
                    continue
                print(f"    the {t} point layout is corrected on its own measurement, "
                      f"{scaled:+.4f} per switch on a {pad_here:+d} base", file=sys.stderr)
            offsets = switch_offsets(n, scaled)
            for acq in group:
                data = np.asarray(acq.data)
                moved = apply_base_and_method(data, pad_here, offsets, scaled, n, t, method)
                acq.data = moved.astype(data.dtype)
                if drop_first:
                    # the first FID point is the integral of the spectrum, which a sum of
                    # Lorentzians plus a constant baseline fits badly; the legacy reconstruction
                    # drops it and that alone accounts for much of its cleaner fit
                    zeroed = np.array(acq.data, copy=True)
                    zeroed[:, :t] = 0
                    acq.data = zeroed
                corrected += 1
        record(header, applied, applied_pad, method)
        print(f"  corrected {corrected} acquisition(s)"
              + (f", leaving {skipped_prescan} prescan acquisition(s) as they are"
                 if skipped_prescan else ""), file=sys.stderr)
        corrected_cube, _ = acquisition_cube([a for a in groups[biggest] if not is_prescan(a)]
                                             or groups[biggest], nswitch, total)
        corrected_signal = np.abs(corrected_cube).sum(axis=(2, 3))

    # Where the echo ends up, against the window the reconstruction reads.  With the pad applied
    # these should agree: the pad is what puts the echo on the position the sequence asks for, and
    # the method is what holds the rest of the train there, so a warning here is now a real finding
    # rather than the standing complaint it was when nothing moved the readout at all.
    landed = int(np.argmax(corrected_signal.sum(axis=0)))
    inside = start <= landed <= start + kept - 1
    print(f"  the echo lands at {landed}, against the {expected} the sequence puts it at "
          f"(ramp {start} + {kept}/2) and the {start}..{start + kept - 1} the reconstruction reads",
          file=sys.stderr)
    if not inside:
        print(f"WARNING the echo sits OUTSIDE the window the reconstruction reads, so the drift is "
              f"straight but the readout is in the wrong place; a base of "
              f"{applied_pad + expected - landed:+d} rather than {applied_pad:+d} would centre it",
              file=sys.stderr)

    if plot:
        # every method over the same padded base, each on a copy, so the figure compares them
        # against one another and against the raw readout without touching the stream.  The base on
        # its own is not drawn: it translates the whole picture by `pad` and leaves the walk exactly
        # as it was, which the raw panel already shows
        offsets = switch_offsets(nswitch, applied)
        panels = [("raw", raw_signal)]
        for name in METHODS:
            panels.append((name, pooled_signal(
                raw_lines, nswitch, total,
                lambda line, m=name: apply_base_and_method(
                    line, applied_pad, offsets, applied, nswitch, total, m))))
        plot_method_matrix(panels, drift['lines'], nswitch, total, start, kept,
                           f"{Path(input_path).name}: {drift['slope']:+.4f} per switch, "
                           f"base {applied_pad:+d}", Path(plot))

    if output_path:
        write_stream(output_path, header, items)
        print(f"wrote {output_path}", file=sys.stderr)
    else:
        print("  no output path, so nothing was written", file=sys.stderr)
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Take the EPSI echo drift out of a converted MRD2 stream")
    parser.add_argument("-i", "--input", default=os.environ.get("INPUT_PIPE"),
                        help="the .mrd2 to correct, a FIFO carrying one, or - for stdin "
                             "(default: $INPUT_PIPE)")
    parser.add_argument("-o", "--output", default=os.environ.get("OUTPUT_PIPE"),
                        help="where to write, or - for stdout (default: $OUTPUT_PIPE)")
    parser.add_argument("--method", choices=METHODS, default=DEFAULT_METHOD,
                        help=f"where to apply the displacement (default: {DEFAULT_METHOD})")
    parser.add_argument("--pad", type=int, default=None, metavar="SAMPLES",
                        help="prepend this many zeros instead of the measured base, moving the "
                             "whole readout that many samples later; 0 turns the base off. The "
                             "measurement reproduces the ramp+3 an older converter hard-coded "
                             "(cirrhrat 7, ischemia 5)")
    parser.add_argument("--drop-first-switch", action="store_true",
                        help="zero the first switch, whose sample is the integral of the spectrum "
                             "and fits a sum of Lorentzians badly")
    parser.add_argument("-p", "--plot", metavar="PNG",
                        help="write a figure of the echo position per switch: the raw readout "
                             "and every method over the zero-pad base. Without -o this is a dry "
                             "run and no stream is written")
    parser.add_argument("--force", action="store_true",
                        help="correct a stream that already records a drift taken out of it")
    args = parser.parse_args()

    if not args.input:
        parser.error("--input is required when $INPUT_PIPE is unset")
    if not args.output and not args.plot:
        parser.error("--output is required when $OUTPUT_PIPE is unset, unless --plot is given")
    # neither a FIFO the sidecar has not created yet nor stdin is a path that exists
    if args.input not in ("-", os.environ.get("INPUT_PIPE")) and not Path(args.input).exists():
        parser.error(f"{args.input} does not exist")

    try:
        return correct_stream(args.input, args.output, args.method, args.force,
                              args.drop_first_switch, args.plot, args.pad)
    except RuntimeError as failure:
        print(f"Cannot read {args.input}: {failure}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
