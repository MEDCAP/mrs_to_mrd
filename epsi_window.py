"""
Report where the usable sampling window sits inside an EPSI gradient switch, and whether the echo
drifts along the switch train, for every EPSI file the input resolves to. Converts nothing - a
diagnostic companion to MRStomrd2, sharing its EPSI switch layout (is_epsi, switch_layout) so the
window this reports and the window a real conversion writes are the same arithmetic.

    -t/--tar     tar archive of one scan directory, or a FIFO carrying one
    -i/--input   a single .MRD file
    -f/--folder  directory to walk for .MRD files

The window is reported one block per file, because the echo can sit at a different position in one
scan directory than in the next and averaging that over a whole experiment is exactly what would
hide it. The drift along the switch train is reported once for all of them together, because that
is the granularity it can be measured at: a single repetition does not carry the signal to place
the echo per switch, and the drift is the same in every repetition anyway. That is why grouping is
deliberately not run here: every file the input resolves to is read and reported on independently,
which is what lets an echo that moves between scan directories show up.

This command corrects nothing, but the correction lives in this module: roll_switches and
shift_rawdata take a measured drift out of a readout, and MRStomrd2 -w is what walks the datasets of
an experiment and applies them. The measurement, the roll and the plots stay here as module
functions so that both sides - this report per file, MRStomrd2 per scan - work from one arithmetic.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Tuple

import numpy as np

import mrs_organize
from MRSreader import MRSdata
from MRStomrd2 import is_epsi, switch_layout
from mrs_tar import read_scan_tar


# How wide a drift the search covers, in samples per switch. The kidney data runs from -0.03 on the
# one experiment that does not drift to +0.27 on the worst that does, so half a sample per switch,
# which sweeps the echo one and a half switch widths over a 64 switch train, leaves ample room. A
# winner sitting on this bound is reported as no measurement rather than as a drift, since it means
# either the real optimum is outside the range or there was no optimum and noise won
DRIFT_SLOPE_LIMIT = 0.5
DRIFT_COARSE_STEP = 0.02        # first pass over that range
DRIFT_FINE_STEP = 0.002         # second pass, one coarse step either side of the winner
# How much sharper removing the drift has to make the readout profile before the drift is believed,
# as a ratio of peak over median after to peak over median before. Measured across five experiments:
# 1.05 on the one with no drift, 1.85 and above on the four with one, so the threshold sits in a
# wide gap rather than on a cliff
DRIFT_MIN_GAIN = 1.5
# peak over median the corrected profile has to reach, so that a large ratio on a scan that had no
# echo to sharpen either way is refused rather than acted on
DRIFT_MIN_SNR = 5.0


# ---------- finding the EPSI readout sampling window ----------------------------------------


def switch_cube(mrs: MRSdata) -> np.ndarray:
    """
    The readout laid out as (switch, position in switch, view, everything else).

    The sample axis comes first in rawdata, so splitting it in C order is exactly the switch-major
    order the samples were acquired in. The slice, sliceview and echo axes are single valued in an
    EPSI scan and are folded in with the repetitions, since for every purpose here they are all just
    repeats. Samples past the last whole switch are dropped, which is the same truncation the
    conversion applies when nsamples does not divide by the switch count
    Args:
        - mrs: one parsed MRS file, read rather than probed
    Returns:
        - (nswitch, total, nviews, repeats) view of the raw data
    Raises:
        - ValueError when called on a file whose data block was never read
    """
    if mrs.rawdata is None:
        raise ValueError("no data to profile; read_from_file rather than probe_from_file")
    nswitch, total, _ = switch_layout(mrs)
    return mrs.rawdata[:nswitch * total].reshape(nswitch, total, mrs.nviews, -1)


def spectral_peak(cube: np.ndarray, view: Optional[int] = None):
    """
    How much signal sits at each position within a readout switch, over the switches given.

    Every sample of the readout falls at one of `total` positions inside a gradient switch, and the
    gradient echo, the point where kx crosses zero, sits at one of them. Summing raw magnitude over
    the switches would locate it only when the object fills the field of view: for a compact object
    |k-space| is close to flat and the echo hides in the phase instead. So this transforms along the
    switch axis and keeps the strongest spectral line, which is a coherent sum over the switches
    given and lifts the echo well clear of the noise. That coherent gain is why the drift
    measurement passes a block of switches rather than one at a time.
    Args:
        - cube: (switch, position, view, repeat), whole or a slice along the switch axis
        - view: phase encode line to read, or None for the one carrying the most signal. Not a knob
          to offer: drift_sharpness has to score every candidate slope on the same view, since
          re-picking the brightest one per candidate would leave the scores incomparable
    Returns:
        - (profile of length total, peak over median of the profile, the view it was read from)
    """
    spectrum = np.fft.fft(cube, axis=0)
    power = (np.abs(spectrum) ** 2).sum(axis=(1, 2, 3))
    band = np.abs(spectrum[int(np.argmax(power))])               # (position, view, repeat)

    if view is None:
        view = int(np.argmax(band.sum(axis=(0, 2))))
    profile = band[:, view, :].sum(axis=1)
    median = np.median(profile)
    return profile, (float(profile.max() / median) if median else np.inf), view


def switch_profile(mrs: MRSdata):
    """
    How much signal sits at each position within a readout switch, over the whole readout.
    Args:
        - mrs: one parsed MRS file, read rather than probed
    Returns:
        - (profile of length total, peak over median of the profile)
    Raises:
        - ValueError when called on a file whose data block was never read
    """
    profile, snr, _ = spectral_peak(switch_cube(mrs))
    return profile, snr


def sliding_window(profile, kept: int):
    """
    Score every candidate sampling window of `kept` consecutive positions.

    Windows wrap, because the last position of one switch is followed by the first position of the
    next, so a window may legitimately straddle the boundary.
    Args:
        - profile: signal per position within the switch, from switch_profile
        - kept: width of the window, i.e. the points the reconstruction keeps per switch
    Returns:
        - (signal summed inside each window, position of the profile peak within each window), both
          indexed by the window's first position
    """
    total = len(profile)
    offsets = np.arange(kept)
    signal = np.array([profile[(start + offsets) % total].sum() for start in range(total)])
    # where the echo lands inside the window; the window is right when this is its middle
    peak_at = (int(np.argmax(profile)) - np.arange(total)) % total
    return signal, peak_at


def window_report(mrs: MRSdata) -> dict:
    """
    Pick the sampling window, both ways, and translate it for the reconstruction.

    Two criteria, because they can disagree and the disagreement is the interesting part: the window
    holding the most signal, and the window whose middle sits on the echo. mrd2recon addresses the
    same window as an offset from discard_pre, so that conversion is reported alongside, as the
    --pad that would select it
    Returns:
        - dict of the profile, the two winning window starts and the pads they imply
    """
    nswitch, total, kept = switch_layout(mrs)
    profile, snr = switch_profile(mrs)
    signal, peak_at = sliding_window(profile, kept)
    middle = (kept - 1) / 2

    by_signal = int(np.argmax(signal))
    # distance from the middle measured the short way round, so a wrapped window is not penalised
    # for having its peak reported as position 27 rather than -1
    from_middle = np.abs((peak_at - middle + total / 2) % total - total / 2)
    by_centre = int(np.argmin(from_middle))
    # the same arithmetic generate_acquisition writes into the header, so the pads below address the
    # window the converted stream actually carries
    discard_pre = (total - kept) // 2

    def pad_for(start: int) -> int:
        """
        The mrd2recon --pad that makes it read this window.

        It addresses a window as discard_pre - pad, so the pad is that difference, wrapped the short
        way round the switch: the positions are cyclic, so a window starting at 29 of 34 is 5 before
        the boundary rather than 29 after it
        """
        return int((discard_pre - start + total // 2) % total - total // 2)

    return dict(profile=profile, signal=signal, peak_at=peak_at, snr=snr,
                nswitch=nswitch, total=total, kept=kept, discard_pre=discard_pre,
                peak=int(np.argmax(profile)),
                by_signal=by_signal, by_centre=by_centre,
                pad_by_signal=pad_for(by_signal),
                pad_by_centre=pad_for(by_centre))


def plot_switch_profile(mrs: MRSdata) -> dict:
    """
    Plot the sliding window scan over the positions within a readout switch.

    Three panels: the signal at each position with the chosen windows drawn on it, the sliding
    window scan itself, and the same profile per switch so that an echo whose position drifts along
    the echo train shows up rather than being averaged away. Opening the window blocks until it is
    closed, so a whole experiment is walked one file at a time.
    Returns:
        - the window_report dict, so a caller can act on the numbers it plotted
    """
    # imported here rather than at module scope, so reporting needs no plotting stack unless plotted
    import matplotlib.pyplot as plt

    report = window_report(mrs)
    total, kept = report['total'], report['kept']
    profile, signal = report['profile'], report['signal']
    middle = (kept - 1) / 2
    positions = np.arange(total)

    figure, axes = plt.subplots(3, 1, figsize=(11, 10))
    title = (f"{mrs.sequence_name or 'unknown sequence'}: {report['nswitch']} switches of "
             f"{total} points, {kept} kept, peak/median {report['snr']:.2f}")
    figure.suptitle(title)

    axes[0].plot(positions, profile, 'o-', color='C0')
    axes[0].axvline(report['peak'], color='C1', lw=2, label=f"echo peak at {report['peak']}")
    for start, pad, colour, name in ((report['by_signal'], report['pad_by_signal'], 'C2', 'most signal'),
                                     (report['by_centre'], report['pad_by_centre'], 'C3', 'echo centred')):
        # drawn as the positions it covers, so a window that wraps appears at both ends
        covered = (start + np.arange(kept)) % total
        axes[0].plot(covered, profile[covered], 'o', ms=11, mfc='none', color=colour,
                     label=f"{name}: start {start}, pad {pad}, "
                           f"echo at {report['peak_at'][start]} of {kept}")
    axes[0].set_xlabel(f"position within the {total} point switch")
    axes[0].set_ylabel("signal")
    axes[0].legend(fontsize=8)

    axes[1].plot(positions, signal / signal.max(), 'o-', color='C0', label='signal in window')
    axes[1].plot(positions, report['peak_at'] / total, 's--', ms=3, color='C4',
                 label='where the peak lands in the window')
    axes[1].axhline(middle / total, color='C3', ls=':', label='window middle')
    axes[1].axvline(report['by_signal'], color='C2', lw=2)
    axes[1].axvline(report['by_centre'], color='C3', lw=2)
    axes[1].set_xlabel("first position of the window")
    axes[1].set_ylabel("normalised")
    axes[1].legend(fontsize=8)

    cube = switch_cube(mrs)
    axes[2].imshow(np.abs(cube).sum(axis=(2, 3)), aspect='auto', origin='lower',
                   interpolation='nearest')
    axes[2].axvline(report['by_centre'], color='C3', lw=1.5)
    axes[2].axvline((report['by_centre'] + kept - 1) % total, color='C3', lw=1.5)
    axes[2].set_xlabel(f"position within the switch")
    axes[2].set_ylabel("switch")

    figure.tight_layout()
    plt.show()
    return report


# ---------- removing the echo drift along the switch train ----------------------------------


def pooled_switch_cube(mrs_list: Sequence[MRSdata]) -> np.ndarray:
    """
    Every file of one scan laid out as a single (switch, position, view, repeats) cube.

    The drift is a property of the gradient timing, so it is the same in every repetition, and one
    repetition on its own does not carry enough signal to measure it. Pooling is what lets the two
    ways a scan arrives be measured identically: an epsi experiment that puts one repetition in each
    subdirectory pools across files, one that puts all 25 in a single file pools within it, and
    switch_cube already folds the repetition axis into its last axis either way.
    Args:
        - mrs_list: the parsed files of one scan, all read rather than probed
    Returns:
        - the pooled cube, or None when there is nothing to pool
    """
    cubes = []
    for mrs in mrs_list:
        cube = switch_cube(mrs)
        if cubes and cube.shape[:3] != cubes[0].shape[:3]:
            print(f"WARNING leaving a {cube.shape[:3]} readout out of a {cubes[0].shape[:3]} drift "
                  f"measurement", file=sys.stderr)
            continue
        cubes.append(cube)
    if not cubes:
        return None
    return cubes[0] if len(cubes) == 1 else np.concatenate(cubes, axis=3)


def roll_switches(block: np.ndarray, shifts: Sequence[int]) -> np.ndarray:
    """
    Move each switch of a readout by its own shift.

    The only arrangement a drift fits in. Every other way of addressing the readout - discard_pre,
    a --pad, the ramp time a reconstruction derives its window from - names one offset for the whole
    readout, and a drift is precisely the case where one offset cannot be right in every switch.

    The loop runs over switches rather than over the positions inside one, since that is where the
    shift varies: within a switch it is a single number, and rolling the switch by it moves every
    position together. Positions are cyclic within a switch, so this wraps rather than padding -
    what leaves one end of a switch arrives at the other, which is where the neighbouring switch's
    samples sat anyway
    Args:
        - block: (nswitch, positions within a switch, ...), a switch_cube or a reshaped readout
        - shifts: integer shift per switch, from switch_shifts, positive moving samples later
    Returns:
        - a new array of the same shape, switch i rolled by shifts[i]
    """
    rolled = np.empty_like(block)
    for iswitch, shift in enumerate(shifts):
        rolled[iswitch] = np.roll(block[iswitch], int(shift), axis=0)
    return rolled


def shift_rawdata(mrs: MRSdata, shifts: Sequence[int]) -> None:
    """
    Take a measured drift out of one file's readout, in place.

    The raw side twin of mrd2shift.roll_acquisition, which does this to a converted stream. Working
    on rawdata rather than on a cube is what makes the correction outlive the measurement: every
    axis is kept, so the file goes on converting, plotting and reconstructing exactly as it did,
    only with its echoes lined up.

    Samples past the last whole switch are left where they are, the same truncation switch_cube and
    the conversion already apply: nsamples // nswitch leaves a remainder on some sequences, and
    those trailing points belong to no switch to be rolled with
    Args:
        - mrs: one parsed MRS file of the epsi family, read rather than probed
        - shifts: integer shift per switch, from switch_shifts
    Returns:
        - None; mrs.rawdata is rewritten
    Raises:
        - ValueError when called on a file whose data block was never read
    """
    if mrs.rawdata is None:
        raise ValueError("no data to shift; read_from_file rather than probe_from_file")
    nswitch, total, _ = switch_layout(mrs)
    used = nswitch * total
    tail = mrs.rawdata.shape[1:]
    # reshaped through the sample axis, which comes first, so this is the switch-major order the
    # samples were acquired in - the same split switch_cube reads the echo position from. Assigned
    # back through a slice rather than rolled in place, since the reshape may be a copy
    body = mrs.rawdata[:used].reshape((nswitch, total) + tail)
    mrs.rawdata[:used] = roll_switches(body, shifts).reshape((used,) + tail)


def drift_sharpness(cube: np.ndarray, nswitch: int, slope: float, view: Optional[int]) -> float:
    """
    How well defined the echo is once a drift of `slope` is taken out of the switch train.

    Scored by peak over median of the readout profile, which is what the echo drifting smears: every
    switch contributes to that profile, so an echo that sits at one position in all of them gives a
    sharp peak and one that walks across the switch gives something close to flat.

    Scored through the same roll that is finally applied, so the sharpening measure_echo_drift
    promises is the sharpening the corrected data has
    Args:
        - cube: (switch, position, view, repeats), as pooled_switch_cube returns
        - nswitch: switches in the readout, i.e. cube.shape[0]
        - slope: candidate drift in samples per switch
        - view: the view to profile, held fixed across candidates so they compare
    Returns:
        - peak over median of the profile after the candidate correction
    """
    rolled = roll_switches(cube, switch_shifts(nswitch, cube.shape[1], slope))
    return spectral_peak(rolled, view=view)[1]


def measure_echo_drift(cube: np.ndarray, nswitch: int, total: int) -> dict:
    """
    Whether the echo moves along the switch train, and by how much per switch.

    The reconstruction reads the same positions out of every switch of a readout: mrd2recon works
    from `iswitch * points_per_switch + discard_pre - pad`, whose offset within the switch does not
    depend on the switch. That is only right while the echo stays put, and on the kidney data it
    does not: it walks about a dozen positions of a twenty position switch across the train, out of
    the window the reconstruction keeps at one end and against its far edge at the other.

    Measured by searching the drift itself rather than by locating the echo and fitting a line
    through where it was found. Each candidate slope is applied and scored by how sharp it leaves the
    readout profile, and the best scoring one wins. Locating the echo first cannot work here, because
    every way of doing it needs a profile of the echo and a drift this large is what destroys that
    profile: summed over the whole train the echo smears across two thirds of the switch and its peak
    over median falls to 2.4, so a measurement that starts by requiring a clean profile refuses
    exactly the scans that need correcting. Accumulating integer lags between neighbouring blocks of
    switches avoids that but is biased by the rounding, reading +0.25 where the truth is +0.232.
    Searching the correction has neither problem, and its acceptance test is the question actually
    worth asking: did the echo get sharper.

    Coarse pass over the whole range, then a fine pass either side of the winner, which is 61 and 21
    candidates for the default settings and about a quarter of a second on a 64 by 20 by 12 by 25
    cube.
    Args:
        - cube: (switch, position, view, repeats), as pooled_switch_cube returns
        - nswitch, total: the switch layout the cube was built on, from switch_layout
    Returns:
        - dict of the measurement: usable, reason, the slope in samples per switch, the drift it
          comes to over the whole train, peak over median before and after and their ratio, whether
          the winner sat on the bound of the search, and the switch period the slope implies
    """
    # picked once and held for every candidate, so that the scores compare
    view = spectral_peak(cube)[2]
    base = drift_sharpness(cube, nswitch, 0.0, view)

    coarse = np.arange(-DRIFT_SLOPE_LIMIT, DRIFT_SLOPE_LIMIT + DRIFT_COARSE_STEP / 2,
                       DRIFT_COARSE_STEP)
    scores = [drift_sharpness(cube, nswitch, slope, view) for slope in coarse]
    around = coarse[int(np.argmax(scores))]
    # clipped to the range, so the fine pass never proposes a slope the coarse one was not allowed
    fine = np.arange(max(around - DRIFT_COARSE_STEP, -DRIFT_SLOPE_LIMIT),
                     min(around + DRIFT_COARSE_STEP, DRIFT_SLOPE_LIMIT) + DRIFT_FINE_STEP / 2,
                     DRIFT_FINE_STEP)
    scores = [drift_sharpness(cube, nswitch, slope, view) for slope in fine]
    best = int(np.argmax(scores))
    slope, snr = float(fine[best]), float(scores[best])

    drift = slope * (nswitch - 1)
    report = dict(usable=False, reason="", slope=slope, drift=drift, snr=snr, base_snr=base,
                  gain=(snr / base if base else 0.0), nswitch=nswitch, total=total, view=view,
                  edge=abs(slope) >= DRIFT_SLOPE_LIMIT - DRIFT_FINE_STEP / 2,
                  period=total + slope, reps=cube.shape[3])

    if report['edge']:
        report['reason'] = (f"the best drift found, {slope:+.3f} samples per switch, sits on the "
                            f"{DRIFT_SLOPE_LIMIT} bound of the search, so either the real one is "
                            f"outside that range or there is no echo here and noise won")
        return report
    if abs(drift) < 1.0:
        report['reason'] = (f"the best drift found, {drift:+.2f} samples over {nswitch} switches, is "
                            f"below the one sample a roll could move")
        return report
    if report['gain'] < DRIFT_MIN_GAIN:
        report['reason'] = (f"removing {drift:+.2f} samples over {nswitch} switches only sharpens "
                            f"the readout profile from {base:.2f} to {snr:.2f} peak/median, a "
                            f"factor of {report['gain']:.2f} against the {DRIFT_MIN_GAIN} it needs "
                            f"to be worth acting on")
        return report
    if snr < DRIFT_MIN_SNR:
        report['reason'] = (f"removing {drift:+.2f} samples over {nswitch} switches leaves a "
                            f"peak/median of only {snr:.2f}, below the {DRIFT_MIN_SNR} that says "
                            f"there was an echo to sharpen rather than noise to line up")
        return report

    report['usable'] = True
    report['reason'] = (f"{drift:+.2f} samples over {nswitch} switches, {slope:+.4f} per switch, "
                        f"sharpening the readout profile from {base:.2f} to {snr:.2f} peak/median, "
                        f"a factor of {report['gain']:.2f}")
    return report


def switch_shifts(nswitch: int, total: int, slope: float) -> np.ndarray:
    """
    How far to move each switch so the echo stops moving along the train.

    Measured against the middle of the train rather than against the sampling window, so this
    removes the drift and leaves the mean echo position exactly where it was. That keeps it
    orthogonal to where the window sits: --pad, the pads this reports and the ramp time mrd2recon
    derives its default from all go on meaning what they mean today, and a scan can be de-drifted
    without its window having to be found again.

    Shifts come back reduced the short way round the switch, since a position is cyclic within one
    Args:
        - nswitch, total: the switch layout, from switch_layout
        - slope: drift to remove, in samples per switch
    Returns:
        - integer shift per switch, length nswitch, positive meaning the samples move later
    """
    middle = (nswitch - 1) / 2
    shifts = np.rint(slope * (middle - np.arange(nswitch))).astype(int)
    return (shifts + total // 2) % total - total // 2


def measure_group_drift(mrs_list: Sequence[MRSdata]
                       ) -> Optional[Tuple[dict, Tuple[int, int, int]]]:
    """
    Measure the echo drift of one dataset, from every repetition of every file of it pooled.

    One call for the three steps a drift measurement always takes together - pool, read the switch
    layout, search the correction - so that a caller walking datasets stays a loop and every caller
    measures the same way. The layout comes back with the measurement because the shifts that take
    the drift out depend on it: the same slope lands on different integers in a 28 point switch than
    in a 20 point one.

    What goes into one call is the caller's decision and it is not a free one. The drift is a
    property of the gradient timing, so it is the same in every repetition of a scan and pooling
    them is what gives it enough signal to be found; an averaged prescan pooled in with 25
    repetitions of real data measures neither of them
    Args:
        - mrs_list: the parsed files of one dataset, read rather than probed
    Returns:
        - (the measurement from measure_echo_drift, (nswitch, total, kept)), or None when the list
          held no EPSI readout to pool
    """
    epsi = [mrs for mrs in mrs_list if is_epsi(mrs) and mrs.rawdata is not None]
    if not epsi:
        return None
    cube = pooled_switch_cube(epsi)
    if cube is None:
        return None
    nswitch, total, kept = switch_layout(epsi[0])
    return measure_echo_drift(cube, nswitch, total), (nswitch, total, kept)


def shift_report(drift: dict, layout: Tuple[int, int, int], label: str = "") -> None:
    """
    Print what a drift measurement found, and what taking it out would cost.

    The shift range is printed against the discarded ramp on purpose: rolling a switch brings in the
    samples of its neighbour, which are ramp points rather than signal, so once the shift exceeds
    the ramp the ends of the train have nothing valid left to move in. That is the number that says
    whether a measured drift can be corrected by moving whole samples at all
    Args:
        - drift: the measurement from measure_echo_drift
        - layout: (nswitch, total, kept), from measure_group_drift
        - label: what to call this dataset, e.g. a meas_id
    Returns:
        - None; everything goes to stdout beside the figures
    """
    nswitch, total, kept = layout
    print(f"\necho drift of {label or 'the pooled readout'}, {drift['reps']} repetitions pooled")
    print(f"  {'drift' if drift['usable'] else 'no usable drift'}: {drift['reason']}")
    print(f"  best slope {drift['slope']:+.4f} samples per switch, a switch period of "
          f"{drift['period']:.2f} samples where nsamples/{nswitch} records {total}")
    if drift['usable']:
        shifts = switch_shifts(nswitch, total, drift['slope'])
        discard = (total - kept) // 2
        print(f"  taking it out moves each switch by {int(shifts.min())} to {int(shifts.max())} "
              f"samples, against {discard} discarded ramp points either side")


def read_inputs(args) -> Iterable[Tuple[str, MRSdata]]:
    """
    Every .MRD the input resolves to, parsed, in acquisition order.

    Grouping is deliberately not run: this reports one file at a time, so nothing here needs to know
    which stream a file would end up in, and both mrs_organize.collect_mrd_paths and read_scan_tar
    already drop the AppleDouble shadows and order what is left the way conversion would read it
    Args:
        - args: the parsed command line, in exactly one of the three input modes
    Returns:
        - (name, parsed file) pairs, the name being a path or a tar member name
    """
    if args.tar:
        with open(args.tar, "rb") as tar_stream:
            _, spr_frequency, members = read_scan_tar(tar_stream)
        for name, payload in members:
            mrs = MRSdata()
            mrs.parse_from_buffer(payload)
            mrs.set_base_frequency(spr_frequency)   # the .MRD may defer its frequency to the sidecar
            yield name, mrs
    elif args.input:
        mrs = MRSdata()
        mrs.read_from_file(args.input)
        yield str(args.input), mrs
    else:
        for filepath in mrs_organize.collect_mrd_paths(args.folder):
            mrs = MRSdata()
            mrs.read_from_file(filepath)
            yield filepath, mrs


def report_windows(named_files: Iterable[Tuple[str, MRSdata]]) -> int:
    """
    Report the EPSI sampling window for every file, converting nothing.
    Args:
        - named_files: (name, parsed file) pairs, as read_inputs yields them
    Returns:
        - how many files were reported on
    """
    # materialised rather than streamed, because the drift is measured from all of them pooled
    epsi: List[Tuple[str, MRSdata]] = []
    for name, mrs in named_files:
        if not is_epsi(mrs):
            print(f"Skipping {name}: {mrs.sequence_name or 'unknown sequence'} is not an EPSI "
                  f"readout, so it has no gradient switches to place a window in", file=sys.stderr)
            continue
        if mrs.rawdata is None or mrs.rawdata.size == 0:
            print(f"Skipping {name}: no raw data was read", file=sys.stderr)
            continue
        epsi.append((name, mrs))
    if not epsi:
        return 0

    # this walks its input file by file and never groups it, so an experiment folder arrives here as
    # its data and its averaged phantom together. The averaged readouts are left out of the drift,
    # because pooling one of them in with the repetitions would measure something neither of them is.
    # Their windows are still reported below like any other file's
    data = [(name, mrs) for name, mrs in epsi if mrs.naverages <= 1]
    phantoms = len(epsi) - len(data)
    if phantoms:
        print(f"Leaving {phantoms} averaged file(s) out of the drift measurement, since pooling one "
              f"in with the repetitions would measure neither", file=sys.stderr)

    # one measurement for the whole input, from every repetition of every unaveraged file pooled
    measured = measure_group_drift([mrs for _, mrs in data]) if data else None
    if measured is None:
        print("\nno unaveraged EPSI readout to measure the echo drift from", file=sys.stderr)
    else:
        shift_report(*measured, label=f"{len(data)} file(s)")

    reported = 0
    for name, mrs in epsi:
        # one figure per file, each blocking until it is closed, so a whole experiment is walked
        # rather than averaged into one picture
        report = plot_switch_profile(mrs)
        print(f"\n{name}")
        print(f"  {report['nswitch']} switches of {report['total']} points, keeping "
              f"{report['kept']}, discard_pre={report['discard_pre']}, profile peak/median "
              f"{report['snr']:.2f}")
        print(f"  echo peak at position {report['peak']}")
        print(f"  most signal:  window starts at {report['by_signal']}, echo lands at "
              f"{report['peak_at'][report['by_signal']]} of {report['kept']}, "
              f"mrd2recon --pad {report['pad_by_signal']}")
        print(f"  echo centred: window starts at {report['by_centre']}, echo lands at "
              f"{report['peak_at'][report['by_centre']]} of {report['kept']}, "
              f"mrd2recon --pad {report['pad_by_centre']}")
        if report['snr'] < 1.5:
            print(f"  WARNING the profile of {name} is nearly flat, so this scan carries too little "
                  f"signal to place the window; check it against one that does", file=sys.stderr)
        reported += 1
    return reported


def main() -> int:
    """
    Report the EPSI sampling window and echo drift for what the input resolves to, converting
    nothing.
    Returns:
        - 0 when at least one file was reported on, 1 when nothing was
    """
    parser = argparse.ArgumentParser(
        description="Report the EPSI sampling window and echo drift for MR Solutions MRS data")
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("-t", "--tar", type=Path,
                      help="tar archive of one scan directory, or a FIFO carrying one")
    mode.add_argument("-i", "--input", type=Path,
                      help="single MRS .MRD file")
    mode.add_argument("-f", "--folder", type=Path,
                      help="directory to walk for MRS .MRD files")
    args = parser.parse_args()

    if args.tar and not args.tar.exists():
        parser.error(f"{args.tar} does not exist")
    if args.input and not args.input.is_file():
        parser.error(f"{args.input} is not a file")
    if args.folder and not args.folder.is_dir():
        parser.error(f"{args.folder} is not a directory")

    reported = report_windows(read_inputs(args))
    if not reported:
        print("No EPSI readout to report a sampling window for", file=sys.stderr)
        return 1
    print(f"Reported the sampling window of {reported} file(s)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
