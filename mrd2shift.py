"""
EVO2 epsi data shows a drift and shift of echo position, relative to the sampling window
To account for the correction, this script calculates the largest peak position and tries to 
align them on the expected peak position so that the mrd2recon stays relevant 

It can take a converted raw file as input.
    -i/--input   a single raw .MRD file
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Tuple

import numpy as np

from MRSreader import MRSdata
from MRSorganize import ScanGroup, organize_folder, read_scan_tar


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

# ndarray -> identify peaks
# ndarray, identified peaks -> plot
# mrd.acquisition -> ndarray
# groups -> ndarray

# helper function to call from MRStomrd2.py for local testing on scan_groups
def check_peak_position(scan_groups: ScanGroup) -> bool:
    """
    Check where the readout and rephasing echoes happen per switch for raw data, not prescan data

    Args:
        - scan_groups: the scan to walk, from mrs_organize. Only its real data is read; the averaged
          prescan beside it calibrates a reconstruction rather than being reconstructed
    Returns:
        - how many files were plotted
    """
    rawdata_across_rep = None
    for filepath in scan_groups.rawdata_file_list:
        mrs = MRSdata()
        mrs.read_from_file(filepath)
        # aggregate each repetition
        if rawdata_across_rep is None:
            rawdata_across_rep = mrs.rawdata
        else:
            rawdata_across_rep = np.concatenate((rawdata_across_rep, mrs.rawdata), axis=-1)
    # revise rawdata shape from (nsamples,...) to (nswitches, nsamples/nswitches...)
    total_points_per_switch = mrs.nsamples // mrs.nswitches
    per_switch_rawdata = rawdata_across_rep.reshape(((mrs.nswitches, total_points_per_switch) + rawdata_across_rep.shape[1:]))
    # sum across views and repetitions
    per_switch_rawdata = np.sum(np.abs(per_switch_rawdata),axis=(2,6)).squeeze()
     
    # ramp up - readout points - ramp down - ramp up rephase - ramp down rephase
    # ramp = mrs.tramp / mrs.sample_period 
    # cirrhrat sample_period=280, no_samples=1792, tramp=112, points_per_switch=12, nswitches=64
    # 112/28=4 - 12 - 4 - 4 - 4 = 1792/64 = 28
    # expected first peak = 4 + 7 (peak sits on the right of midpoint of 12 pts) = 11 or idx=10
    
    # units: sample_period in 10us and tramp in us
    ramp_points = mrs.tramp // (mrs.sample_period/10)
    ntotal_points_per_switch = mrs.nsamples // mrs.nswitches
    npoints_per_switch = mrs.npoints_per_switch
    # check the tramp and points_per_switch adds up to ntotal_points_per_switch
    assert ntotal_points_per_switch == 4 * ramp_points + npoints_per_switch, \
            f"Expect ntotal_pts={ntotal_points_per_switch} with tramp={ramp_points}"
    # plot pre-corrected data
    plot_sample_window(per_switch_rawdata, ramp_points, npoints_per_switch, mrs.nswitches)

    corrected_rawdata = np.zeros_like(per_switch_rawdata)
    
    for i in range(mrs.nswitches):
        # shift = 7 for cirrhrat
        # shift = 5 for ischemia
        shift = int(ramp_points) + 3 # shift first echo by 3 points to the right 
        if i==0:
            corrected_rawdata[i] = np.concatenate((np.zeros(shift), per_switch_rawdata[0,:ntotal_points_per_switch-shift]))
        else:
            corrected_rawdata[i] = per_switch_rawdata.flatten()[i*ntotal_points_per_switch-shift:(i+1)*ntotal_points_per_switch-shift]
    plot_sample_window(corrected_rawdata, ramp_points, npoints_per_switch, mrs.nswitches)

    # calculate the drift based on the largest peaks on each switch
    expected_peak_idx = ramp_points + npoints_per_switch//2  # for cirrhrat data, idx=6 out of 12points are used  
    drift_corrected_rawdata = np.zeros_like(per_switch_rawdata)
    for i in range(mrs.nswitches):
        # every 3 points
        if i<3:
            shift=0
        elif i<6:
            shift=1
        elif i<9:
            shift=2
        # every 5 points
        elif i<14:
            shift=3
        elif i<19:
            shift=4
        elif i<24:
            shift=5
        elif i<29:
            shift=6
        elif i<34:
            shift=7
        elif i<39:
            shift=8
        elif i<44:
            shift=9
        elif i<49:
            shift=10
        elif i<54:
            shift=11
        elif i<59:
            shift=12
        elif i<63:
            shift=13
        elif i == 63:
            print(corrected_rawdata[i,shift:].shape)
            drift_corrected_rawdata[i] = np.concatenate((corrected_rawdata[i,shift:],np.zeros(shift)))
            continue
        drift_corrected_rawdata[i] = np.concatenate((corrected_rawdata[i,shift:], corrected_rawdata[i+1,:shift]))
    plot_sample_window(drift_corrected_rawdata, ramp_points, npoints_per_switch, mrs.nswitches) 
    return True

def plot_sample_window(rawdata: np.ndarray, ramp_points: int, npoints_per_switch: int, nswitches: int) -> bool:
    """
    Plot rawdata sample points across switches and overlay the actual sampling window
    Discard points count is calculated by adding up ramp_points on the side of npoints_per_switch
    ramp_points --- npoints_per_switch --- ramp_points --- ramp_points(rephase) --- ramp_points(rephase)
    Args:
        rawdata: shape in nswitches, ntotal_points_per_switch=(nsamples/nswitches)
        ramp_points: sample points during ramp time: tramp // sample_period
        npoints_per_switch:             
    """
    import matplotlib.pyplot as plt
    figure, axes = plt.subplots(figsize=(6, 9))
    # plot the largest peak as x
    largest_peak = np.argmax(rawdata, axis=1)
    axes.plot(largest_peak, range(nswitches), 'xr', label='largest peak')
    # plot magnitude signal for each switch
    img = axes.imshow(rawdata, vmin=0, origin='lower', aspect='auto')
    # ramp -> npoints_per_switch -> ramp -> rephasing (two ramps)
    axes.axvspan(ramp_points-0.5, ramp_points+npoints_per_switch-0.5, color='C7', alpha=0.5,
                 label=f'sample window of {npoints_per_switch} readout')
    # echo should happen at the middle of npoints_per_switch e.g.) 12points_per_switch->echo at idx7 
    expected_echo_position = ramp_points + (npoints_per_switch // 2)
    axes.axvspan(expected_echo_position-1, expected_echo_position, color='r', alpha=0.3, 
                 label='expected echo position')
    axes.set_xlabel(f"position within the {npoints_per_switch + 4 * ramp_points} point switch")
    axes.set_ylabel("switch")
    axes.legend(fontsize=8, loc='upper right')
    figure.colorbar(img, ax=axes, label='signal aggregated over views and repetitions')
    figure.tight_layout()
    plt.show()
    plt.close(figure)
    return True


# def read_mrd_acq(input_file: BinaryIO, output_file: BinaryIO):
#     """
#     Read a mrd file acquisition field, correct for echo position, and rewrite back a corrected raw file
#     """
#     with mrd.BinaryMrdReader(input_file) as reader:
#         header = reader.read_header()
#         for item in reader.read_data():
#             if not isinstance(item, mrd.StreamItem.Acquisition):
#                 continue
#             # no correction for phantom prescan
#             elif item.value.head.falgs & mrd.AcquisitionFlags.IS_NAVIGATION_DATA:
#                 continue
#             else:



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

    # figure, axes = plt.subplots(3, 1, figsize=(16, 9))
    title = (f"{mrs.sequence_name or 'unknown sequence'}: {report['nswitch']} switches of "
             f"{total} points, {kept} kept, peak/median {report['snr']:.2f}")
    figure.suptitle(title)

    axes[0].plot(positions, profile, 'o-', color='C0')
    # the sequence's own answer to where the window goes, on the same axes as the two measured ones:
    # from the end of the leading ramp through the flat top, which is where the signal was meant to be
    window = sequence_window(mrs)
    if window:
        axes[0].axvspan(window[0] - 0.5, window[1] + 0.5, color='C7', alpha=0.15,
                        label=f"the sequence samples {window[0]}..{window[1]}")
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
    fill_screen(figure)
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


def ramp_samples(mrs: MRSdata) -> Optional[float]:
    """
    How many samples the readout gradient spends on one ramp, or None when the file cannot say.

    tramp is recorded in microseconds and sample_period in units of 100 ns, hence the /10. The result
    is deliberately left fractional: the kidney sequences run a 2.5 sample ramp, and rounding it here
    would push both the echo position and the sampling window half a sample off before anything has
    used them
    """
    if not mrs.tramp or not mrs.sample_period:
        return None
    return mrs.tramp / (mrs.sample_period / 10.0)


def expected_echo_positions(mrs: MRSdata) -> Optional[Tuple[int, int]]:
    """
    Where the sequence says its two echoes sit inside a switch, from its own parameters.

    Not a measurement. The readout loop runs ramp -> points_per_switch -> ramp -> ramp -> ramp and
    repeats, so the switch holds the readout echo at the middle of the flat top and, after the
    gradient rephases, a second crossing of k-space centre

        readout   = tramp / sample_period + points_per_switch / 2
        rephasing = 2 * tramp / sample_period + points_per_switch

    On cirrhrat the loop closes exactly: tramp 112 over a 28 us sample period is a 4 sample ramp, and
    4 ramps of 4 plus the 12 kept points fill a 28 point switch, putting the readout echo at 10 and
    the rephasing one at 20. On the kidney data tramp 100 over 40 us is 2.5 samples, four of which
    overshoot the 8 non-flat samples a 20 point switch has, and the two come out at 8.5 and 17 -
    truncated here, since a roll moves whole samples.

    The first index is the one to trust. Measured against the data the readout echo lands where this
    says it does, but the rephasing echo reads a sample later than 17 on ischemia_179 and the two peak
    families sit 9 to 13.5 positions apart where this predicts 8.5 to 10, so the second index is a
    check to report rather than a position to move anything onto.

    Read only off a sequence the formula was established against, by exact name: a variant that lays
    its ramps out differently would be moved onto a position that is not its own, and being left
    alone is the better failure. That is why ischemia_121_1, which records epsigre43_FB_13C, is
    de-drifted but never moved
    Args:
        - mrs: one parsed or probed MRS file of the epsi family
    Returns:
        - (readout echo, rephasing echo) within a switch, or None when this file carries no
          established formula
    """
    ramp = ramp_samples(mrs)
    if mrs.sequence_name != 'epsigre' or ramp is None:
        return None
    _, _, kept = switch_layout(mrs)
    return int(ramp + kept / 2), int(2 * ramp + kept)


def sequence_window(mrs: MRSdata) -> Optional[Tuple[int, int]]:
    """
    The stretch of a switch the sequence actually samples on the gradient plateau.

    From the first position clear of the leading ramp through the end of the flat top, so
    ceil(tramp / sample_period) .. that + points_per_switch - 1. Rounded up rather than down because
    the sample the ramp is still running through is not a plateau sample: on the kidney layout that
    puts the window at 3..14, whose centre is the 8.5 expected_echo_positions reports, where rounding
    down would centre it on 7.5.

    This is not the window a conversion writes. generate_acquisition centres its discard_pre in the
    switch instead - 4..15 on the kidney layout, 8..19 on cirrhrat - and window_report's pads are all
    anchored to that. The difference between the two is exactly the mrd2recon --pad that would move
    the reconstruction onto the plateau: 1 sample on the kidney data, 4 on cirrhrat, where the centred
    window sits well past the flat top
    Args:
        - mrs: one parsed or probed MRS file of the epsi family
    Returns:
        - (first position, last position) inclusive, or None when the file carries no ramp time
    """
    ramp = ramp_samples(mrs)
    # npoints_per_switch read straight off the file rather than through MRStomrd2.switch_layout,
    # which is the width that layout would report anyway: this module is imported by MRStomrd2, so
    # importing back would be circular, and the window is undrawable without a recorded width in
    # any case
    if mrs.sequence_name != 'epsigre' or ramp is None or not mrs.npoints_per_switch:
        return None
    kept = mrs.npoints_per_switch
    start = int(np.ceil(ramp))
    return start, start + kept - 1


def aligned_echo_position(cube: np.ndarray, view: Optional[int] = None) -> Tuple[int, float]:
    """
    Where the echo actually sits inside a switch, over the switches given.

    Read off the coherent profile spectral_peak builds rather than off a magnitude sum, because this
    position is what a constant move is measured against and so it has to be the sharper of the two
    estimators: on ischemia_179 the coherent profile reaches 13.26 peak over median where the
    magnitude sum manages 1.84, and a position read off the flatter one is a position read off noise.

    The sharpness comes back with it so that a move resting on a smeared profile is visible in the
    report rather than silent
    Args:
        - cube: (switch, position, view, repeats), already rolled straight if it is going to be
        - view: the view to read, or None for the one carrying the most signal
    Returns:
        - (position within the switch, peak over median of the profile it was read from)
    """
    profile, sharpness, _ = spectral_peak(cube, view=view)
    return int(np.argmax(profile)), sharpness


def peak_families(signal: np.ndarray, slope: float, anchor: Optional[int] = None,
                  fallback_separation: Optional[float] = None) -> dict:
    """
    Sort the per-switch peaks into the two echoes a switch carries, against switch 0.

    The brightest position in one switch is whichever of its two echoes happened to win, so a plot of
    those peaks looks like scatter when the echo is in fact perfectly orderly: on ischemia_179 they
    read 3,4,5... and 14,15,16..., two families each stepping one position every ~4.3 switches exactly
    as its +0.232 drift should. Splitting them is what makes that legible, and it is why nothing here
    needs a switch to peak on the readout echo - a switch peaking on the rephasing echo still says
    where the train is.

    Which family is which is decided by switch 0 and nothing else. The signal decays along the train -
    peak magnitude per switch runs 0.70 to 1.00 relative at the front of ischemia_179 against 0.17 to
    0.27 at the back - so the first switch is the most reliable reading there is, and the drift is
    linear from it. Sorting by "whichever mode holds more switches" instead gets the front of the train
    backwards: with the drift left in the residuals it labelled switches 0 to 3 of ischemia_179, all
    peaking at 3 on the readout echo, as the second one.

    The separation between the families is measured rather than taken from the ramp formula, which
    predicts 8.5 to 10 where the data reads 9 to 11 on the kidney scans and about 14 on cirrhrat
    Args:
        - signal: (switch, position) magnitude summed over views and repetitions
        - slope: the drift still present in `signal`, in samples per switch, so that a switch on the
          readout echo has a residual near zero whatever the drift. Zero for a corrected readout
        - anchor: the readout echo's position in switch 0, or None to read it off switch 0's own peak
        - fallback_separation: the spacing the sequence gives its two echoes, used when the data is
          too smeared to measure one; None falls back to half a switch
    Returns:
        - dict of the per-switch peaks, the anchor, the separation and whether it was measured or
          fallen back on, which family each switch fell in and how many that is, and the slope the
          peaks refit to on their own
    """
    nswitch, total = signal.shape
    peaks = np.argmax(signal, axis=1)
    switches = np.arange(nswitch)
    if anchor is None:
        anchor = int(peaks[0])
    # distance from the line the anchor and the slope draw, the short way round the switch
    residual = (peaks - anchor - slope * switches + total / 2) % total - total / 2

    # where the switches that did not peak on the readout echo cluster. Taken as the most populated
    # whole position rather than as a mean, since the set holds noise as well as the second echo and
    # one stray residual drags a mean a long way - on ischemia_179 the mean reads -5.2 where the
    # cluster sits at +9. A quarter of a switch either side of zero is left out because the readout
    # family itself spreads that far: on ischemia_179 its residuals walk from 0 to -3 across the
    # train, and at an exclusion of 3 that tail held 9 switches and won the mode outright, reading a
    # separation of -3 where the second echo plainly sits at +9
    census = np.zeros(total, dtype=int)
    for value in residual:
        census[int(round(value)) % total] += 1
    radius = max(3, total // 4)
    apart = [position for position in range(total)
             if min(position % total, (-position) % total) >= radius]
    mode = max(apart, key=lambda position: census[position]) if apart else 0
    # counted over the mode and its neighbours, since the family spreads across two or three whole
    # positions as the drift walks it: on ischemia_179 its 12 switches sit at +8, +9 and +10, and the
    # strongest single position holds only 5 of them
    population = int(sum(census[(mode + offset) % total] for offset in (-1, 0, 1)))
    # believed only when enough switches sit there. On cirrhrat_43_1 the census is flat, because that
    # readout is too smeared for a per switch argmax to mean anything, and a separation read off noise
    # would then decide which echo switch 0 peaked on. The sequence's own spacing is the fallback, and
    # the report says which was used
    trusted = population >= max(6, nswitch // 8)
    if trusted:
        separation = float((mode + total / 2) % total - total / 2)
    else:
        separation = float(fallback_separation if fallback_separation else total / 2)

    # each switch to whichever line is nearer, so every one lands in a family
    to_readout = np.abs(residual)
    to_second = np.abs((residual - separation + total / 2) % total - total / 2)
    readout_family = to_readout <= to_second
    on_readout = int(readout_family.sum())

    # the slope the peaks alone say, with the or-condition: a switch counts when it sits within a
    # sample of either line, so switches that peaked on the rephasing echo constrain it too
    candidates = np.arange(-DRIFT_SLOPE_LIMIT, DRIFT_SLOPE_LIMIT + DRIFT_FINE_STEP / 2,
                           DRIFT_FINE_STEP)
    inliers = []
    for candidate in candidates:
        offset = (peaks - anchor - candidate * switches + total / 2) % total - total / 2
        near_readout = np.abs(offset) <= 1
        near_second = np.abs((offset - separation + total / 2) % total - total / 2) <= 1
        inliers.append(int((near_readout | near_second).sum()))
    best = int(np.argmax(inliers))
    return dict(peaks=peaks, anchor=int(anchor), separation=separation, measured=bool(trusted),
                population=population, readout_family=readout_family, on_readout=on_readout,
                on_second=nswitch - on_readout,
                fitted_slope=float(candidates[best]), inliers=int(inliers[best]))


def rephasing_peak(profile: np.ndarray, readout: int, gap: int) -> Tuple[int, float]:
    """
    The strongest position of the readout profile that is not the readout echo.

    Searched outside readout +- gap so the answer cannot be the flank of the echo it is being compared
    against. The amplitude comes back as a fraction of the readout peak rather than being tested
    against a threshold, because that fraction is what tells a real second echo from a plateau: it
    reads 0.36 on ischemia_179 and 0.52 on ischemia_121_1, both a clear peak half a switch away from
    the readout one, against 0.83 on cirrhrat_43_1 - which is not a second echo at all but the
    shoulder of a profile too smeared to have two of anything
    Args:
        - profile: signal per position within the switch, as spectral_peak returns
        - readout: the position of the readout echo
        - gap: how far either side of it to exclude
    Returns:
        - (position, its amplitude as a fraction of the profile's peak)
    """
    total = len(profile)
    searchable = np.ones(total, dtype=bool)
    for offset in range(-gap, gap + 1):
        searchable[(readout + offset) % total] = False
    if not searchable.any():
        return readout, 1.0
    positions = np.arange(total)[searchable]
    found = int(positions[np.argmax(profile[searchable])])
    return found, float(profile[found] / profile.max()) if profile.max() else 0.0


def measure_group_drift(mrs_list: Sequence[MRSdata]
                       ) -> Optional[Tuple[dict, Tuple[int, int, int], np.ndarray]]:
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
    The pooled cube comes back with the measurement rather than being built and dropped, because
    every caller reads the echo position out of it straight afterwards and pooling 25 repetitions of
    12 views twice to answer two questions about the same samples is work for nothing
    Args:
        - mrs_list: the parsed files of one dataset, read rather than probed
    Returns:
        - (the measurement from measure_echo_drift, (nswitch, total, kept), the pooled cube), or None
          when the list held no EPSI readout to pool
    """
    epsi = [mrs for mrs in mrs_list if is_epsi(mrs) and mrs.rawdata is not None]
    if not epsi:
        return None
    cube = pooled_switch_cube(epsi)
    if cube is None:
        return None
    nswitch, total, kept = switch_layout(epsi[0])
    return measure_echo_drift(cube, nswitch, total), (nswitch, total, kept), cube


def echo_alignment(mrs_list: Sequence[MRSdata], slope: Optional[float] = None) -> Optional[dict]:
    """
    Everything needed to put one dataset's echo where it belongs: the drift to take out, and the
    constant move that lands the result on the position the sequence asks for.

    Two corrections, worked out together and returned as one shift per switch, so applying them is a
    single roll and neither can be applied without the other:

      - the drift, which is the echo walking along the switch train. Removed by switch_shifts, whose
        middle-of-the-train anchor is left exactly as it is: the anchor decides nothing here, since
        the constant below moves the aligned echo onto its target whatever the anchor chose, and
        measure_echo_drift's accept gates are calibrated against the scores that anchor produces.
      - the constant, which is the whole readout sitting at the wrong position within the switch.
        expected_echo_position says where the echo belongs, aligned_echo_position measures where it
        is, and the difference is the move.

    The two are independent, so the constant is worked out whether or not the drift was usable: a
    scan whose echo does not walk still has it somewhere other than where the sequence puts it. What
    says whether that move can be trusted is the profile sharpness reported beside it.

    On the two datasets this was established against the constant reproduces, from the sequence
    parameters alone, the per-tramp prepend the conversion path hard-codes: switch 0 comes out moved
    by +5 on the tramp 100 kidney data and +7 on tramp 112 cirrhrat
    Args:
        - mrs_list: the parsed files of one dataset, read rather than probed
        - slope: roll rate in samples per switch to apply instead of the measured one, for correcting
          a scan case by case. Zero is a real answer, meaning take no drift out but still make the
          constant move; None leaves the measurement in charge
    Returns:
        - dict of the drift, the layout, the expected and measured positions, the constant, the total
          shift per switch, where the echo ends up, the two peak families, the rephasing check and
          the sequence's own sampling window, or None when there was nothing to pool
    """
    measured = measure_group_drift(mrs_list)
    if measured is None:
        return None
    drift, layout, cube = measured
    nswitch, total, kept = layout
    reference = next(mrs for mrs in mrs_list if is_epsi(mrs) and mrs.rawdata is not None)

    # a slope given by hand wins over the measured one, and the measurement still runs and is still
    # reported: the operator correcting a scan case by case wants to see what the data said next to
    # what they asked for. Given zero is a real answer and not a missing one, hence the `is not None`
    given = slope is not None
    # zero, not the measured slope, when the measurement refused it: what is applied has to be what
    # the shifts below actually carry, or every number reported against it describes something else
    applied = float(slope) if given else (drift['slope'] if drift['usable'] else 0.0)
    # zeros rather than a skip when there is no drift to take out, so that the constant below is
    # applied to a readout that was left alone as readily as to one that was straightened
    shifts = (switch_shifts(nswitch, total, applied) if applied
              else np.zeros(nswitch, dtype=int))
    straightened = roll_switches(cube, shifts)
    aligned, sharpness = aligned_echo_position(straightened)

    # the two echoes as the data shows them, read in the de-drifted frame the constant is added in.
    # The sequence's own spacing goes in as the fallback for a readout too smeared to measure one
    positions = expected_echo_positions(reference)
    expected, expected_rephasing = positions if positions else (None, None)
    by_sequence = (expected_rephasing - expected) if positions else None
    signal = np.abs(straightened).sum(axis=(2, 3))
    families = peak_families(signal, 0.0, fallback_separation=by_sequence)
    separation = families['separation']

    # switch 0 is the anchor: the signal decays along the train, so the first switch is the most
    # reliable reading of where the echo is, and the drift is linear from it. But switch 0 peaks on
    # whichever of its two echoes was brighter, and on cirrhrat_43_1 that is the rephasing one - its
    # peak reads 17 where its neighbours read 3, 4, 5 - so which echo it landed on is decided against
    # the pooled profile before the anchor is believed, and the separation backed out when it was the
    # second. Anchoring on it blindly would move that scan's readout echo to position 1 of 28
    anchor = int(families['anchor'])
    to_readout = abs((anchor - aligned + total / 2) % total - total / 2)
    to_second = abs((anchor - aligned - separation + total / 2) % total - total / 2)
    anchor_on_readout = to_readout <= to_second
    if not anchor_on_readout:
        anchor = int(round(anchor - separation)) % total

    # measured from switch 0 rather than from the pooled profile, so the switch carrying the most
    # signal is the one that lands exactly where the sequence says. Reduced the short way round the
    # switch, since a position is cyclic within one: moving an echo from 18 to 2 of 20 is 4 samples
    # later, not 16 earlier
    constant = int((expected - anchor + total // 2) % total - total // 2) if expected is not None else 0
    shifts = shifts + constant

    # the same families read off the readout as acquired, with the applied slope taken out of the
    # residuals rather than out of the samples: that leaves the slope this refits an estimate from the
    # peaks alone rather than one the correction has already been baked into
    raw_families = peak_families(np.abs(cube).sum(axis=(2, 3)), applied,
                                 fallback_separation=by_sequence)
    corrected = roll_switches(cube, shifts) if constant else straightened
    rephasing, rephasing_amplitude = (rephasing_peak(spectral_peak(corrected)[0],
                                                     (anchor + constant) % total,
                                                     max(int(round(abs(separation) / 2)), 1))
                                      if expected is not None else (None, 0.0))

    discard_pre = (total - kept) // 2
    landed = (anchor + constant) % total
    # judged against the window the sequence samples where there is one, since that is where the echo
    # was meant to land; a conversion's centred window is the fallback and the --pad above is the gap
    window = sequence_window(reference)
    window_start = window[0] if window else discard_pre
    return dict(drift=drift, layout=layout, shifts=shifts, cube=cube,
                applied=applied, given=given,
                expected=expected, expected_rephasing=expected_rephasing,
                aligned=aligned, sharpness=sharpness, constant=constant,
                anchor=anchor, anchor_on_readout=bool(anchor_on_readout), separation=separation,
                separation_measured=bool(families['measured']),
                rephasing=rephasing, rephasing_amplitude=rephasing_amplitude,
                families=raw_families, window=window, window_start=window_start,
                landed=landed, discard_pre=discard_pre,
                inside=bool(((landed - window_start) % total) < kept),
                sequence_name=reference.sequence_name, tramp=reference.tramp,
                sample_period=reference.sample_period)


def shift_report(alignment: dict, label: str = "") -> None:
    """
    Print what a dataset's echo does, and what correcting it comes to.

    The shift range is printed against the discarded ramp on purpose: rolling a switch brings in the
    samples of its neighbour, which are ramp points rather than signal, so once the shift exceeds
    the ramp the ends of the train have nothing valid left to move in. That is the number that says
    whether a measured drift can be corrected by moving whole samples at all
    Args:
        - alignment: the dict echo_alignment returns
        - label: what to call this dataset, e.g. a meas_id
    Returns:
        - None; everything goes to stdout beside the figures
    """
    drift, (nswitch, total, kept) = alignment['drift'], alignment['layout']
    shifts, expected = alignment['shifts'], alignment['expected']
    families = alignment['families']

    print(f"\necho of {label or 'the pooled readout'}, {drift['reps']} repetitions pooled")
    print(f"  {'drift' if drift['usable'] else 'no usable drift'}: {drift['reason']}")
    print(f"  best slope {drift['slope']:+.4f} samples per switch, a switch period of "
          f"{drift['period']:.2f} samples where nsamples/{nswitch} records {total}")
    if alignment['given']:
        print(f"  applying {alignment['applied']:+.4f} per switch as asked for, where the search "
              f"measured {drift['slope']:+.4f}")
    # the same train read a second way: peaks per switch, sorted into the two echoes a switch holds
    how = ("measured" if alignment['separation_measured']
           else "too smeared to measure, taken from the sequence")
    print(f"  per switch peaks: {families['on_readout']} of {nswitch} on the readout echo and "
          f"{families['on_second']} on the second, {alignment['separation']:+.0f} positions apart "
          f"({how}), which refits a slope of {families['fitted_slope']:+.4f} "
          f"({families['inliers']} of {nswitch} switches within a sample) against the "
          f"{alignment['applied']:+.4f} applied")
    # what the constant below is measured from: switch 0, which carries the most signal of any switch
    landed_on = 'the readout echo' if alignment['anchor_on_readout'] else 'the rephasing echo'
    print(f"  switch 0 peaks on {landed_on}, so the readout echo of switch 0 sits at "
          f"{alignment['anchor']} once the drift is out"
          + ('' if alignment['anchor_on_readout']
             else f", backed out by the {alignment['separation']:.1f} between the two"))
    if alignment['window']:
        start, end = alignment['window']
        pad = alignment['discard_pre'] - start
        print(f"  the sequence samples {start}..{end} of the switch, where a conversion keeps "
              f"{alignment['discard_pre']}..{alignment['discard_pre'] + kept - 1}: mrd2recon --pad "
              f"{pad} reads the plateau rather than a window {pad} sample(s) past it")
    if expected is None:
        print(f"  no expected position for sequence '{alignment['sequence_name']}', so the echo is "
              f"left where it was acquired: the ramp layout is only established for epsigre")
    else:
        ramp_us = alignment['sample_period'] / 10.0
        print(f"  the sequence puts the echo at {expected} = tramp {alignment['tramp']}us / "
              f"{ramp_us:.0f}us sample period + {kept}/2, against the {alignment['anchor']} switch 0 "
              f"reads (the pooled profile agrees at {alignment['aligned']}, peak/median "
              f"{alignment['sharpness']:.2f})")
        start = alignment['window_start']
        print(f"  so the whole readout moves {alignment['constant']:+d}, landing the echo at "
              f"{alignment['landed']}, {'inside' if alignment['inside'] else 'OUTSIDE'} the "
              f"{start}..{start + kept - 1} the sequence samples")
        # the second echo checks the first rather than moving anything: the readout echo is what the
        # constant above was measured from, and this says whether the switch looks like the sequence
        # says it should once that move is made
        off = alignment['rephasing'] - alignment['expected_rephasing']
        print(f"  the rephasing echo lands at {alignment['rephasing']} against the "
              f"{alignment['expected_rephasing']} the sequence puts it at, {off:+d} off, at "
              f"{alignment['rephasing_amplitude']:.2f} of the readout peak")
        if abs(off) > 1:
            print(f"WARNING the rephasing echo of {label or 'this readout'} sits {off:+d} samples "
                  f"from where the sequence puts it, so the ramp layout this correction was worked "
                  f"out from does not describe this scan; the move was applied anyway",
                  file=sys.stderr)
    discard = (total - kept) // 2
    print(f"  correcting it moves switch 0 by {int(shifts[0]):+d} and switch {nswitch - 1} by "
          f"{int(shifts[-1]):+d}, a range of {int(shifts.min())} to {int(shifts.max())} samples "
          f"against {discard} discarded ramp points either side")


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
    alignment = echo_alignment([mrs for _, mrs in data]) if data else None
    if alignment is None:
        print("\nno unaveraged EPSI readout to measure the echo drift from", file=sys.stderr)
    else:
        shift_report(alignment, label=f"{len(data)} file(s)")

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
    mode.add_argument("-i", "--input", type=Path,
                      help="single MRS .MRD file")
    args = parser.parse_args()

    if args.input and not args.input.is_file():
        parser.error(f"{args.input} is not a file")
    
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
