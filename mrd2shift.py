"""
Check a converted MRD v2 stream, and optionally roll its EPSI echoes straight.

MRStomrd2 -w asks where the echo sits from the raw .MRD side. This asks the same question of the
converted file, from the reconstruction's point of view: does the stream carry the switch layout,
the discard counts, the dwell time and the ramp time mrd2recon needs, do the acquisitions present
agree with the encoding limits the header promises, and where does the echo actually peak in each
switch as written.

Reporting is the default and the point. mrd2recon reads every switch of a readout at one fixed
offset, `iswitch * totalppswitch + discard_pre - pad`, so an echo that walks along the switch train
leaves that window and nothing downstream says so. With --output the same stream is written back
with each switch rolled by its own shift, which is the one arrangement that can take a drift out:
an acquisition carries every switch of the readout but only one discard_pre, so a header field
cannot express a per-switch correction.

    python mrd2shift.py -i raw.mrd2                      # report only
    python mrd2shift.py -i raw.mrd2 -o straight.mrd2     # report, and write the rolled stream

Both default to the Tyger buffer pipes ($INPUT_PIPE, $OUTPUT_PIPE), so a codespec needs no
arguments. The whole stream is held in memory: the drift is measured from every acquisition before
the first one can be written, and a buffer FIFO can only be read once.

What is measured, and whether it is worth acting on, is epsi_window's decision - measure_echo_drift
searches the correction and scores it, and refuses a slope that only lines noise up. Nothing here
second-guesses that: an unusable measurement copies the stream through untouched.
"""

from __future__ import annotations

import argparse
import math
import os
import sys
from collections import defaultdict
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

import mrd
from epsi_window import measure_echo_drift, switch_shifts

# What this program took out of a stream, recorded in the output header so a second pass can refuse
# rather than roll an already straightened readout twice
SLOPE_PARAMETER = "echo_drift_slope_applied"


# ---------- reading what the converter wrote -------------------------------------------------


def is_epsi_acquisition(acq: mrd.Acquisition) -> bool:
    """
    Whether an acquisition records an EPSI switch layout.

    The twin of mrd2recon.is_epsi_acquisition, re-implemented rather than imported: mrd2recon
    references mrd.ArrayType, which the pinned mrd fork revision does not carry, so importing it
    fails outright. Keep the two readings in step
    """
    user_int = list(acq.head.user_int or [])
    return len(user_int) >= 2 and int(user_int[0]) > 0


def switch_layout(acq: mrd.Acquisition) -> Tuple[int, int, int]:
    """
    How one EPSI readout divides into gradient switches, read off the acquisition.

    MRStomrd2.generate_acquisition writes [nswitch, points_per_switch] into user_int, because
    idx.contrast carries the sequence's own echo index rather than the switches packed into one
    readout. The kept width is what is left once the discard points come off both ends. The twin of
    mrd2recon.switch_layout
    Args:
        - acq: one acquisition of the readout
    Returns:
        - (nswitch, total points in one switch, points kept per switch)
    Raises:
        - ValueError on an acquisition that records no layout
    """
    user_int = list(acq.head.user_int or [])
    if len(user_int) < 2 or int(user_int[0]) <= 0:
        raise ValueError("this acquisition records no switch layout in user_int")
    nswitch, total = int(user_int[0]), int(user_int[1])
    return nswitch, total, total - (acq.head.discard_pre or 0) - (acq.head.discard_post or 0)


def header_user_long(header: mrd.Header, name: str) -> Optional[int]:
    """The named integer user parameter of a header, or None when it carries no such record"""
    user = getattr(header, "user_parameters", None)
    for param in (getattr(user, "user_parameter_long", None) or []) if user else []:
        if param.name == name:
            return int(param.value)
    return None


def header_user_double(header: mrd.Header, name: str) -> Optional[float]:
    """The named float user parameter of a header, or None when it carries no such record"""
    user = getattr(header, "user_parameters", None)
    for param in (getattr(user, "user_parameter_double", None) or []) if user else []:
        if param.name == name:
            return float(param.value)
    return None


def recon_pad(header: mrd.Header, acq: mrd.Acquisition) -> Tuple[Optional[int], str]:
    """
    The pad mrd2recon will read this file with, and where it comes from.

    The ramp at the head of a switch is not yet on the gradient plateau, so the usable window opens
    ceil(tramp / dwell) samples in, while discard_pre assumes the discarded points split evenly
    across both ends; the pad is the difference. This is mrd2recon.epsi_leading_pad's arithmetic,
    except that where that function falls back to EPSIGRE_DEFAULT_PAD - a name it uses twice and
    defines nowhere, so the fallback raises NameError - this reports that there is no usable pad,
    which is the finding about the file
    Args:
        - header: the stream header, read for its tramp record
        - acq: any acquisition of the readout, read for its dwell and discard count
    Returns:
        - (pad, a one line account of it), the pad being None when it cannot be derived
    """
    tramp_us = header_user_long(header, "tramp")
    if tramp_us is None:
        return None, ("no pad: the header records no ramp time, and mrd2recon's fallback for that "
                      "raises NameError rather than defaulting; convert again to record tramp")
    if not acq.head.sample_time_ns:
        return None, ("no pad: the acquisition records no dwell time, and mrd2recon's fallback for "
                      "that raises NameError rather than defaulting")
    # ceil, because a window opening part way through the last ramp sample still includes it
    ramp = int(math.ceil(tramp_us * 1000.0 / acq.head.sample_time_ns))
    discard_pre = acq.head.discard_pre or 0
    return discard_pre - ramp, (f"pad {discard_pre - ramp} from tramp={tramp_us}us over a "
                                f"{acq.head.sample_time_ns}ns dwell, i.e. a {ramp} sample ramp "
                                f"against discard_pre={discard_pre}")


def read_stream(path: str) -> Tuple[mrd.Header, List[mrd.StreamItem]]:
    """
    Read a whole MRD v2 stream into memory.

    Materialised rather than streamed for two reasons: the drift is measured from every acquisition
    before the first one can be written, and an input that is a Tyger buffer FIFO can only be read
    once. Non-acquisition items are kept in order so a rewrite passes them through untouched
    Args:
        - path: an .mrd2 file, or a FIFO carrying one
    Returns:
        - (header, every stream item in order)
    """
    # opening a FIFO for read blocks until the buffer sidecar opens the write end
    with open(path, "rb") as stream:
        with mrd.BinaryMrdReader(stream) as reader:
            header = reader.read_header()
            items = list(reader.read_data())
    return header, items


def acquisitions_of(items: Sequence[mrd.StreamItem]) -> List[mrd.Acquisition]:
    """The acquisitions of a stream, in order. They are the items this program reads and rewrites"""
    return [item.value for item in items if isinstance(item, mrd.StreamItem.Acquisition)]


def is_navigation(acq: mrd.Acquisition) -> bool:
    """
    Whether an acquisition is one of the averaged prescans MRStomrd2 folds into the same stream,
    which it flags rather than writing to a file of its own
    """
    return bool(int(acq.head.flags) & int(mrd.AcquisitionFlags.IS_NAVIGATION_DATA))


# ---------- what the file says about itself --------------------------------------------------


def report_header(header: mrd.Header) -> None:
    """Print the header fields a reconstruction reads, and say which are missing"""
    meas = header.measurement_information
    print("header")
    print(f"  measurement    {getattr(meas, 'measurement_id', None) or '(none)'}")
    print(f"  sequence       {getattr(meas, 'sequence_name', None) or '(none)'}")
    print(f"  frequency      {header.experimental_conditions.h1resonance_frequency_hz} Hz")
    tramp = header_user_long(header, "tramp")
    print(f"  tramp          {f'{tramp} us' if tramp is not None else '(none) - see the pad below'}")
    applied = header_user_double(header, SLOPE_PARAMETER)
    if applied is not None:
        print(f"  {SLOPE_PARAMETER} {applied:+.4f} samples per switch already taken out of this "
              f"stream")
    if not header.encoding:
        print("  encoding       (none) - nothing records the matrix this was acquired at")
        return
    space = header.encoding[0].encoded_space
    print(f"  matrix         {space.matrix_size.x} x {space.matrix_size.y} x {space.matrix_size.z}")


def limit_maximum(header: mrd.Header, name: str) -> Optional[int]:
    """The maximum index the header's encoding limits promise for one axis, if it records that axis"""
    if not header.encoding:
        return None
    limit = getattr(header.encoding[0].encoding_limits, name, None)
    return None if limit is None else int(limit.maximum)


def report_census(header: mrd.Header, items: Sequence[mrd.StreamItem],
                  acqs: Sequence[mrd.Acquisition]) -> int:
    """
    Print what the stream actually holds, checked against what its header promises.

    An index the acquisitions use beyond the limit the header records is the conversion bug this
    program exists to catch: a reader sizes its k-space off those limits, so an acquisition outside
    them is data the reconstruction has nowhere to put
    Args:
        - header, items, acqs: the stream, its items and just its acquisitions
    Returns:
        - how many axes disagree with the header
    """
    others = len(items) - len(acqs)
    epsi = sum(1 for acq in acqs if is_epsi_acquisition(acq))
    navigation = sum(1 for acq in acqs if is_navigation(acq))
    print(f"\nacquisitions: {len(acqs)}"
          f"{f', plus {others} other stream item(s)' if others else ''}")
    print(f"  {epsi} record an EPSI switch layout, {len(acqs) - epsi} do not")
    print(f"  {navigation} flagged IS_NAVIGATION_DATA, i.e. the prescans folded into this stream")

    disagreements = 0
    for label, index, limit_name in (
            ("kspace_encode_step_1", lambda i: i.kspace_encode_step_1, "kspace_encoding_step_1"),
            ("kspace_encode_step_2", lambda i: i.kspace_encode_step_2, "kspace_encoding_step_2"),
            ("slice", lambda i: i.slice, "slice"),
            ("contrast", lambda i: i.contrast, "contrast"),
            ("repetition", lambda i: i.repetition, "repetition")):
        used = sorted({int(index(acq.head.idx) or 0) for acq in acqs})
        if not used:
            continue
        promised = limit_maximum(header, limit_name)
        note = ""
        if promised is None:
            note, disagreements = "  <- the header records no limit for this axis", disagreements + 1
        elif used[-1] > promised:
            note = f"  <- beyond the header's maximum of {promised}"
            disagreements += 1
        print(f"  {label:<20} {len(used)} value(s), 0..{used[-1]}, header maximum "
              f"{'(none)' if promised is None else promised}{note}")
    return disagreements


def layout_key(acq: mrd.Acquisition) -> Tuple[int, int, int, int]:
    """
    What makes two acquisitions share a readout geometry: the switch layout and the discards.

    A stream holds the data and the prescans that calibrate it, acquired at different matrices, so
    the shifts one needs are not the shifts the other needs even though the drift is the same
    """
    nswitch, total, _ = switch_layout(acq)
    return nswitch, total, int(acq.head.discard_pre or 0), int(acq.head.discard_post or 0)


def group_layouts(acqs: Sequence[mrd.Acquisition]) -> Dict[Tuple[int, int, int, int],
                                                           List[mrd.Acquisition]]:
    """Every EPSI acquisition of a stream, bucketed by readout geometry, biggest bucket first"""
    layouts: Dict[Tuple[int, int, int, int], List[mrd.Acquisition]] = defaultdict(list)
    for acq in acqs:
        if is_epsi_acquisition(acq):
            layouts[layout_key(acq)].append(acq)
    return dict(sorted(layouts.items(), key=lambda item: len(item[1]), reverse=True))


def report_layout(header: mrd.Header, key: Tuple[int, int, int, int],
                  acqs: Sequence[mrd.Acquisition]) -> None:
    """Print one readout geometry and the window a reconstruction would read out of it"""
    nswitch, total, discard_pre, discard_post = key
    kept = total - discard_pre - discard_post
    navigation = sum(1 for acq in acqs if is_navigation(acq))
    print(f"\nlayout {nswitch} switches of {total} points, keeping {kept}: {len(acqs)} "
          f"acquisition(s), {navigation} of them prescan")
    print(f"  user_int=[{nswitch}, {total}], discard_pre={discard_pre}, "
          f"discard_post={discard_post}, dwell={acqs[0].head.sample_time_ns}ns")
    pad, account = recon_pad(header, acqs[0])
    if pad is None:
        print(f"  {account}")
    else:
        start = discard_pre - pad
        print(f"  {account}")
        print(f"  mrd2recon reads positions {start}..{start + kept - 1} of every switch")


# ---------- where the echo actually is ------------------------------------------------------


def acquisition_cube(acqs: Sequence[mrd.Acquisition], nswitch: int,
                     total: int) -> Tuple[np.ndarray, int]:
    """
    One layout's acquisitions laid out as (switch, position, view, repetition).

    The shape epsi_window measures from, rebuilt from the stream rather than from MRSdata: each
    acquisition is one view of one repetition, and its samples split switch-major exactly as they
    were acquired. Samples past the last whole switch are dropped, the same truncation conversion
    applies
    Args:
        - acqs: the acquisitions of one readout geometry
        - nswitch, total: that geometry, from switch_layout
    Returns:
        - (cube, how many grid points no acquisition filled)
    """
    views = sorted({int(acq.head.idx.kspace_encode_step_1 or 0) for acq in acqs})
    reps = sorted({int(acq.head.idx.repetition or 0) for acq in acqs})
    view_at = {view: i for i, view in enumerate(views)}
    rep_at = {rep: i for i, rep in enumerate(reps)}

    cube = np.zeros((nswitch, total, len(views), len(reps)), dtype=np.complex64)
    filled = np.zeros((len(views), len(reps)), dtype=bool)
    for acq in acqs:
        samples = np.asarray(acq.data)
        # MRS acquires on one channel; summed rather than assumed, so a multi coil stream still reads
        line = samples[0] if samples.shape[0] == 1 else samples.sum(axis=0)
        iview = view_at[int(acq.head.idx.kspace_encode_step_1 or 0)]
        irep = rep_at[int(acq.head.idx.repetition or 0)]
        cube[:, :, iview, irep] = line[:nswitch * total].reshape(nswitch, total)
        filled[iview, irep] = True
    return cube, int((~filled).sum())


def echo_positions(cube: np.ndarray, total: int, kept: int, discard_pre: int) -> dict:
    """
    Where the pooled signal peaks in each switch, and how much of it the kept window catches.

    Magnitude summed over views and repetitions, which is the same measure
    MRStomrd2.shift_echo_position uses, so the two sides of a conversion can be compared directly
    Args:
        - cube: (switch, position, view, repetition), from acquisition_cube
        - total, kept, discard_pre: the readout geometry
    Returns:
        - dict of the pooled signal, the peak per switch and how many switches peak inside the window
    """
    signal = np.abs(cube).sum(axis=(2, 3))
    peaks = np.argmax(signal, axis=1)
    # counted cyclically, the way generate_acquisition addresses the window, so a window running off
    # the end of a switch is not mistaken for one that peaks outside it
    inside = int((((peaks - discard_pre) % total) < kept).sum())
    # the profile of the whole readout, and how sharp it is. This is what a drift smears and what
    # rolling it out restores, and the number moves when the per switch peak count does not: an
    # argmax can land on a second feature in a few switches whether the echo drifts or not
    profile = signal.sum(axis=0)
    median = np.median(profile)
    return dict(signal=signal, peaks=peaks, inside=inside, profile=profile,
                profile_peak=int(np.argmax(profile)),
                sharpness=float(profile.max() / median) if median else np.inf)


def report_echo(geometry: dict, drift: Optional[dict], nswitch: int, total: int, kept: int) -> None:
    """Print where the echoes are and what removing their drift would take"""
    peaks = geometry['peaks']
    print(f"  the peaks visit positions {int(peaks.min())} to {int(peaks.max())}, at "
          f"{int(peaks[0])} in the first switch and {int(peaks[-1])} in the last")
    print(f"  {geometry['inside']} of {nswitch} switches peak inside the kept window")
    print(f"  pooled over every switch the profile peaks at {geometry['profile_peak']}, "
          f"peak/median {geometry['sharpness']:.2f}")
    if drift is None:
        return
    print(f"  {'drift' if drift['usable'] else 'no usable drift'}: {drift['reason']}")
    print(f"  best slope {drift['slope']:+.4f} samples per switch, a switch period of "
          f"{drift['period']:.2f} samples where the layout records {total}")
    if drift['usable']:
        shifts = switch_shifts(nswitch, total, drift['slope'])
        discard = (total - kept) // 2
        print(f"  taking it out moves each switch by {int(shifts.min())} to {int(shifts.max())} "
              f"samples, against {discard} discarded ramp points")


def plot_echo(geometry: dict, key: Tuple[int, int, int, int], name: str) -> None:
    """
    Draw the switch by position map with the peaks marked, the same picture MRStomrd2 -w draws from
    the raw side. Imported here rather than at module scope so the Tyger path needs no plotting stack
    """
    import matplotlib.pyplot as plt

    nswitch, total, discard_pre, discard_post = key
    kept = total - discard_pre - discard_post
    signal, peaks = geometry['signal'], geometry['peaks']
    figure, axes = plt.subplots(figsize=(9, 7))
    axes.imshow(signal, aspect='auto', origin='lower', interpolation='nearest',
                extent=(-0.5, total - 0.5, -0.5, nswitch - 0.5))
    axes.axvspan(discard_pre - 0.5, discard_pre + kept - 0.5, color='w', alpha=0.15,
                 label=f"kept window {discard_pre}..{discard_pre + kept - 1}")
    axes.plot(peaks, np.arange(nswitch), 'x', color='C3', ms=6, label='peak of each switch')
    axes.set_xlabel(f"position within the {total} point switch")
    axes.set_ylabel("switch")
    axes.set_title(f"{name}: echo position per switch, read off the converted stream")
    axes.legend(fontsize=8, loc='upper right')
    figure.colorbar(axes.images[0], ax=axes, label='signal summed over views and repetitions')
    figure.tight_layout()
    plt.show()


# ---------- rolling the echoes straight -----------------------------------------------------


def roll_acquisition(acq: mrd.Acquisition, shifts: np.ndarray, nswitch: int, total: int) -> None:
    """
    Move each switch of one readout by its own shift, in place.

    The only representation a drift fits in: one acquisition carries every switch of the readout but
    a single discard_pre, so no header field can say "later in this switch than in that one".
    Samples past the last whole switch are left alone, and no header field is touched - switch_shifts
    is measured against the middle of the train, so the mean echo position, and with it discard_pre
    and every --pad already established, go on meaning what they mean today
    Args:
        - acq: the acquisition to rewrite
        - shifts: integer shift per switch, from epsi_window.switch_shifts
        - nswitch, total: the readout geometry this acquisition was grouped on
    Returns:
        - None; acq.data, and acq.phase where it carries anything, are replaced
    """
    data = np.array(acq.data, copy=True)                # (coils, nsamples)
    used = nswitch * total
    # reshaped through an explicit copy and written back, since a slice of a multi coil acquisition
    # is not contiguous and an in place roll on a reshape of it would be lost
    body = data[:, :used].reshape(data.shape[0], nswitch, total).copy()
    for iswitch, shift in enumerate(shifts):
        if shift:
            body[:, iswitch, :] = np.roll(body[:, iswitch, :], int(shift), axis=1)
    data[:, :used] = body.reshape(data.shape[0], used)
    acq.data = data

    phase = np.asarray(acq.phase) if acq.phase is not None else None
    # conversion writes zeros here, which roll to zeros, so this only runs on a stream that carries
    # a real phase - and then it has to move with the samples it describes
    if phase is not None and phase.size == data.shape[1] and np.any(phase):
        rolled = np.array(phase, copy=True)
        block = rolled[:used].reshape(nswitch, total).copy()
        for iswitch, shift in enumerate(shifts):
            if shift:
                block[iswitch] = np.roll(block[iswitch], int(shift))
        rolled[:used] = block.reshape(used)
        acq.phase = rolled


def record_slope(header: mrd.Header, slope: float) -> None:
    """
    Note in the header what was taken out, so a second pass refuses rather than rolling twice.

    Accumulated rather than overwritten: a stream rolled twice with --force has had the sum of both
    slopes taken out of it, and that is what the record should say
    """
    if header.user_parameters is None:
        header.user_parameters = mrd.UserParametersType()
    for param in header.user_parameters.user_parameter_double:
        if param.name == SLOPE_PARAMETER:
            param.value = float(param.value) + slope
            return
    header.user_parameters.user_parameter_double.append(
        mrd.UserParameterDoubleType(name=SLOPE_PARAMETER, value=float(slope)))


def write_stream(path: str, header: mrd.Header, items: Sequence[mrd.StreamItem]) -> None:
    """
    Write the header and every item back out as one MRD v2 stream.

    Opened even when nothing was changed: on Tyger this is the output buffer's FIFO, and a stream
    that never opens leaves the sidecar blocked rather than telling it the job is done
    """
    with open(path, "wb") as output:
        # the writer must be closed to emit the end-of-stream sentinel, hence the with block
        with mrd.BinaryMrdWriter(output) as writer:
            writer.write_header(header)
            writer.write_data(items)


def shift_stream(header: mrd.Header, acqs: Sequence[mrd.Acquisition], slope: float) -> int:
    """
    Roll every EPSI acquisition of a stream so its echoes stop moving along the switch train.

    The drift is a property of the gradient timing, so the one slope applies to the prescans as well
    as to the data; the shifts are recomputed per geometry, since the same slope lands on different
    integers in a 28 point switch than in a 34 point one
    Args:
        - header: the stream header, which records the slope taken out
        - acqs: every acquisition of the stream
        - slope: drift to remove, in samples per switch
    Returns:
        - how many acquisitions were rolled
    """
    shifts_for: Dict[Tuple[int, int, int, int], np.ndarray] = {}
    rolled = 0
    for acq in acqs:
        if not is_epsi_acquisition(acq):
            continue
        key = layout_key(acq)
        nswitch, total = key[0], key[1]
        if key not in shifts_for:
            shifts_for[key] = switch_shifts(nswitch, total, slope)
            print(f"  {nswitch} switches of {total}: shifts {int(shifts_for[key].min())} to "
                  f"{int(shifts_for[key].max())} samples", file=sys.stderr)
        roll_acquisition(acq, shifts_for[key], nswitch, total)
        rolled += 1
    record_slope(header, slope)
    return rolled


# ---------- driver --------------------------------------------------------------------------


def check_stream(input_path: str, output_path: str = "", force: bool = False,
                 plot: bool = False) -> int:
    """
    Report what a converted stream carries, and write it back rolled straight when asked.
    Args:
        - input_path: the .mrd2 to read, or a FIFO carrying one
        - output_path: where to write the rolled stream, or "" to report only
        - force: roll a stream that already records a slope taken out of it
        - plot: draw the switch by position map of the measured layout
    Returns:
        - 0 when the file was read and reported on, 1 when it carries nothing to report on
    """
    header, items = read_stream(input_path)
    acqs = acquisitions_of(items)
    report_header(header)
    disagreements = report_census(header, items, acqs)

    layouts = group_layouts(acqs)
    if not layouts:
        sys.stdout.flush()
        print(f"\nno acquisition in {input_path} records an EPSI switch layout, so there is no "
              f"sampling window to check and nothing to shift", file=sys.stderr)
        if acqs:
            print(f"a spectral readout converts without one; an EPSI scan converted before "
                  f"MRStomrd2 recorded user_int has to be converted again", file=sys.stderr)
        return 1

    # the drift is measured from the biggest layout's real data: pooling the averaged prescans in
    # with the repetitions would measure neither of them
    measured_key, measured_acqs = next(iter(layouts.items()))
    data_acqs = [acq for acq in measured_acqs if not is_navigation(acq)] or list(measured_acqs)

    geometries: Dict[Tuple[int, int, int, int], dict] = {}
    drift = None
    for key, layout_acqs in layouts.items():
        report_layout(header, key, layout_acqs)
        nswitch, total, discard_pre, discard_post = key
        kept = total - discard_pre - discard_post
        pool = [acq for acq in layout_acqs if not is_navigation(acq)] or list(layout_acqs)
        cube, missing = acquisition_cube(pool, nswitch, total)
        if missing:
            print(f"  WARNING {missing} view/repetition slot(s) no acquisition filled, read as zero")
        geometry = echo_positions(cube, total, kept, discard_pre)
        geometries[key] = geometry
        if key == measured_key:
            drift = measure_echo_drift(cube, nswitch, total)
        report_echo(geometry, drift if key == measured_key else None, nswitch, total, kept)

    sys.stdout.flush()
    if disagreements:
        print(f"\nWARNING {disagreements} axis/axes disagree with the header's encoding limits, "
              f"which is what a reader sizes its k-space from", file=sys.stderr)
    if plot:
        name = getattr(header.measurement_information, 'measurement_id', None) or input_path
        plot_echo(geometries[measured_key], measured_key, str(name))

    if not output_path:
        print(f"\nreport only; pass --output to write this stream back with the echoes rolled "
              f"straight", file=sys.stderr)
        return 0

    already = header_user_double(header, SLOPE_PARAMETER)
    if drift is None or not drift['usable']:
        print(f"\nnothing to take out, copying the stream through unchanged", file=sys.stderr)
    elif already is not None and not force:
        print(f"\nthis stream already records {already:+.4f} samples per switch taken out of it, so "
              f"it is copied through unchanged; pass --force to roll it again", file=sys.stderr)
    else:
        print(f"\nrolling {drift['slope']:+.4f} samples per switch out of the stream",
              file=sys.stderr)
        rolled = shift_stream(header, acqs, drift['slope'])
        print(f"  rolled {rolled} acquisition(s)", file=sys.stderr)
    write_stream(output_path, header, items)
    print(f"wrote {output_path}", file=sys.stderr)
    return 0


def main() -> int:
    """
    Check a converted .mrd2, and with --output write it back with its echoes rolled straight.

    Every check that can stop the run lives here; past this point a stream that cannot be read is
    reported and the run ends with a status rather than a traceback
    Returns:
        - 0 when the file was read and reported on, 1 when it could not be read or carries no EPSI
          switch layout
    """
    parser = argparse.ArgumentParser(
        description="Check a converted MRD2 stream and optionally roll its EPSI echoes straight")
    parser.add_argument("-i", "--input", default=os.environ.get("INPUT_PIPE"),
                        help="the .mrd2 to check, or a FIFO carrying one (default: $INPUT_PIPE)")
    parser.add_argument("-o", "--output", default=os.environ.get("OUTPUT_PIPE"),
                        help="write the stream back with the echoes rolled straight (default: "
                             "$OUTPUT_PIPE; without either, this reports and writes nothing)")
    parser.add_argument("--force", action="store_true",
                        help="roll a stream that already records a drift taken out of it")
    parser.add_argument("--plot", action="store_true",
                        help="draw the switch by position map with the peak of each switch marked")
    args = parser.parse_args()

    if not args.input:
        parser.error("--input is required when $INPUT_PIPE is unset")
    if args.input != os.environ.get("INPUT_PIPE") and not os.path.exists(args.input):
        parser.error(f"{args.input} does not exist")

    try:
        return check_stream(args.input, args.output or "", args.force, args.plot)
    except RuntimeError as failure:
        # what the reader raises on a stream written against another schema, which is most of the
        # .mrd2 on disk: converted before the current mrd fork
        print(f"Cannot read {args.input}: {failure}", file=sys.stderr)
        print("this stream was written against a different mrd schema; convert the scan again with "
              "the current MRStomrd2", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
