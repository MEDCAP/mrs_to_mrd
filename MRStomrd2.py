"""
Convert MR Solutions .MRD raw data to MRD v2.

Three input modes, one per way a scan arrives:

    -t/--tar     Single experiment epsi folder wrapped in a tar file, where each subdirectory folder 
                 represent single repetition and pre-scan data that do not have the same dimension 
                 as rawdata. 

    -i/--input   A single .MRD filepath, which is how spectral fid data usually arrives: one file
                 already holds every repetition, so there is nothing to collect or group
    
    -f/--folder  Legacy method for local testing. One experiment folder, walked for its .MRD files

Both -t and -f represent one experiment and take their name from it, so mrs_organize reads the
directories inside only for their scan ids. Which files belong to one scan it decides from the
sequence and the acquisition matrix. All three converge on convert_group_to_mrd, so header choice
and repetition numbering have exactly one implementation.

Every file converts to acquisitions the same way whatever kind of scan it is. mrs_organize decides
which files make up a stream and whether that stream is data or an averaged prescan, and that answer
shows in the output name rather than in how the acquisitions are flagged on the way in.

A group hands over paths rather than parsed files, so convert_group_to_mrd reads one .MRD at a time
and releases it before the next. The repetition count the header needs comes from the group, which
recorded it when it probed the files to group them.

Where the usable sampling window sits inside an EPSI gradient switch, and whether the echo drifts
along the switch train, is reported by epsi_window.py rather than here - a diagnostic companion
that converts nothing, sharing is_epsi and switch_layout so its window is the same arithmetic as
what a real conversion writes.

-w stays here, because it is a question about a scan rather than about a file: it plots where the
echo peaks inside each gradient switch, from every repetition of one group pooled, so a peak that
walks along the switch train shows up against the window a conversion would keep. It then corrects
the files it read and plots the same thing again, which is what turns "the echo moves" into "and this
is what it looks like put right": the drift rolled out, and the readout moved onto the position the
sequence puts the echo at, which epsi_window derives from tramp and the sample period rather than
from the data. It groups its input exactly as a conversion does
and converts nothing: the roll lives in memory, and the conversion path reads its own copy of every
.MRD, so what -w reports on still converts from the samples as acquired. Rolling a converted stream
is mrd2shift.py's job.

Errors are raised only in main(). Past that point a file that cannot be converted is reported on
stderr and skipped, so one unreadable or stray file does not lose the rest of the scan.
"""

from __future__ import annotations

import argparse
import os
import sys
from itertools import product
from pathlib import Path
from typing import BinaryIO, Callable, Iterable, List, Optional, Sequence, Tuple

import numpy as np

# mrd python package
import mrd
import mrs_organize
from MRSreader import MRSdata
from mrs_organize import ScanGroup, is_prescan
from mrs_tar import read_scan_tar


def is_epsi(mrs: MRSdata) -> bool:
    """
    Whether a file's readout is split into gradient switches, read directly off the header rather
    than matched against the sequence name: no_switches only appears in the parameter block of a
    sequence that actually switched, so a file recording more than the unswitched default of one
    is what makes this a generic, name-agnostic check
    """
    return mrs.nswitch > 1


def switch_layout(mrs: MRSdata) -> Tuple[int, int, int]:
    """
    How one EPSI readout divides into gradient switches.

    no_switches and no_pts_switch are recorded separately from the sample count, so the acquisition
    holds more points per switch than the sequence calls usable: the extra ones are the gradient
    ramps either side of the flat top. nswitch is a divisor both here and in reconstruction, so a
    file recording 0 switches is read as 1 rather than being allowed to raise part way through a
    scan, and a file recording no usable width is read as keeping the whole switch rather than half
    of it
    Args:
        - mrs: one parsed MRS file of the epsi family
    Returns:
        - (nswitch, total points in one switch, points the sequence keeps per switch)
    """
    nswitch = max(mrs.nswitch, 1)
    total = mrs.nsamples // nswitch
    return nswitch, total, (mrs.npoints_per_switch or total)


def keep_mrs(mrs: MRSdata, name: str, accepted: Sequence[MRSdata] = ()) -> bool:
    """
    Whether a parsed file can join a group. Reported rather than raised, so that one bad file does
    not cost the rest of the scan.

    Grouping ran off a header-only probe, so this is where the full parse gets to disagree with it.
    A group holds two kinds of file - real data and the prescan that calibrates it - so the
    comparison has to be against a file of the same kind: a prescan disagreeing with the real data
    beside it is the group working as designed, not a mismatch to report
    Args:
        - mrs: the parsed file
        - name: what to call it on stderr, a path or a tar member name
        - accepted: files already accepted into this group, to find one of the same kind to compare
          dimensions against; empty for the first file of either kind
    Returns:
        - True when the file parsed
    """
    # MRSdata.read_from_file reports a parse failure and returns with rawdata left unset
    if mrs.rawdata is None or mrs.rawdata.size == 0:
        print(f"Skipping {name}: no raw data was read", file=sys.stderr)
        return False
    prescan = is_prescan(mrs.naverages, mrs.nrepetitions)
    reference = next((accepted_mrs for accepted_mrs in accepted
                      if is_prescan(accepted_mrs.naverages, accepted_mrs.nrepetitions) == prescan),
                     None)
    if reference is not None:
        # the fields grouping bucketed on, so a difference here is a difference grouping did not see
        differences = ", ".join(
            f"{attribute} {getattr(mrs, attribute)} vs {getattr(reference, attribute)}"
            for attribute in mrs_organize.SIGNATURE_FIELDS
            if getattr(mrs, attribute) != getattr(reference, attribute))
        if differences:
            print(f"WARNING {name} was grouped on a header reading {differences}, which its data "
                  f"block does not agree with; converting it anyway", file=sys.stderr)
    return True


def generate_acquisition(mrs: MRSdata, rep_base: int, rep_count: int) -> Iterable[mrd.StreamItem]:
    """
    Emit one acquisition per point of the encoding grid of one MRS file.

    Every axis of rawdata but the samples one is walked, each bounded by the MRSdata field that
    names it rather than by a position in rawdata.shape, so a file using an axis converts rather
    tha
        - EPSI split one repetition per file:    nrepetitions=1,  nviews=8 ->    8 acquisitions
        - EPSI acquired into a single file:      nrepetitions=N,  nviews=8 ->  N*8 acquisitions
        - spectral, repetitions on the nex axis: nrepetitions=40, nviews=1 ->   40 acquisitions
    rep_base is where this file's repetitions start within the group, so all three number
    identically from the reader's point of view.

    Each axis takes the MRD index that means the same thing
        nviews       -> kspace_encode_step_1 
        nsliceviews  -> kspace_encode_step_2
        nslices      -> slice
        nechoes      -> contrast
        nrepetitions -> repetition, offset by rep_base
    Args:
        - mrs: one parsed MRS file, rawdata indexed
               (nsamples, nviews, nsliceviews, nslices, nechoes, nrepetitions)
        - rep_base: repetition index this file's first repetition maps to
        - rep_count: repetitions in the whole group, for the LAST_IN_REPETITION flag
    Returns:
        - Iterable of mrd.StreamItem.Acquisition
    """
    epsi = is_epsi(mrs)
    prescan = is_prescan(mrs.naverages, mrs.nrepetitions)
    if epsi:
        # samples hold nswitch echoes, each npoints_per_switch long with a gradient ramp either
        # side
        # example cirrhrat_43_1: 1792 samples / 64 switches = 28, (28 - 12) / 2 = 8 points per ramp.
        # switch_layout is what -w reports against, so the window written here and the window that
        # report describes are the same arithmetic
        nswitch, points_per_switch, kept = switch_layout(mrs)
        discard = (points_per_switch - kept) // 2
    # one repetition's worth of acquisitions, which is every axis inside the repetition one. The
    # product walks them in the order the axes are listed, repetition slowest and view fastest, so a
    # repetition's acquisitions stay contiguous in the stream and views stay contiguous within it
    per_repetition = mrs.nechoes * mrs.nslices * mrs.nsliceviews * mrs.nviews
    grid = product(range(mrs.nrepetitions), range(mrs.nechoes), range(mrs.nslices),
                   range(mrs.nsliceviews), range(mrs.nviews))
    for counter, (irep, iecho, islice, isliceview, iview) in enumerate(grid):
        within = counter % per_repetition           # where this sits inside its own repetition
        acq = mrd.Acquisition()
        # MRS acquires on one channel, so add the coil axis: acq.data.shape=(coils=1, samples)
        acq.data = np.expand_dims(mrs.rawdata[:, iview, isliceview, islice, iecho, irep], axis=0)
        acq.head.acquisition_time_stamp_ns = np.uint64(mrs.acquisition_timestamp * 100) # 100ns -> ns units
        acq.head.sample_time_ns = mrs.sample_period * 100    # sample_period in units of 100ns
        acq.head.idx.average = mrs.naverages    # mrs already collapses averaged samples into single sample
        repetition = rep_base + irep
        acq.head.idx.repetition = repetition    # index of repetition
        acq.head.idx.kspace_encode_step_1 = iview
        acq.head.idx.kspace_encode_step_2 = isliceview
        acq.head.idx.slice = islice
        acq.head.idx.contrast = iecho           # index of echoes
        # unique and increasing across the whole group, since repetition already carries rep_base
        acq.head.scan_counter = repetition * per_repetition + within
        if repetition == 0 and within == 0:
            acq.head.flags |= mrd.AcquisitionFlags.FIRST_IN_REPETITION
        if repetition == rep_count - 1 and within == per_repetition - 1:
            acq.head.flags |= mrd.AcquisitionFlags.LAST_IN_REPETITION
        # both hold on a single view acquisition, so these are two ifs rather than if/elif
        if iview == 0:
            acq.head.flags |= mrd.AcquisitionFlags.FIRST_IN_PHASE
        if iview == mrs.nviews - 1:
            acq.head.flags |= mrd.AcquisitionFlags.LAST_IN_PHASE
        # the switch layout the discard points were worked out from. Only an EPSI readout has one,
        # and nswitch and points_per_switch are only defined for one, so this stays inside the
        # branch. It cannot ride on idx.contrast, which carries the sequence's own echo index
        if epsi:
            acq.head.user_int = [nswitch, points_per_switch]
            # if there is discard, encode it into acq.head
            if discard:
                acq.head.discard_pre = discard
                acq.head.discard_post = discard
        # prescan acquisitions stay in the same stream as the data they calibrate rather than a
        # file of their own, told apart downstream by this flag
        if prescan:
            acq.head.flags |= mrd.AcquisitionFlags.IS_NAVIGATION_DATA
        acq.phase = np.zeros(mrs.nsamples, dtype=np.float32)   # phase is not recorded
        yield mrd.StreamItem.Acquisition(acq)


def make_header(mrs: MRSdata, meas_id: str, rep_count: int) -> mrd.Header:
    """
    Fill in the MRD header from one file's parameters. Every file in a group was acquired at the
    same matrix, which is what let them be combined, so any of them describes the geometry
    Args:
        - mrs: the file the header describes, chosen by convert_group_to_mrd
        - meas_id: measurement id, e.g. cirrhrat_43_1
        - rep_count: repetitions in the group
    Returns:
        - mrd.Header
    """
    header = mrd.Header()

    subject = mrd.SubjectInformationType()
    subject.patient_id = meas_id            # e.g.) cirrhrat_43_1, KIC_Huh7msps5_08-15-2025.mrs
    header.subject_information = subject

    # t_r, t_e and flip_angle_deg are lists in this schema, one entry per value the sequence used,
    # and an MRS scan records a single one of each. Recorded only where the file carried a value,
    # for the reason tramp is below: writing the 0 default would read as a measured TR of zero
    # rather than as a record the file never held, and te is genuinely 0 on these sequences
    measured = {name: [float(value)]
                for name, value in (("t_r", mrs.tr), ("t_e", mrs.te),
                                    ("flip_angle_deg", mrs.flip_angle)) if value}
    if measured:
        header.sequence_parameters = mrd.SequenceParametersType(**measured)

    meas = mrd.MeasurementInformationType()
    meas.sequence_name = mrs.sequence_name
    meas.measurement_id = meas_id
    meas.protocol_name = meas_id.split("_")[0]
    meas.relative_table_position = mrd.ThreeDimensionalFloat(x=mrs.FOVoffset[0] * 1e3,
                                                            y=mrs.FOVoffset[1] * 1e3,
                                                            z=mrs.FOVoffset[2] * 1e3)  # m -> mm
    header.measurement_information = meas

    header.experimental_conditions.h1resonance_frequency_hz = mrs.base_frequency

    # tramp is the readout gradient ramp time in us. It is what places the EPSI sampling window:
    # the ramp is the leading stretch of a switch that is not yet on the gradient plateau, so a
    # reconstruction needs it to work out which points of each switch are usable. Recorded only when
    # the file carried it, since writing the 0 default would read as a measured ramp of zero.
    # UserParameterLongType is this schema's integer parameter, there is no int arm of its own
    if mrs.tramp:
        if header.user_parameters is None:
            header.user_parameters = mrd.UserParametersType()
        header.user_parameters.user_parameter_long.append(
            mrd.UserParameterLongType(name="tramp", value=int(mrs.tramp)))

    encoded_space = mrd.EncodingSpaceType()
    # the encoded matrix is the k-space one, so its third axis is the second phase encode rather
    # than the slice count: slices are separate acquisitions carrying an index, not an encoded axis
    encoded_space.matrix_size = mrd.MatrixSizeType(x=mrs.nsamples, y=mrs.nviews, z=mrs.nsliceviews)
    encoded_space.field_of_view_mm = mrd.FieldOfViewMm(x=mrs.FOV * 1e3, y=mrs.FOV * 1e3, z=0)

    # each limit is the size of one rawdata dimension, as a maximum index, and every dimension
    # generate_acquisition walks has one, so a reader can size any axis it finds indexed
    limits = mrd.EncodingLimitsType()
    limits.kspace_encoding_step_0 = mrd.LimitType(maximum=mrs.nsamples - 1)
    limits.kspace_encoding_step_1 = mrd.LimitType(maximum=mrs.nviews - 1)
    limits.kspace_encoding_step_2 = mrd.LimitType(maximum=mrs.nsliceviews - 1)
    # reconstruction sizes its k-space off the phase limit, so it has to carry the view count too
    limits.phase = mrd.LimitType(maximum=mrs.nviews - 1)
    limits.slice = mrd.LimitType(maximum=mrs.nslices - 1)
    limits.contrast = mrd.LimitType(maximum=mrs.nechoes - 1)     
    # clamped, so that a file reporting no repetitions still gets a valid header
    limits.repetition = mrd.LimitType(minimum=0, maximum=max(rep_count - 1, 0))

    encoding = mrd.EncodingType()
    encoding.encoded_space = encoded_space
    encoding.encoding_limits = limits
    header.encoding.append(encoding)
    return header


def read_mrs_group(filepaths: Sequence[str],
                   load: Optional[Callable[[MRSdata, str], None]] = None) -> List[MRSdata]:
    """
    Parse a group of .MRD files, dropping the ones that cannot contribute
    Args:
        - filepaths: paths (or tar member names) in acquisition order
        - load: fills one MRSdata from one entry of filepaths; defaults to reading from disk, tar
          mode passes one that parses an in-memory buffer instead
    Returns:
        - parsed MRSdata in the same order, possibly shorter than filepaths
    """
    load = load or (lambda mrs, filepath: mrs.read_from_file(filepath))
    mrs_list: List[MRSdata] = []
    for filepath in filepaths:
        mrs = MRSdata()                                 # one instance per file, they are all kept
        load(mrs, filepath)
        if not keep_mrs(mrs, str(filepath), mrs_list):
            continue
        mrs_list.append(mrs)
    return mrs_list


def convert_group_to_mrd(group: ScanGroup, output: BinaryIO,
                         load: Optional[Callable[[MRSdata, str], None]] = None) -> bool:
    """
    Write one group - the real data of one experiment plus the prescans beside it - as one MRD
    stream.

    One file is parsed at a time and released before the next is read. 

    The rawdata paths convert first, so the header is built from an acquisition rather than from a
    prescan, and the prescans follow into the same stream numbered after them. Nothing in the
    acquisitions says which list a file came from - they are told apart by IS_NAVIGATION_DATA
    Args:
        - group: the files to convert (real data and prescan) and the meas_id to record
        - output: writable binary stream. Must be a file object, not a path: BinaryMrdWriter
                  special-cases str only, so a Path would be mistaken for a stream
        - load: fills one MRSdata from one of the group's paths; defaults to reading from disk, tar
          mode passes one that parses an in-memory buffer instead
    Returns:
        - True when a stream was written
    """
    load = load or (lambda mrs, filepath: mrs.read_from_file(filepath))
    filepaths = group.rawdata_file_list + group.prescan_file_list
    rep_count = group.stream_repetitions
    print(f"Converting {group.meas_id}: {len(filepaths)} files, {rep_count} repetitions",
          file=sys.stderr)

    # opened on the first file that parses rather than up front, since a writer closed without a
    # header would leave an unreadable stream behind. Closing it is what emits the end-of-stream
    # sentinel, hence the finally
    writer: Optional[mrd.BinaryMrdWriter] = None
    rep_base = 0
    try:
        # first convert raw data
        for filepath in group.rawdata_file_list:
            mrs = MRSdata()                     # one at a time, released once written
            # through the loader rather than off disk directly, which is the seam the whole function
            # is built around: tar mode passes one that parses an in-memory buffer, and -w passes one
            # that applies its echo correction as each file is read
            load(mrs, filepath)

            if writer is None:
                writer = mrd.BinaryMrdWriter(output)
                writer.write_header(make_header(mrs, group.meas_id, rep_count))
            writer.write_data(generate_acquisition(mrs, rep_base, rep_count))
            rep_base += mrs.nrepetitions
    finally:
        if writer is not None:
            writer.close()
    if writer is None:
        print(f"No data to convert for {group.meas_id}", file=sys.stderr)
        return False
    return True


def shift_echo_position(mrs_list: List[MRSdata], name: str = "", slope: float = 0.0,
                        anchor: Optional[int] = None) -> Optional[dict]:
    """
    Where the echo peaks inside each gradient switch, for one scan's repetitions pooled.

    nsamples of an EPSI readout is nswitch stretches of nsamples // nswitch points, and the
    gradient echo sits at one position of each stretch. The reconstruction reads every switch at
    the same offset, so the question this answers is whether that one offset can be right: a peak
    that sits at the same position in all nswitch switches says yes, a peak that walks along the
    train says no single window suits the whole readout.

    Signal is summed over the views and over every repetition of every file in the group, which is
    what gives one switch enough to peak on: one repetition of one view is mostly noise. That
    pooling is also why this takes a list rather than a file, and why an averaged phantom is
    reported as a scan of its own rather than added in - it would be pooling a different scan.

    This reads the data it is given and changes nothing, which is what lets correct_echo_position
    call it twice - once on the readout as acquired and once after rolling the drift out - and have
    the two blocks and the two figures be comparable. For the sampling window itself, and for the
    drift measured rather than read off an argmax, see epsi_window.py
    Args:
        - mrs_list: the parsed files of one scan, read rather than probed, in acquisition order
        - name: what to call this scan in the printed block and the figure title, e.g. a meas_id
        - slope: the drift still in these samples, so the two echoes can be told apart against the
          line it draws. Zero once the drift has been rolled out, which is why the two figures of one
          scan are drawn with different values of it
        - anchor: where switch 0's readout echo sits, or None to read it off switch 0's own peak
    Returns:
        - dict of the pooled signal, the peak position per switch, how many switches peak inside the
          kept window, the readout profile with its peak/median and the switch layout, or None when
          the list held no EPSI readout to profile
    """
    # epsi_window imports is_epsi and switch_layout from here, so this import is function-level to
    # keep that from being circular. matplotlib is deferred for the reason it is everywhere else in
    # this file: converting needs no plotting stack
    import epsi_window
    import matplotlib.pyplot as plt

    label = name or "scan"
    signal, layout, pooled, reference = None, None, 0, None
    # loop through raw data mrs_list to identify the peak
    for mrs in mrs_list:
        # switch_cube maps rawdata (nsamples, nviews, ...) to (nswitch, nsamples // nswitch, nviews,
        # everything else), e.g. 1280 samples of 12 views -> 64 switches of 20 positions, the 20
        # including the gradient ramp either side of the flat top. Summed to magnitude over the views
        # and repetitions, since one repetition of one view is mostly noise
        cube = np.abs(epsi_window.switch_cube(mrs)).sum(axis=(2, 3))         # (nswitch, total)
        if signal is not None and cube.shape != signal.shape:
            print(f"WARNING leaving a {cube.shape} readout of {label} out of a {signal.shape} "
                  f"pool", file=sys.stderr)
            continue
        signal = cube if signal is None else signal + cube
        layout = layout or switch_layout(mrs)
        # any pooled file describes the sequence, since they all share the acquisition matrix
        reference = reference or mrs
        pooled += 1
    if signal is None:
        print(f"No EPSI readout in {label} to find an echo position in", file=sys.stderr)
        return None

    nswitch, total, kept = layout
    discard_pre = (total - kept) // 2       # the window generate_acquisition writes into the header
    peaks = np.argmax(signal, axis=1)
    switches = np.arange(nswitch)
    # the window the sequence actually samples: from the end of the leading ramp through the flat top.
    # That is what the peaks are worth reading against, since it is where the signal was meant to be,
    # and it is not where a conversion currently keeps its points - the difference between the two is
    # the mrd2recon --pad that would move a reconstruction onto the plateau
    window = epsi_window.sequence_window(reference) if reference is not None else None
    start = window[0] if window else discard_pre
    # counted the same cyclic way generate_acquisition addresses a window, so one running off the end
    # of the switch is not mistaken for one that peaks outside it
    inside = int((((peaks - start) % total) < kept).sum())
    # which of the two echoes each switch peaked on. A switch peaks on whichever of them was brightest
    # in it, so without this split the marks read as scatter on a train that is perfectly orderly
    positions = epsi_window.expected_echo_positions(reference) if reference is not None else None
    families = epsi_window.peak_families(signal, slope, anchor,
                                         fallback_separation=(positions[1] - positions[0]
                                                              if positions else None))
    readout_family = families['readout_family']
    # the readout profile, which is every switch summed on top of each other, and how sharp it is.
    # This is what a drifting echo smears and what rolling the drift out restores, and it is the
    # number that moves when the per switch argmax does not: an argmax can land on a second feature
    # in a few switches whether the echo drifts or not
    profile = signal.sum(axis=0)
    median = np.median(profile)
    sharpness = float(profile.max() / median) if median else float('inf')

    sampled = (f"the sequence samples {start}..{start + kept - 1}" if window
               else f"a conversion keeps {start}..{start + kept - 1}")
    print(f"\n{label}: {pooled} file(s) pooled, {nswitch} switches of {total} positions, "
          f"{sampled}")
    print(f"  the peaks visit positions {int(peaks.min())} to {int(peaks.max())}, at {int(peaks[0])} "
          f"in the first switch and {int(peaks[-1])} in the last, "
          f"{int(readout_family.sum())} of them on the readout echo and "
          f"{int((~readout_family).sum())} on the second")
    print(f"  {inside} of {nswitch} switches peak inside that window, and pooled over every "
          f"switch the profile peaks at {int(np.argmax(profile))} with peak/median "
          f"{sharpness:.2f}")
    # no slope is printed on purpose. A position is cyclic within a switch and the argmax jumps
    # between whichever features are momentarily strongest - on cirrhrat_43_1 two of them 14 apart -
    # so a line through these points describes neither. epsi_window.measure_echo_drift searches the
    # correction and scores it, which is the measurement worth quoting
    print(f"  these are argmaxes; the drift itself is searched and scored by "
          f"epsi_window.measure_group_drift")

    figure, axes = plt.subplots(figsize=(9, 7))
    axes.imshow(signal, aspect='auto', origin='lower', interpolation='nearest',
                extent=(-0.5, total - 0.5, -0.5, nswitch - 0.5))
    # the span the sequence samples, so the peaks are read against where the signal was meant to be
    # rather than in the abstract
    axes.axvspan(start - 0.5, start + kept - 0.5, color='w', alpha=0.15,
                 label=f"{'sequence samples' if window else 'kept window'} "
                       f"{start}..{start + kept - 1}")
    if window and discard_pre != start:
        # where a conversion keeps its points instead, and so what a --pad would have to make up
        for edge in (discard_pre, discard_pre + kept - 1):
            axes.axvline(edge, color='w', ls=':', lw=1, alpha=0.6)
        axes.plot([], [], color='w', ls=':', lw=1, alpha=0.6,
                  label=f"a conversion keeps {discard_pre}..{discard_pre + kept - 1}, "
                        f"--pad {discard_pre - start}")
    # split by echo, so an outlying mark reads as "that is the rephasing echo" rather than as scatter
    axes.plot(peaks[readout_family], switches[readout_family], 'x', color='C3', ms=6,
              label=f"readout echo ({int(readout_family.sum())} switches)")
    if not readout_family.all():
        axes.plot(peaks[~readout_family], switches[~readout_family], '+', color='C1', ms=7,
                  label=f"second echo ({int((~readout_family).sum())} switches)")
    axes.set_xlabel(f"position within the {total} point switch")
    axes.set_ylabel("switch")
    axes.set_title(f"{label}: echo position per switch, {pooled} file(s) pooled")
    axes.legend(fontsize=8, loc='upper right')
    figure.colorbar(axes.images[0], ax=axes, label='signal summed over views and repetitions')
    figure.tight_layout()
    plt.show()
    return dict(signal=signal, peaks=peaks, nswitch=nswitch, total=total, kept=kept,
                discard_pre=discard_pre, pooled=pooled, inside=inside, profile=profile,
                sharpness=sharpness, window=window, families=families,
                readout_family=readout_family)


def correct_echo_position(group: Optional[ScanGroup],
                          load: Optional[Callable[[MRSdata, str], None]] = None,
                          slope: Optional[float] = None,
                          output: Optional[Path] = None) -> int:
    """
    Show where the echo sits in each gradient switch of one experiment's data, and what correcting it
    does: the drift along the switch train taken out, and the readout moved onto the position the
    sequence puts the echo at. Converts nothing.

    Only the real data of the group. The averaged prescan beside it calibrates a reconstruction rather
    than being reconstructed, so where its echo sits is nobody's question - and it could not have
    been pooled in with the data anyway, since averaged samples placed beside 27 repetitions measure
    neither of them.

    Every repetition of every file goes into one measurement, which is what gives the drift enough
    signal to be found, and the dataset gets two blocks and two figures - the readout as acquired and
    the readout corrected - so the correction reads as a before and after rather than as a claim.

    The roll happens on the parsed files in memory and goes no further unless `output` names a file to
    write, which is what makes the correction testable: a shifted stream can be reconstructed and the
    Lorentzian fit of the result is the honest verdict on whether the shift was right. Without it
    nothing is written, and the conversion path reads its own copy of each .MRD, so a scan reported on
    here converts from the samples as acquired exactly as it did before. Worth knowing when that path
    is finished: the
    constant part of what epsi_window.echo_alignment works out is the same quantity as the per tramp
    prepend it hard-codes, derived from the sequence rather than tabulated, so conversion wants one
    or the other and never both. For correcting a stream that has already been converted, see
    mrd2shift.py
    Args:
        - group: the experiment to report on, from mrs_organize, or None when nothing grouped
        - load: how to turn one of a group's paths into a parsed MRSdata; see read_mrs_group
        - slope: roll rate in samples per switch to apply instead of the measured one, from --slope,
          for correcting a scan case by case. Zero takes no drift out but still makes the constant
          move; None leaves the measurement in charge
        - output: where to write the corrected stream, or None to report and write nothing
    Returns:
        - 1 when the data was reported on, 0 when there was nothing to report on
    """
    # imported here rather than at module scope for two reasons: epsi_window imports this module, so
    # a module scope import back would be circular, and it pulls in a plotting stack that converting
    # has no need of
    import epsi_window

    if group is None or not group.rawdata_file_list:
        return 0
    label = group.meas_id
    # every file is held at once, because the drift is measured from all of their repetitions pooled:
    # one repetition of 12 views cannot place the echo per switch
    mrs_list = read_mrs_group(group.rawdata_file_list, load)

    # measured before anything is plotted, because the figure needs the slope to tell the two echoes
    # of a switch apart: against a flat line a drifting readout's peaks smear across the switch and
    # the front of the train gets labelled as the second echo. Whether the drift is worth acting on
    # is measure_echo_drift's decision and nothing here second-guesses it; the constant move is
    # independent of it, so a readout with no drift is still put where it belongs
    alignment = epsi_window.echo_alignment(mrs_list, slope=slope)
    if alignment is None:
        # no EPSI readout to align, but the map is still worth drawing
        return 1 if shift_echo_position(mrs_list, f"{label} as acquired") is not None else 0

    # nothing is applied before this figure: it is the readout exactly as the scanner wrote it, only
    # labelled with what is about to be taken out of it
    as_acquired = shift_echo_position(mrs_list, f"{label} as acquired, no correction applied",
                                      slope=alignment['applied'],
                                      anchor=alignment['families']['anchor'])
    if as_acquired is None:
        return 0
    epsi_window.shift_report(alignment, label=label)

    shifts = alignment['shifts']
    if np.any(shifts):
        for mrs in mrs_list:
            if is_epsi(mrs):
                epsi_window.shift_rawdata(mrs, shifts)
        # the same measurement over the same files, so the peaks, the count inside the sampled window
        # and the profile peak/median are read against the block printed above. The slope is zero
        # here: it has been taken out of the samples, so there is none left in the residuals
        expected = alignment['expected']
        corrected = (f"{label} corrected: {alignment['applied']:+.4f} per switch"
                     + (f" {'given' if alignment['given'] else 'measured'}")
                     + (f", echo moved onto {expected}" if expected is not None else ""))
        shift_echo_position(mrs_list, corrected, slope=0.0,
                            anchor=alignment['expected'] if expected is not None else None)

    if output is not None:
        # written through the same conversion the real path uses, with the shift applied as each file
        # is read: convert_group_to_mrd already takes a loader for the tar case, so the header, the
        # acquisitions and the repetition numbering are exactly what a conversion would write
        read = load or (lambda mrs, filepath: mrs.read_from_file(filepath))

        def shifting_load(mrs: MRSdata, filepath: str) -> None:
            read(mrs, filepath)
            if is_epsi(mrs) and np.any(shifts):
                epsi_window.shift_rawdata(mrs, shifts)

        print(f"Writing the corrected stream to {output}", file=sys.stderr)
        with open(output, "wb") as stream:
            convert_group_to_mrd(group, stream, load=shifting_load)
    return 1


def convert_folder_to_mrd(folder: Path,
                          meas_id_override: str = "",
                          dry_run: bool = False,
                          window_check: bool = False,
                          slope: Optional[float] = None,
                          output: Optional[Path] = None) -> bool:

    """
    Walk one experiment folder for .MRD files and convert each scan to its own stream inside it
        spectral: {experiment_folder}/{KIC_Huh7msps5_08-15-2025.MRD}
            -> KIC_Huh7msps5_08-15-2025.mrs/KIC_Huh7msps5_08-15-2025.mrs_1puls_extrf_KIC.mrd2
        EPSI:     {experiment_folder}/{modal}/{scan_id}/{24804_000_0.MRD}, one repetition per scan directory
            -> cirrhrat_43_1/cirrhrat_43_1_epsigre_combined.mrd2
        EPSI:     a subdirectory among those has navg>1 and nrep=1 instead - an averaged prescan,
                  converted into the same stream as the data it calibrates rather than a file of
                  its own, its acquisitions flagged IS_NAVIGATION_DATA
    Args:
        - folder: the experiment folder to walk
        - meas_id_override: measurement id to record instead of the folder's name
        - dry_run: report the grouping and what looks wrong with it, converting nothing
        - window_check: report the sampling window and correct the echo position, converting nothing
        - slope: with window_check, the roll rate in samples per switch to apply instead of the
          measured one
    Returns:
        - True when at least one scan was written, or when a dry run found something to convert
    """
    # scan subfolders of provided folder and group them into rawdata and prescan data as object ScanGroup
    groups = mrs_organize.organize_folder(folder, meas_id_override)
    if dry_run:
        mrs_organize.report(groups)
        return bool(groups)
    # -w is a question about the scans rather than a conversion, so it answers and returns. The roll
    # it applies lives in the parsed files it read, and nothing downstream of here sees it - unless
    # --output names a file, which writes the corrected samples as a stream to reconstruct
    if window_check:
        return bool(correct_echo_position(groups, slope=slope, output=output))
    written = False
    
    # convert .MRD file in groups to mrd.acquisition
    rawdata_across_repetition = None 
    for i, filepath in enumerate(groups.rawdata_file_list):
        mrs = MRSdata()
        mrs.read_from_file(filepath)
        # shift the echo center for EVO2 EPSI sequences only as the gradient is on before the first sample point
        # discard last n points and prepend n zeropoints at the beginning of nsamples dimension
        if mrs.sequence_name == 'epsigre':
            # If EVO2 epsi on tramp=100us, readout grad started 7points before the first sample e.g.) kidney data
            if mrs.tramp == 100:
                print(f'Detected EVO2{mrs.sequence_name} and tramp={mrs.tramp}, prepending 5 zeros', file=sys.stderr)
                nprepend_zeros = 5
            # If EVO2 epsi on tramp=112us, readout grad started 5points before the first sample e.g.) cirrhrat data
            elif mrs.tramp == 112:
                print(f'Detected EVO2{mrs.sequence_name} and tramp={mrs.tramp}, prepending 7 zeros', file=sys.stderr)
                nprepend_zeros = 7
            if nprepend_zeros:
                mrs.rawdata = np.concatenate((np.zeros(((nprepend_zeros,) + mrs.rawdata.shape[1:])),mrs.rawdata[:-nprepend_zeros,...]),axis=0)
                print(f"shape after prepend zeros={mrs.rawdata.shape}")

    return True
    
def convert_file_to_mrd(input_path: Path, output_path: Optional[Path] = None,
                        meas_id_override: str = "") -> bool:
    """
    Convert a single .MRD file. Spectral fid data arrives this way: one file already holds every
    repetition on its nex axis, so there is nothing to collect or group
    Args:
        - input_path: the .MRD file
        - output_path: where to write, defaulting to the name its scan would have been given
        - meas_id_override: measurement id to record instead of the directory holding the file,
          which is what names a file sitting outside any experiment folder
    Returns:
        - True when a stream was written
    """
    groups = mrs_organize.organize_folder(input_path, meas_id_override)
    if not groups:
        print(f"Nothing to convert in {input_path}", file=sys.stderr)
        return False
    group = groups[0]
    destination = output_path or Path(group.output_path)
    print(f"Converting {input_path} to {destination}", file=sys.stderr)
    with open(destination, "wb") as output:
        return convert_group_to_mrd(group, output)


def convert_tar_to_mrd(tar_path: Path, output_path: Path, meas_id_override: str = "") -> bool:
    """
    Convert tar filed directory of single experiment epsi MRS data folder

    Tyger hands a job its input buffer as a named FIFO, which is strictly sequential, so the
    archive is read forward once and held in memory: the parameter block sits after the raw data at
    EOF, the .SPR sidecar can follow the .MRD members in tar order, and the header needs the whole
    group before the first acquisition can be written.

    A Tyger job has one output buffer, so this writes one stream. The members are grouped exactly
    as a folder is, which turns "the caller tarred one scan" from an assumption into a checked
    precondition: more than one scan in the archive is an error, unless --output names a directory
    to write them all into
    Args:
        - tar_path: tar archive of one scan directory, or a FIFO carrying one
        - output_path: where to write the stream, which may also be a FIFO, or a directory when the
          archive holds more than one scan
        - meas_id_override: measurement id to record instead of the tar's root directory name
    Returns:
        - True when a stream was written
    """
    # opening a FIFO for read blocks until the buffer sidecar opens the write end
    with open(tar_path, "rb") as tar_stream:
        tar_meas_id, spr_frequency, members = read_scan_tar(tar_stream)
    fallback = meas_id_override or tar_meas_id or Path(output_path).stem
    if not (meas_id_override or tar_meas_id):
        # nothing to name the scan after: the archive holds no single root directory
        print(f"No single root directory in the tar, calling this scan {fallback}", file=sys.stderr)
    if not spr_frequency:
        print("No base frequency from a .SPR sidecar in the tar", file=sys.stderr)

    payloads = dict(members)
    groups = mrs_organize.organize_members(members, fallback_meas_id=fallback)
    report_warnings(groups)
    into_directory = output_path.is_dir()
    if len(groups) > 1 and not into_directory:
        raise ValueError(f"the tar holds {len(groups)} scans "
                         f"({', '.join(group.output_name for group in groups)}) but --output names "
                         f"a single stream: tar one scan directory, or point --output at a "
                         f"directory")

    def load(mrs: MRSdata, name: str) -> None:
        mrs.parse_from_buffer(payloads[name])
        mrs.set_base_frequency(spr_frequency)   # the .MRD may defer its frequency to the sidecar

    written = False
    for group in groups:
        # meas_id_override names the scan itself, so it wins over what grouping inferred
        group.meas_id = meas_id_override or group.meas_id or fallback
        destination = output_path / group.output_name if into_directory else output_path
        # opened even with nothing to write: on Tyger this is the output buffer's FIFO, and closing
        # it empty is what tells the sidecar the job produced nothing, rather than leaving it
        # blocked on a stream that never opens
        with open(destination, "wb") as output:
            written |= convert_group_to_mrd(group, output, load=load)
    if not groups and not into_directory:
        with open(output_path, "wb"):       # nothing grouped, but the buffer still has to close
            pass
    return written




def main() -> int:
    """
    Convert MRS data to MRD2, in exactly one of three input modes, or with -w report where the echo
    peaks inside each gradient switch of what those modes resolve to and convert nothing.

    Every check that can stop the run lives here. Past this point a file that cannot be converted is
    reported and skipped, so the run ends with a status rather than a traceback
    Returns:
        - 0 when at least one stream was written, or one scan was reported on with -w. 1 when
          nothing was
    """
    parser = argparse.ArgumentParser(description="Convert MR Solutions MRS data to MRD2 format")
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("-t", "--tar", type=Path,
                      help="tar archive of one scan directory, or a FIFO carrying one")
    mode.add_argument("-i", "--input", type=Path,
                      help="single MRS .MRD file")
    mode.add_argument("-f", "--folder", type=Path,
                      help="directory to walk for MRS .MRD files")
    parser.add_argument("-o", "--output", type=Path,
                        help="file or FIFO to write the MRD2 stream to. Required with --tar "
                             "(default: $OUTPUT_PIPE), optional with --input, and with --folder -w "
                             "writes the corrected stream so it can be reconstructed")
    parser.add_argument("--meas-id", default="",
                        help="measurement id to record instead of the name of the folder, tar root "
                             "or directory the data came from")
    parser.add_argument("-n", "--dry-run", action="store_true",
                        help="with --folder only: report how the files group them, without converting")
    parser.add_argument("-w", "--window", action="store_true",
                        help="with --folder only: plot where the echo peaks in each gradient switch, "
                             "then take the drift out, move the readout onto the position the "
                             "sequence puts the echo at, and plot it again, converting nothing")
    parser.add_argument("--slope", type=float, default=None, metavar="PER_SWITCH",
                        help="with -w only: roll every switch at this rate in samples per switch, "
                             "instead of the rate measured from the data. 0 takes no drift out but "
                             "still moves the readout onto the sequence's echo position")
    args = parser.parse_args()

    if args.tar and not args.tar.exists():
        parser.error(f"{args.tar} does not exist")
    if args.input and not args.input.is_file():
        parser.error(f"{args.input} is not a file")
    if args.folder and not args.folder.is_dir():
        parser.error(f"{args.folder} is not a directory")
    if args.tar and (args.window or args.dry_run):
        parser.error(f"-t can only used to convert file")
    if args.slope is not None and not args.window:
        parser.error("--slope sets the roll rate -w applies, so it needs -w")

    if args.folder:
        # convert folder allows
        print(f"Converting folder of single experiment folder {args.folder}", file=sys.stderr)
        written = convert_folder_to_mrd(args.folder, args.meas_id, args.dry_run, args.window,
                                        args.slope, args.output)
    elif args.tar:
        output = args.output or Path(os.environ.get("OUTPUT_PIPE", ""))
        if not str(output):
            parser.error("--tar needs --output, or $OUTPUT_PIPE set")
        print(f"Converting tar of single experiment folder {args.tar}", file=sys.stderr)
        written = convert_tar_to_mrd(args.tar, output, args.meas_id)
    else:
        print(f"Converting single file {args.input}", file=sys.stderr)
        written = convert_file_to_mrd(args.input, args.output, args.meas_id)
    if not written:
        print("Nothing was converted", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
