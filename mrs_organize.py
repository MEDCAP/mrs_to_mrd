"""
Context:
EPSI requires organizing a series of folders of .MRD file, each representing single repetition.
Hence, the acquisition in the same experiment folder should be wrapped as a single raw file.

Raw data folder structure:
MR Solutions writes one directory per acquisition, named with the scanner's numeric scan id and
holding <id>_000_0.MRD next to its .SPR, .dat and .SQL. What sits above that directory varies:

    EPSI        experiment_name/[anything/]scan_id/raw file
                example: cirrhrat_43_1/epsi/24804/24804_000_0.MRD
                one repetition per file, so the scan directories of one experiment combine

    spectral    KIC_Huh7msps5_08-15-2025.mrs/KIC_huh7_5.MRD
                one file holding all its repetitions on the nrepetitions axis, no scan directory

There is nothing here that reads which of those two shapes a file arrived in: every file in one
experiment folder that shares a sequence name converts into one stream together, whatever kind of
scan it is. That's what makes the converter generic rather than needing to special case EPSI.

How to combine:
Every input names exactly one experiment - -f is pointed at one folder, a tar holds one, -i is a
lone file - so there is one group to return rather than a list of them. Nothing is keyed on to
decide which group a file joins, because there is only ever the one: a file is either real
acquisition data that concatenates on the repetition axis, or the averaged prescan beside it, and
both convert into the same stream.

That is why group_experiment returns a single ScanGroup, holding two lists of paths:
rawdata_file_list, the real acquisition data, and prescan_file_list, the averaged calibration scan
beside it (naverages>1 at a single repetition, MR Solutions' way of recording a scan that was
averaged rather than repeated). convert_group_to_mrd converts rawdata_file_list first, so it builds
the header from real data, and prescan_file_list after it into the same stream rather than a file of
their own - they're told apart downstream by IS_NAVIGATION_DATA, not by which list or file they came
from.

Pointing -f at a folder of several experiments therefore converts them into one stream rather than
one each, which is what the guards below are for: it is the caller's job to name one experiment.

A group holds paths and the few scalars a conversion needs, never parsed files: the dimensions are
read once by a header-only probe to decide the grouping, and everything past that point re-reads
one file at a time. That is what lets a 27 file series convert without ever holding more than one
of them in memory, and it is why nrepetitions is recorded here rather than recomputed downstream.
rawdata_shape is recorded for the same reason: the probe already knows the axes one file is reshaped
onto, so the shape the whole series comes to is known before a single file is read back.

Which files are real data and which are prescan is decided structurally, by the acquisition matrix.
The two shapes a scan arrives in are told apart by nothing more than how many paths the rawdata list
ends up holding:

    case 1  many files, each nrepetitions=1, one per scan directory - the EPSI series, combined
    case 2  one file whose nrepetitions axis already holds them all - spectral, and EPSI acquired
            into a single file

Both cases resolve to one repetition count without either being recognised, by summing nrepetitions
over the rawdata files: case 1 contributes 1 per file, so the count is how many files arrived, and
case 2 contributes the whole axis from the one file. That sum is ScanGroup.nrepetitions, and it is
the only repetition count a group carries. The prescans are outside it - an averaged calibration
scan is not a repetition of the acquisition it calibrates, so however many prescan files a group
holds, none of them lengthens the repetition axis the header declares.

Every rawdata file must share one acquisition matrix, since the whole point of the group is that
they concatenate on the repetition axis; a second distinct matrix among them is misfiled data and
group_experiment raises rather than silently lengthening the series with it. That raise is the check
standing where the grouping key used to: with one group per input there is nowhere for an odd file
to go, so it has to be said rather than absorbed. The prescan list carries no such rule - a scan can
be calibrated by any number of differently shaped prescans, and none of them is concatenated onto
anything.

experiment_name is not read out of the path, it is what the caller pointed at: -f names one
experiment folder and a tar holds one, so its root is the experiment. Whatever sits between the
experiment and the scan directory, and whether anything sits there at all, never has to be
recognised - which is the point, since that level is named for the sequence, the nucleus or the
contrast and no list of those names stays complete.

Only -i names a lone file with no folder to take, and then the path is walked up past the scan id
the scanner always writes. If nothing above it names the experiment, use the tar filename

A tar names one experiment the same way a folder does, and read_scan_tar is what turns one into the
members organize_members groups. Tyger hands a job its input buffer as a named FIFO, which is
strictly sequential: lseek() on it returns ESPIPE. That rules out the three filesystem lookups the
folder walk relies on - the .SPR sidecar found via Path.parent.iterdir(), the directory names
meas_id is derived from, and the grouping of files belonging to one scan - so the experiment
directory is wrapped in a tar to put all three into one sequential stream:

    tar cf - -C /data cirrhrat_0_1 | tyger buffer write $input_buffer

One tar is one experiment and converts to one stream, exactly as one folder does, so it is the
caller's job to tar exactly one experiment directory

The stream is named for what it holds:
    <experiment>_<sequence>_prescan.mrd2    every file in the group is an averaged prescan
    <experiment>_<sequence>_combined.mrd2   more than one file converted into one acquisition
    <experiment>_<sequence>_<scan_id>.mrd2  a single file converting on its own
"""

from __future__ import annotations

import argparse
import os
import re
import sys
import tarfile
from dataclasses import dataclass, field
from pathlib import Path, PurePosixPath
from typing import BinaryIO, Callable, List, NamedTuple, Optional, Sequence, Tuple

from MRSreader import MRSdata

MRD_SUFFIX = ".MRD"
SPR_SUFFIX = ".SPR"
OUTPUT_SUFFIX = ".mrd2"


class Signature(NamedTuple):
    """
    The acquisition matrix of one file: everything that distinguishes one kind of scan from another.

    Two files whose signatures differ cannot concatenate on the repetition axis - a calibration scan
    at 2176x12 beside a 1792x8 series, or a file that already carries 25 repetitions of its own. The
    leading fields are the axes rawdata is reshaped onto, in that order
    (MRSreader.MRSdata.rawdata), so a difference reads in the same order as the shape it describes
    """
    nsamples: int
    nviews: int
    nsliceviews: int
    nslices: int
    nechoes: int
    nrepetitions: int
    naverages: int
    datatype: int


def signature_of(mrs: MRSdata) -> Signature:
    """
    The acquisition matrix of one probed or fully parsed file
    Args:
        - mrs: the file, probed or fully parsed
    Returns:
        - its Signature
    """
    # a probe that failed on a malformed header leaves the dimensions at their zero default, and a
    # scan of no repetitions still holds one
    return Signature(nsamples=int(mrs.nsamples), nviews=int(mrs.nviews),
                     nsliceviews=int(mrs.nsliceviews), nslices=int(mrs.nslices),
                     nechoes=int(mrs.nechoes), nrepetitions=max(int(mrs.nrepetitions), 1),
                     naverages=int(mrs.naverages), datatype=int(mrs.datatype))


def describe_signature(signature: Signature) -> str:
    """The acquisition matrix as a line, so a report or error can show what a group holds"""
    return ", ".join(f"{name}={value}" for name, value in signature._asdict().items())


def scan_id_of(path: str) -> Optional[int]:
    """
    The scanner's scan id for a file, taken from the filename prefix (24804_000_0.MRD) or failing
    that from a numeric parent directory. Ids are assigned in acquisition order, which is what makes
    them worth recovering
    Args:
        - path: filesystem path or tar member name
    Returns:
        - the id, or None for a file that carries neither
    """
    posix = _posix(path)
    match = re.match(r"^(\d+)_", posix.name)
    if match:
        return int(match.group(1))
    parent = posix.parent.name
    return int(parent) if parent.isdigit() else None


def natural_key(path: str) -> tuple:
    """
    Sort key that orders embedded numbers by value, so scan id 999 precedes 1000. A file's position
    in its group becomes its repetition index, so this ordering is load bearing, and plain lexical
    order gets it wrong the first time a scan id gains a digit
    """
    return tuple(int(token) if token.isdigit() else token.lower()
                 for token in re.split(r"(\d+)", str(path)))


def sanitize(component: str) -> str:
    """Make a path component safe to use in a filename"""
    return re.sub(r"[^A-Za-z0-9._-]", "_", component).strip("_") or "unnamed"


def _posix(path) -> PurePosixPath:
    """
    Read a path as posix regardless of where it came from. Tar members are always posix, and a
    sequence name in the parameter block is a Windows path, so backslashes are normalised
    """
    return PurePosixPath(str(path).replace("\\", "/"))


def _is_junk(path) -> bool:
    """
    macOS writes an AppleDouble ._<name> shadow beside every real file, and tars it alongside them.
    Reading one as a .MRD yields garbage, so anything with a ._ path component is dropped, whether
    it was walked off the filesystem or read out of a tar
    """
    return any(part.startswith("._") for part in _posix(path).parts)


def experiment_dir_for(path, experiment_root: str = "", resolve: bool = True) -> Tuple[str, str]:
    """
    The directory identifying the experiment a file belongs to.

    The caller names the experiment rather than the path being read for it: -f is pointed at one
    experiment folder and a tar holds one, so experiment_root is the answer and no directory name
    has to be recognised on the way. This is why there is no list of modality directory names -
    whatever sits between the experiment and the scan directory, and whether anything sits there at
    all, stops mattering.

    Only a lone file arrives with no root to take, and then the scan directory the scanner always
    writes is walked past to reach the folder holding it
    Args:
        - path: filesystem path or tar member name
        - experiment_root: the experiment the caller named, '' to take it from the path
        - resolve: make the path absolute first, so a file named on its own resolves against the
          filesystem. Off for tar members, which have no filesystem to resolve against
    Returns:
        - (experiment_dir, meas_id). meas_id is '' when there is nothing above the scan directory to
          name the experiment after, which leaves the caller to supply one
    """
    if experiment_root:
        return experiment_root, _posix(experiment_root).name
    posix = _posix(os.path.abspath(str(path)) if resolve else path)
    parts = list(posix.parent.parts)
    if not parts:
        return "", ""
    index = len(parts) - 1
    # a scan directory is one the scanner named with its numeric id, e.g. 24804. parts[0] is the
    # anchor ('/') on an absolute path and is never an experiment, hence index > 0
    while index > 0 and parts[index].isdigit():
        index -= 1
    experiment_dir = str(PurePosixPath(*parts[:index + 1]))
    meas_id = parts[index]
    if meas_id in ("/", "") or meas_id.isdigit():
        meas_id = ""
    return experiment_dir, meas_id


@dataclass
class ScanGroup:
    """
    The files converting to one MRD stream, and where that stream goes.

    Paths and scalars only, never parsed files: the probe that decided the grouping is the only
    time these files are read before conversion streams them back one at a time
    """
    meas_id: str
    sequence_name: str
    experiment_dir: str
    output_dir: str
    # real acquisition data, in acquisition order. One path per repetition in case 1, a single path
    # already holding every repetition in case 2, which is what len() > 1 tells apart
    rawdata_file_list: List[str] = field(default_factory=list)
    # averaged prescans, converting into the same stream after the data above rather than a file of
    # their own. The header is never built from one of these
    prescan_file_list: List[str] = field(default_factory=list)
    # repetitions the acquisition holds, counted across rawdata_file_list alone, and the only
    # repetition count a group carries. Recorded so conversion can write the header without having
    # read every file first. A prescan is a calibration beside the acquisition rather than a
    # repetition of it, so it is never counted here however many prescan files arrived
    nrepetitions: int = 0
    signature: Optional[Signature] = None   # the rawdata acquisition matrix, for reporting
    # the shape the whole group's rawdata comes to once its files are concatenated: the shared
    # acquisition matrix's axes carrying nrepetitions above rather than one file's. Recorded here so
    # a caller can size what a group holds without probing its files itself, which is the one number
    # the two cases (a series of files, or one file already holding the axis) do not read the same
    # way. () for a group whose signature was never filled in
    rawdata_shape: tuple = ()
    output_name: str = ""

    @property
    def reference_path(self) -> str:
        """
        The file convert_group_to_mrd builds the header from: real data when there is any, else the
        prescan, since a group always holds at least one file of one kind or the other
        """
        return (self.rawdata_file_list or self.prescan_file_list)[0]

    @property
    def is_combined(self) -> bool:
        """Whether several data files were combined into this group, rather than one converting alone"""
        return len(self.rawdata_file_list) > 1

    @property
    def is_prescan(self) -> bool:
        """
        Whether every file in this group is an averaged prescan, rather than the usual mix of real
        acquisition data plus the prescan that calibrates it
        """
        return bool(self.prescan_file_list) and not self.rawdata_file_list

    @property
    def output_path(self) -> str:
        return os.path.join(self.output_dir, self.output_name)


def collect_mrd_paths(root) -> List[str]:
    """
    Every .MRD file under a directory, recursively
    Args:
        - root: directory to walk, or a single .MRD file
    Returns:
        - paths in natural order. Ordering is settled again once prescan files are sorted to the end
          of the group, but walking in a stable order keeps the output reproducible
    """
    root = Path(root)
    if root.is_file():
        return [str(root)]
    paths = [str(path) for path in root.rglob("*")
             if path.is_file() and path.suffix == MRD_SUFFIX and not _is_junk(path)]
    return sorted(paths, key=natural_key)


def group_experiment(paths: Sequence[str],
                     probe: Callable[[str], MRSdata],
                     root: str = "",
                     experiment_root: str = "",
                     fallback_meas_id: str = "",
                     resolve: bool = True) -> Optional[ScanGroup]:
    """
    Sort the .MRD files of one experiment into the single stream they convert to. The one place the
    rule lives.

    One input is one experiment, so this returns one group rather than a list of them: there is no
    key to sort files under and no bucket for a file to fall into, only the two lists of the group -
    real acquisition data, and the averaged prescan beside it. What the files are is still read off
    the acquisition matrix, and it is still read off a header-only probe, so nothing here parses a
    file it does not have to.

    Two things the grouping key used to absorb quietly are therefore said out loud instead: a second
    acquisition matrix among the real data raises, and a second sequence name anywhere in the
    experiment is reported. Neither has been seen inside one real experiment folder, so both read as
    the input having named more than one
    Args:
        - paths: candidate .MRD paths or tar member names
        - probe: reads one path far enough to give sequence_name, naverages and the dimensions
        - root: directory the caller pointed at. Output stays inside it even when the experiment
          resolves above it, so pointing at a scan directory never writes somewhere unexpected
        - experiment_root: the experiment every path belongs to, see experiment_dir_for. Only a
          lone file leaves this empty and has its experiment read off the path
        - fallback_meas_id: used when the path gives no experiment name, e.g. a tar of a bare scan
          directory
        - meas_id_override: names the experiment outright, ahead of the folder it was read from
        - resolve: see experiment_dir_for
    Returns:
        - the group, its files in acquisition order and its output name assigned, or None when there
          was no .MRD file to convert
    """
    if not paths:
        return None
    root_dir = str(_posix(os.path.abspath(str(root)))) if root and resolve else str(_posix(root))
    # probed in acquisition order, since a file's position in the rawdata list becomes its
    # repetition index. (path, signature, sequence name) is everything the group is built from
    entries: List[Tuple[str, tuple, str]] = []
    for path in sorted(paths, key=natural_key):
        mrs = probe(path)
        entries.append((str(path), signature_of(mrs), mrs.sequence_name))
    # MR Solutions collapses averages into a single acquisition, so an averaged scan arrives as
    # naverages>1 at one repetition, where the acquisition it calibrates arrives as naverages=1
    # repeated - across scan directories, or on the nrepetitions axis of one file
    rawdata = [entry for entry in entries
               if not (entry[1].naverages > 1 and entry[1].nrepetitions == 1)]
    prescan = [entry for entry in entries
               if entry[1].naverages > 1 and entry[1].nrepetitions == 1]
    # real data first so write_header builds the header from an acquisition rather than a prescan,
    # unless the experiment is nothing but prescan
    reference_path, reference_signature, reference_sequence = (rawdata or prescan)[0]
    # every file belongs to the one experiment the caller named, so this is read once rather than
    # per file. The caller naming it outright wins over the folder it was read from, which is what
    # names a file that arrived without an experiment folder around it
    experiment_dir, path_meas_id = experiment_dir_for(reference_path, experiment_root,
                                                     resolve=resolve)
    meas_id = path_meas_id or fallback_meas_id

    # the rawdata files concatenate on the repetition axis, so they have to be the same acquisition;
    # a second matrix among them is misfiled data rather than a longer series. The prescans are held
    # to nothing, since none of them is concatenated onto anything
    matrices = {signature for _, signature, _ in rawdata}
    if len(matrices) > 1:
        described = "; ".join(sorted(describe_signature(matrix) for matrix in matrices))
        raise ValueError(f"{meas_id} holds {len(matrices)} distinct acquisition matrices of real "
                         f"data, which cannot be repetitions of one scan: {described}. Point -f or "
                         f"-t at one experiment, or convert the odd scan on its own with -i")
    # the sequence name is the ppl the scanner ran, so one experiment reads the same one on every
    # file, prescan included - it is not what tells data and calibration apart. It no longer decides
    # anything, so a file carrying another one now joins this stream instead of starting its own,
    # and the header can only record the one
    sequences = {sequence_name for _, _, sequence_name in entries}
    if len(sequences) > 1:
        print(f"WARNING {meas_id} holds {len(sequences)} sequence names "
              f"({', '.join(sorted(sequences))}), converting them into one stream recorded as "
              f"{reference_sequence!r}", file=sys.stderr)

    # output stays inside the directory the caller pointed at. Pointing -f at a scan directory
    # resolves the experiment above it, which is right for the header and wrong for where to write
    inside_root = (not root_dir or experiment_dir == root_dir
                   or experiment_dir.startswith(root_dir.rstrip("/") + "/"))

    group = ScanGroup(
        meas_id=meas_id,
        sequence_name=reference_sequence,
        experiment_dir=experiment_dir,
        output_dir=experiment_dir if inside_root else root_dir,
        rawdata_file_list=[path for path, _, _ in rawdata],
        prescan_file_list=[path for path, _, _ in prescan],
        # summing over the rawdata files resolves both arrival shapes without asking which one this
        # is: one file per repetition contributes 1 each, so the count is how many files there are,
        # and a file already carrying the axis contributes all of them at once. The prescans are
        # summed over by nothing, since none of them is a repetition of the acquisition
        nrepetitions=sum(signature.nrepetitions for _, signature, _ in rawdata),
        signature=reference_signature)
    # the shape the group's rawdata comes to once its files are concatenated: the shared acquisition
    # matrix's axes, carrying the group's repetition count on the last one rather than one file's,
    # since that is the only axis the files concatenate on
    group.rawdata_shape = (reference_signature.nsamples, reference_signature.nviews,
                           reference_signature.nsliceviews, reference_signature.nslices,
                           reference_signature.nechoes, group.nrepetitions)
    # one group per input, so a name cannot collide with a sibling's and nothing has to be broken
    group.output_name = base_name(group) + OUTPUT_SUFFIX
    return group


def base_name(group: ScanGroup) -> str:
    """
    What to call a group's stream.

    A group is named for what it holds, and every name carries the sequence. An averaged prescan says
    so ahead of anything else, since what it is matters more than how many files it arrived in, and
    the acquisition beside it runs the same ppl so the sequence name alone would not tell them apart.
    Files that combined say so, since no one scan id names the set. A file converting on its own adds
    its scan id, which is unique within an experiment. A file with no scan id to add, which is how
    spectral data arrives, is the sequence on its own
    Args:
        - group: a grouped scan with at least one file
    Returns:
        - the name without its suffix
    """
    name = f"{sanitize(group.meas_id)}_{sanitize(group.sequence_name)}"
    if group.is_prescan:
        return f"{name}_prescan"
    if group.is_combined:
        return f"{name}_combined"
    scan_id = scan_id_of(group.reference_path)
    return f"{name}_{scan_id}" if scan_id is not None else name


def organize_folder(root)-> Optional[ScanGroup]:
    """
    Group every .MRD file under one experiment folder into the one stream it converts to.

    The folder is the experiment and names the stream written out of it, so the directories under it
    are read for nothing but their scan ids. Pointing this at a folder of several experiments makes
    them one experiment: whatever the acquisition matrix check lets through concatenates into a
    single series
    Args:
        - root: the experiment folder to walk, or a single .MRD file
    Returns:
        - the group, or None when the folder held no .MRD file
    """
    def probe(path: str) -> MRSdata:
        mrs = MRSdata()
        mrs.probe_from_file(path)
        return mrs

    paths = collect_mrd_paths(root)
    # a lone file names no experiment, so its own path is read for one and output goes beside it
    root_path = Path(root)
    is_file = root_path.is_file()
    clamp_root = str(root_path.parent if is_file else root_path)
    experiment_root = "" if is_file else str(_posix(os.path.abspath(clamp_root)))
    group = group_experiment(paths, probe, root=clamp_root, experiment_root=experiment_root)
    print(f"Grouped {len(paths)} files into {describe_group(group)}", file=sys.stderr)
    return group


def read_scan_tar(stream: BinaryIO) -> Tuple[str, int, List[Tuple[str, bytes]]]:
    """
    Drain a tar stream holding one experiment directory and return its contents in memory, as the
    (member_name, file_bytes) pairs organize_members groups.

    Everything is buffered because none of it can be acted on incrementally: the MR Solutions
    format appends its ASCII parameter block after the data at EOF, the .SPR sidecar may follow
    the .MRD members in tar order and streaming mode cannot rewind to it, and the MRD header needs
    the repetition count across every file in the experiment before the first acquisition can be
    written.

    Args:
        - stream: readable binary stream positioned at the start of a tar archive. May be
                  non-seekable (a FIFO, a socket, sys.stdin.buffer)
    Returns:
        - meas_id: the experiment directory name, i.e. the single top level path component shared by
                   the .MRD members. Empty if the archive has no single root directory, in which
                   case the caller must supply one
        - basefreq: base frequency in Hz from the .SPR sidecar, or 0 if the archive carries none
        - members: (member_name, file_bytes) for each .MRD member, sorted by member name
    """
    mrd_members: List[Tuple[str, bytes]] = []
    basefreq = 0
    roots = set()

    # 'r|*' is tarfile's stream mode: it reads strictly forward and never seeks, unlike 'r'/'r:*'
    # which probe the file and fail on a FIFO with "OSError: [Errno 29] Illegal seek".
    with tarfile.open(fileobj=stream, mode="r|*") as tar:
        for member in tar:
            if not member.isfile() or _is_junk(member.name):
                continue
            payload = tar.extractfile(member)   # valid only until the next iteration in stream mode
            if payload is None:
                continue
            if member.name.endswith(SPR_SUFFIX):
                freq = MRSdata.parse_spr(payload.read())
                if freq:                        # keep an earlier hit if this SPR has no FREQ entry
                    basefreq = freq
            elif member.name.endswith(MRD_SUFFIX):
                mrd_members.append((member.name, payload.read()))
                parts = _posix(member.name).parts
                if len(parts) > 1:
                    roots.add(parts[0])

    if not mrd_members:
        raise ValueError(f"input tar contains no *{MRD_SUFFIX} members")

    # Member order in a tar is filesystem order and is not guaranteed. The position of a file in
    # this list becomes its acquisition repetition index, so sort it explicitly. Scan filenames are
    # zero padded (12345_000_0.MRD, 12345_001_0.MRD) so lexical order is acquisition order.
    mrd_members.sort(key=lambda item: item[0])

    meas_id = roots.pop() if len(roots) == 1 else ""
    return meas_id, basefreq, mrd_members


def organize_members(members: Sequence[Tuple[str, bytes]],
                     fallback_meas_id: str = "")-> Optional[ScanGroup]:
    """
    Group the .MRD members of a tar, whose names are paths but whose contents are already in memory.

    A tar holds one experiment folder, so its root is that experiment and names the stream out of
    it, exactly as the folder given to organize_folder does. A tar sharing no single root directory
    still converts as one stream, named by fallback_meas_id, since a Tyger job has one output buffer
    to write into
    Args:
        - members: (member_name, file_bytes) as read_scan_tar returns them
        - fallback_meas_id: measurement id for members that carry no experiment directory, normally
          the tar's root directory
    Returns:
        - the group, or None when the tar held no .MRD member. Its files are member names, to be
          looked up in members
    """
    payloads = {name: payload for name, payload in members}

    def probe(name: str) -> MRSdata:
        mrs = MRSdata()
        mrs.probe_from_buffer(payloads[name])
        return mrs

    # a tar carries one experiment folder, so nothing in the member paths has to be recognised: the
    # one top level directory they share names the experiment, whatever the directories between are
    # called. '' when they share no single one, and then fallback_meas_id names the scan instead
    roots = {_posix(name).parts[0] for name in payloads if len(_posix(name).parts) > 1}
    group = group_experiment(list(payloads), probe,
                             experiment_root=next(iter(roots)) if len(roots) == 1 else "",
                             fallback_meas_id=fallback_meas_id, resolve=False)
    print(f"Grouped {len(payloads)} members into {describe_group(group)}", file=sys.stderr)
    return group


def describe_group(group: Optional[ScanGroup]) -> str:
    """One line naming what an input resolved to, for the line each entry point logs"""
    if group is None:
        return "nothing to convert"
    return (f"{group.output_name}: {len(group.rawdata_file_list)} data, "
            f"{len(group.prescan_file_list)} prescan")


def report(group: Optional[ScanGroup]) -> None:
    """
    Print what one experiment resolved to, converting nothing
    Args:
        - group: as returned by organize_folder or organize_members
    """
    if group is None:
        print("No .MRD file to convert", file=sys.stderr)
        return
    all_files = group.rawdata_file_list + group.prescan_file_list
    print(f"\n{group.output_name}")
    print(f"  meas_id   {group.meas_id}")
    print(f"  sequence  {group.sequence_name}")
    print(f"  files     {len(group.rawdata_file_list)} data, "
          f"{len(group.prescan_file_list)} prescan")
    print(f"  matrix    {describe_signature(group.signature)}")
    print(f"  shape     {group.rawdata_shape}")
    print(f"  reps      {group.nrepetitions}")
    print(f"  first     {group.reference_path}")
    if len(all_files) > 1:
        print(f"  last      {all_files[-1]}")
    print(f"  output    {group.output_path}")


def main() -> int:
    """
    Print the stream one experiment folder would convert to, without converting anything
    - -f/--folder: the experiment folder to walk
    """
    parser = argparse.ArgumentParser(
        description="Show which .MRD files of one experiment convert into one stream")
    parser.add_argument("-f", "--folder", type=Path, required=True,
                        help="directory containing MRS data files")
    args = parser.parse_args()
    if not args.folder.exists():
        raise SystemExit(f"{args.folder} does not exist")
    report(organize_folder(args.folder))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
