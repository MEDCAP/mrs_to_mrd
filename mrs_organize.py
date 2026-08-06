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

How to combine:
Combining is the only thing decided here, and it is decided structurally: EPSI files in one
experiment folder combine when their arrays can concatenate on the repetition axis, which means
every dimension of the acquisition agrees

Edge case:
If there is a file that does not match the dimension, it will be converted on its own
e.g. a calibration scan recorded at another matrix size, or an acquisition
that already holds its repetitions on the nrepetitions axis 

experiment_name is not read out of the path, it is what the caller pointed at: -f names one
experiment folder and a tar holds one, so its root is the experiment. Whatever sits between the
experiment and the scan directory, and whether anything sits there at all, never has to be
recognised - which is the point, since that level is named for the sequence, the nucleus or the
contrast and no list of those names stays complete.

Only -i names a lone file with no folder to take, and then the path is walked up past the scan id
the scanner always writes. If nothing above it names the experiment, use the tar filename

Each group is named for what it holds:
    <experiment>_<sequence>_combined.mrd2   files that combined into one acquisition
    <experiment>_<sequence>_<scan_id>.mrd2  a file converting on its own out of a scan directory
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path, PurePosixPath
from typing import Callable, Dict, List, Optional, Sequence, Tuple

from MRSreader import MRSdata

MRD_SUFFIX = ".MRD"
OUTPUT_SUFFIX = ".mrd2"

# substrings that identify the sequence family. Spectral sequence names take several forms across
# scanner generations and protocols: 1puls_extrf_KIC, 1pulsch_phase_clin, one_pulse_clin_gating
EPSI_MARKERS = ("epsi",)
SPECTRAL_MARKERS = ("1pul", "one_pulse", "fid")

# families that record one repetition per scan directory, so sibling directories holding the same
# acquisition matrix are one acquisition split across files
UNIFIED_FAMILIES = frozenset({"epsi"})

# the dimensions that have to agree for two files to be repetitions of one acquisition, which is the
# whole of what decides grouping. Ordered as the raw data is, so a reported difference reads in the
# same order as the shape it describes
SIGNATURE_FIELDS = ("nsamples", "nviews", "nsliceviews", "nslices", "nechoes", "nrepetitions",
                    "naverages", "datatype")


def sequence_family(sequence_name: str) -> str:
    """
    Which conversion family a sequence belongs to
    Args:
        - sequence_name: as read from the SEQUENCE or PPL record, e.g. 'epsigre43_FB_13C'
    Returns:
        - 'epsi', 'spectral', or '' when the name matches neither
    """
    name = (sequence_name or "").lower()
    if any(marker in name for marker in EPSI_MARKERS):
        return "epsi"
    if any(marker in name for marker in SPECTRAL_MARKERS):
        return "spectral"
    return ""


def is_scan_id(name: str) -> bool:
    """A directory the scanner named with its numeric scan id, e.g. 24804"""
    return name.isdigit()


def tar_root(names: Sequence[str]) -> str:
    """
    The experiment a tar holds, which is the one top level directory its members share.

    A tar carries one experiment folder, so nothing in the member paths has to be recognised: the
    root names the experiment and everything below it is that experiment's, whatever the directories
    between are called
    Args:
        - names: tar member names, which are always relative posix paths
    Returns:
        - the root directory name, '' when the members share no single one and the caller has to
          supply a name instead
    """
    roots = {_posix(name).parts[0] for name in names if len(_posix(name).parts) > 1}
    return next(iter(roots)) if len(roots) == 1 else ""


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
    # parts[0] is the anchor ('/') on an absolute path and is never an experiment, hence index > 0
    while index > 0 and is_scan_id(parts[index]):
        index -= 1
    experiment_dir = str(PurePosixPath(*parts[:index + 1]))
    meas_id = parts[index]
    if meas_id in ("/", "") or is_scan_id(meas_id):
        meas_id = ""
    return experiment_dir, meas_id


@dataclass(frozen=True)
class ScanFile:
    """One .MRD file, with everything grouping and checking need and nothing else"""
    path: str
    scan_dir: str
    experiment_dir: str
    meas_id: str
    sequence_name: str
    family: str
    nsamples: int
    nviews: int
    nsliceviews: int
    nslices: int
    nechoes: int
    nrepetitions: int
    naverages: int
    datatype: int
    scan_id: Optional[int]

    @property
    def signature(self) -> tuple:
        """
        Everything that has to agree for two files to be repetitions of one acquisition. Two files
        whose signatures differ cannot concatenate on the repetition axis, so they are separate
        acquisitions however they are filed - a calibration scan at 2176x12 beside a 1792x8 series,
        or a file that already carries 25 repetitions of its own
        """
        return tuple(getattr(self, name) for name in SIGNATURE_FIELDS)

    def describe_signature(self) -> str:
        """
        The acquisition matrix as a line, so a report can show what a group holds. Written with the
        same field names differs_from uses, so the two read against each other
        """
        return ", ".join(f"{name}={getattr(self, name)}" for name in SIGNATURE_FIELDS)

    def differs_from(self, other: "ScanFile") -> str:
        """
        Which dimensions kept this file out of the group beside it, as 'name this vs that' terms.
        Empty when the two agree, which is when they would have combined
        """
        return ", ".join(f"{name} {getattr(self, name)} vs {getattr(other, name)}"
                         for name in SIGNATURE_FIELDS
                         if getattr(self, name) != getattr(other, name))


@dataclass
class ScanGroup:
    """The files converting to one MRD stream, and where that stream goes"""
    meas_id: str
    sequence_name: str
    family: str
    experiment_dir: str
    output_dir: str
    scan_files: List[ScanFile] = field(default_factory=list)
    output_name: str = ""

    @property
    def files(self) -> List[str]:
        """The paths, in acquisition order"""
        return [scan_file.path for scan_file in self.scan_files]

    @property
    def is_combined(self) -> bool:
        """Whether several files were combined into this group, rather than one converting alone"""
        return len(self.scan_files) > 1

    @property
    def nrepetitions(self) -> int:
        """Repetitions the stream will hold, which is what the group is for"""
        return sum(scan_file.nrepetitions for scan_file in self.scan_files)

    @property
    def output_path(self) -> str:
        return os.path.join(self.output_dir, self.output_name)


def collect_mrd_paths(root) -> List[str]:
    """
    Every .MRD file under a directory, recursively
    Args:
        - root: directory to walk, or a single .MRD file
    Returns:
        - paths in natural order. Ordering is settled again per group once families are known, but
          walking in a stable order keeps the output reproducible
    """
    root = Path(root)
    if root.is_file():
        return [str(root)]
    paths = [str(path) for path in root.rglob("*")
             if path.is_file() and path.suffix == MRD_SUFFIX and not _is_junk(path)]
    return sorted(paths, key=natural_key)


def _is_junk(path) -> bool:
    """
    macOS writes an AppleDouble ._<name> shadow beside every real file, in tars and on disk alike.
    Reading one as a .MRD yields garbage
    """
    return any(part.startswith("._") for part in _posix(path).parts)


def describe(path: str, mrs: MRSdata, experiment_root: str = "", resolve: bool = True) -> ScanFile:
    """
    Build the grouping record for one probed file
    Args:
        - path: filesystem path or tar member name
        - mrs: the file, probed or fully parsed
        - experiment_root: see experiment_dir_for
        - resolve: see experiment_dir_for
    Returns:
        - ScanFile
    """
    experiment_dir, meas_id = experiment_dir_for(path, experiment_root, resolve=resolve)
    return ScanFile(path=str(path),
                    scan_dir=str(_posix(path).parent),
                    experiment_dir=experiment_dir,
                    meas_id=meas_id,
                    sequence_name=mrs.sequence_name,
                    family=sequence_family(mrs.sequence_name),
                    nsamples=int(mrs.nsamples),
                    nviews=int(mrs.nviews),
                    nsliceviews=int(mrs.nsliceviews),
                    nslices=int(mrs.nslices),
                    nechoes=int(mrs.nechoes),
                    # a probe that failed on a malformed header leaves the dimensions at their zero
                    # default, and a scan of no repetitions still holds one
                    nrepetitions=max(int(mrs.nrepetitions), 1),
                    naverages=int(mrs.naverages),
                    datatype=int(mrs.datatype),
                    scan_id=scan_id_of(path))


def check_group(group: ScanGroup) -> List[str]:
    """
    Look for signs that a group is not the intact scan it appears to be. Files are sometimes misfiled
    on upload, and the directory structure that grouping trusts cannot report that by itself.

    Every check is a warning. Each has a benign explanation often enough that refusing to convert
    would be wrong, and a converted stream with a warning beside it is more useful than neither
    Args:
        - group: a grouped scan, files in acquisition order
    Returns:
        - warning lines, empty when the group looks intact
    """
    warnings: List[str] = []

    # the scan ids of one acquisition run consecutively, so a hole is a file that never arrived. A
    # calibration scan sits outside the run, and is already in a group of its own by then
    acquired = sorted(f.scan_id for f in group.scan_files if f.scan_id is not None)
    if len(acquired) > 1:
        missing = [scan_id for previous, scan_id in zip(acquired, acquired[1:])
                   for scan_id in range(previous + 1, scan_id)]
        if missing:
            warnings.append(f"{group.meas_id} is missing scan id{'s' if len(missing) > 1 else ''} "
                            f"{_summarize(missing)} between {acquired[0]} and {acquired[-1]}, so "
                            f"{len(missing)} file(s) did not arrive")

    sequence_names = {f.sequence_name for f in group.scan_files}
    if len(sequence_names) > 1:
        warnings.append(f"{group.meas_id} mixes sequences {sorted(sequence_names)}")
    return warnings


def check_groups(groups: Sequence[ScanGroup]) -> List[str]:
    """
    Look for signs of misfiling that only show up between groups rather than inside one.

    A file that agrees with nothing around it becomes a group of one, which is the normal shape of a
    calibration scan and is not reported. Two combined series in one experiment directory is the
    reportable case: it is legitimate if a rig renamed its ppl mid session or reconfigured the
    matrix, and is a misfiled upload the rest of the time, so it is reported rather than merged
    Args:
        - groups: every group from one run
    Returns:
        - warning lines, empty when nothing looks misfiled
    """
    warnings: List[str] = []
    by_experiment: Dict[tuple, List[ScanGroup]] = {}
    for group in groups:
        if group.is_combined:
            by_experiment.setdefault((group.experiment_dir, group.meas_id), []).append(group)
    for (_, meas_id), siblings in by_experiment.items():
        if len(siblings) > 1:
            detail = ", ".join(f"{g.sequence_name} ({len(g.scan_files)} files)"
                               for g in sorted(siblings, key=lambda g: g.sequence_name))
            warnings.append(f"{meas_id} holds {len(siblings)} combined series, converted "
                            f"separately: {detail}")
    return warnings


def _summarize(numbers: Sequence[int], limit: int = 6) -> str:
    """Render a list of ids without letting a large hole fill the log"""
    shown = ", ".join(str(number) for number in numbers[:limit])
    return shown if len(numbers) <= limit else f"{shown}, … (+{len(numbers) - limit} more)"


def group_files(paths: Sequence[str],
                probe: Callable[[str], MRSdata],
                root: str = "",
                experiment_root: str = "",
                fallback_meas_id: str = "",
                resolve: bool = True) -> List[ScanGroup]:
    """
    Sort .MRD files into the streams they convert to. The one place the rule lives
    Args:
        - paths: candidate .MRD paths or tar member names
        - probe: reads one path far enough to give sequence_name, naverages and the timestamp
        - root: directory the caller pointed at. Output stays inside it even when the experiment
          resolves above it, so pointing at a scan directory never writes somewhere unexpected
        - experiment_root: the experiment every path belongs to, see experiment_dir_for. Only a
          lone file leaves this empty and has its experiment read off the path
        - fallback_meas_id: used when the path gives no experiment name, e.g. a tar of a bare scan
          directory
        - resolve: see experiment_dir_for
    Returns:
        - groups, each with its files in acquisition order and its output name assigned
    """
    root_dir = str(_posix(os.path.abspath(str(root)))) if root and resolve else str(_posix(root))
    # EPSI buckets by acquisition matrix, so scan directories recording one repetition each collect
    # together and anything acquired differently stays out. Everything else converts a file at a
    # time: its repetitions are already inside it, and two scans can share a sequence name
    buckets: Dict[tuple, List[ScanFile]] = {}
    alone: List[ScanFile] = []
    for path in sorted(paths, key=natural_key):
        scan_file = describe(path, probe(path), experiment_root, resolve=resolve)
        if not scan_file.family:
            print(f"{path}: sequence {scan_file.sequence_name!r} is neither EPSI nor spectral, "
                  f"converting it on its own", file=sys.stderr)
        if scan_file.family in UNIFIED_FAMILIES:
            key = (scan_file.experiment_dir, scan_file.family, scan_file.sequence_name,
                   scan_file.signature)
            buckets.setdefault(key, []).append(scan_file)
        else:
            alone.append(scan_file)

    members = _combine(buckets) + [[scan_file] for scan_file in alone]
    groups = [ScanGroup(meas_id=scan_files[0].meas_id or fallback_meas_id,
                        sequence_name=scan_files[0].sequence_name,
                        family=scan_files[0].family,
                        experiment_dir=scan_files[0].experiment_dir,
                        output_dir=_clamp(scan_files[0].experiment_dir, root_dir),
                        scan_files=scan_files)
              for scan_files in members]
    # tree order, which is what a reader walking the folder alongside the report expects
    groups.sort(key=lambda group: natural_key(group.scan_files[0].path))
    assign_output_names(groups)
    return groups


def _combine(buckets: Dict[tuple, List[ScanFile]]) -> List[List[ScanFile]]:
    """
    Decide which buckets convert as one stream and which convert a file at a time.

    An experiment folder records one repetition series, so within one experiment and sequence the
    largest bucket is that series and combines. The rest were acquired differently and only happen
    to be filed beside it - a calibration scan, a repeat at another matrix - so each of their files
    converts on its own and keeps a scan id to be named by. Without that, two calibration scans run
    alike would read as one acquisition of two repetitions: cirrhrat_43_1's 24792 and 24793 share a
    matrix and consecutive ids, and are two scans rather than one.

    Buckets tied for the largest all combine, since nothing here can say which of them is the
    series. That is the misfiled upload check_groups reports
    Args:
        - buckets: files sharing an experiment, sequence and acquisition matrix
    Returns:
        - the file lists each becoming one group
    """
    largest: Dict[tuple, int] = {}
    for key, scan_files in buckets.items():
        experiment = key[:-1]                       # the key without the acquisition matrix
        largest[experiment] = max(largest.get(experiment, 0), len(scan_files))
    members: List[List[ScanFile]] = []
    for key, scan_files in buckets.items():
        if len(scan_files) > 1 and len(scan_files) == largest[key[:-1]]:
            members.append(scan_files)
        else:
            members.extend([scan_file] for scan_file in scan_files)
    return members


def _clamp(experiment_dir: str, root_dir: str) -> str:
    """
    Keep output inside the directory the caller pointed at. Pointing -f at a scan directory resolves
    the experiment above it, which is right for the header and wrong for where to write
    """
    if not root_dir:
        return experiment_dir
    if experiment_dir == root_dir or experiment_dir.startswith(root_dir.rstrip("/") + "/"):
        return experiment_dir
    return root_dir


def base_name(group: ScanGroup) -> str:
    """
    What to call a group's stream, before ties are broken.

    A group is named for what it holds, and every name carries the sequence. Files that combined say
    so, since no one scan id names the set. A file converting on its own adds its scan id, which is
    unique within an experiment and is what distinguishes it from the series beside it - the sequence
    name alone would not, because a calibration scan runs the same ppl as the acquisition it
    calibrates. A file with no scan id to add, which is how spectral data arrives, is the sequence on
    its own
    Args:
        - group: a grouped scan with at least one file
    Returns:
        - the name without its suffix
    """
    name = f"{sanitize(group.meas_id)}_{sanitize(group.sequence_name)}"
    if group.is_combined:
        return f"{name}_combined"
    scan_id = group.scan_files[0].scan_id
    return f"{name}_{scan_id}" if scan_id is not None else name


def assign_output_names(groups: Sequence[ScanGroup]) -> None:
    """
    Name each group's stream and break ties.

    Two groups can still propose one name: two combined series in an experiment differ by their
    acquisition matrix rather than by anything in the name. When that happens every colliding group
    takes its scan directory into its name, rather than only the later ones, so a name does not
    depend on the order the tree was walked
    Args:
        - groups: assigned in place
    """
    proposed: Dict[tuple, List[ScanGroup]] = {}
    for group in groups:
        proposed.setdefault((group.output_dir, base_name(group)), []).append(group)
    for (_, name), colliding in proposed.items():
        if len(colliding) == 1:
            colliding[0].output_name = name + OUTPUT_SUFFIX
            continue
        for group in colliding:
            tag = sanitize(_posix(group.scan_files[0].path).parent.name)
            group.output_name = f"{name}_{tag}{OUTPUT_SUFFIX}"
        print(f"{len(colliding)} scans in {colliding[0].output_dir} share the name {name}, "
              f"distinguishing them by scan directory", file=sys.stderr)


def organize_folder(root, quiet: bool = True) -> List[ScanGroup]:
    """
    Group every .MRD file under one experiment folder.

    The folder is the experiment and names every stream written out of it, so the directories under
    it are read for nothing but their scan ids. Pointing this at a folder of several experiments
    puts them all in one experiment, and two series acquired at the same matrix would combine
    Args:
        - root: the experiment folder to walk, or a single .MRD file
        - quiet: suppress MRSreader's per record reporting, which repeats once per file
    Returns:
        - groups, ordered and named
    """
    def probe(path: str) -> MRSdata:
        mrs = MRSdata()
        mrs.probe_from_file(path, quiet=quiet)
        return mrs

    paths = collect_mrd_paths(root)
    # a lone file names no experiment, so its own path is read for one and output goes beside it
    root_path = Path(root)
    is_file = root_path.is_file()
    clamp_root = str(root_path.parent if is_file else root_path)
    experiment_root = "" if is_file else str(_posix(os.path.abspath(clamp_root)))
    groups = group_files(paths, probe, root=clamp_root, experiment_root=experiment_root)
    print(f"Grouped {len(paths)} files into {len(groups)} scans", file=sys.stderr)
    return groups


def organize_members(members: Sequence[Tuple[str, bytes]],
                     fallback_meas_id: str = "",
                     quiet: bool = True) -> List[ScanGroup]:
    """
    Group the .MRD members of a tar, whose names are paths but whose contents are already in memory.

    A tar holds one experiment folder, so its root is that experiment and names every stream out of
    it, exactly as the folder given to organize_folder does
    Args:
        - members: (member_name, file_bytes) as read_scan_tar returns them
        - fallback_meas_id: measurement id for members that carry no experiment directory, normally
          the tar's root directory
        - quiet: suppress MRSreader's per record reporting
    Returns:
        - groups, ordered and named. Their files are member names, to be looked up in members
    """
    payloads = {name: payload for name, payload in members}

    def probe(name: str) -> MRSdata:
        mrs = MRSdata()
        mrs.probe_from_buffer(payloads[name], quiet=quiet)
        return mrs

    groups = group_files(list(payloads), probe, experiment_root=tar_root(list(payloads)),
                         fallback_meas_id=fallback_meas_id, resolve=False)
    print(f"Grouped {len(payloads)} members into {len(groups)} scans", file=sys.stderr)
    return groups


def reference_group(groups: Sequence[ScanGroup]) -> Optional[ScanGroup]:
    """
    The group the others in an experiment are a departure from: the acquisition the experiment was
    for. That is the combined series when there is one, and otherwise the group holding the most
    repetitions, which is the same reasoning _combine uses to pick the series in the first place
    Args:
        - groups: the groups of one experiment
    Returns:
        - the reference, or None when there is only one group and nothing to compare it against
    """
    if len(groups) < 2:
        return None
    combined = [group for group in groups if group.is_combined]
    return max(combined or groups, key=lambda group: (len(group.scan_files), group.nrepetitions))


def report(groups: Sequence[ScanGroup]) -> int:
    """
    Print the grouping and its integrity warnings.

    A file that converts on its own is the part of a report worth explaining, because the folder it
    sits in says it belongs with the rest and only its dimensions say otherwise. Each one is printed
    against the acquisition it did not join, naming the dimensions that differ, so the reason it was
    left out is on the page rather than something to go and work out
    Args:
        - groups: as returned by organize_folder or organize_members
    Returns:
        - number of warnings raised across all groups
    """
    total = 0
    for warning in check_groups(groups):
        total += 1
        print(f"WARNING {warning}")
    by_experiment: Dict[tuple, List[ScanGroup]] = {}
    for group in groups:
        by_experiment.setdefault((group.experiment_dir, group.meas_id), []).append(group)
    references = {id(group): reference_group(siblings)
                  for siblings in by_experiment.values() for group in siblings}
    for group in groups:
        first = group.scan_files[0]
        print(f"\n{group.output_name}")
        print(f"  meas_id   {group.meas_id}")
        print(f"  sequence  {group.sequence_name} ({group.family})")
        print(f"  matrix    {first.describe_signature()}")
        reference = references.get(id(group))
        if reference is not None and reference is not group:
            difference = first.differs_from(reference.scan_files[0])
            print(f"  separate  did not combine into {reference.output_name}: "
                  f"{difference or 'same dimensions, so it is a second series'}")
        print(f"  first     {first.path}")
        if len(group.scan_files) > 1:
            print(f"  last      {group.scan_files[-1].path}")
        print(f"  output    {group.output_path}")
        for warning in check_group(group):
            total += 1
            print(f"  WARNING   {warning}")
    print(f"\n{len(groups)} scans, {total} warning(s)", file=sys.stderr)
    return total


def main() -> int:
    """
    Print how a folder would be grouped, without converting anything
    - -f/--folder: directory to walk
    """
    parser = argparse.ArgumentParser(
        description="Show which .MRD files would be grouped into one experiment")
    parser.add_argument("-f", "--folder", type=Path, required=True,
                        help="directory containing MRS data files")
    args = parser.parse_args()
    if not args.folder.exists():
        raise SystemExit(f"{args.folder} does not exist")
    report(organize_folder(args.folder))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
