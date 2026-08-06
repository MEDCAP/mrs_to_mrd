"""
Read a tar stream containing one MR Solutions scan directory.

Tyger hands a job its input buffer as a named FIFO, which is strictly sequential: lseek() on it
returns ESPIPE. That rules out the three filesystem lookups the folder-based converter relies on
(the .SPR sidecar found via Path.parent.iterdir(), the directory names meas_id is derived from, and
the grouping of files belonging to one scan). Wrapping the scan directory in a tar puts all three
into a single sequential stream:

    tar cf - -C /data cirrhrat_0_1 | tyger buffer write $input_buffer

One tar is one scan and produces one MRD v2 stream. Grouping across scans stays in MRStomrd2's
folder mode; the caller is responsible for tarring exactly one scan directory.
"""

from __future__ import annotations

import tarfile
from pathlib import PurePosixPath
from typing import BinaryIO, List, Tuple

from MRSreader import MRSdata

MRD_SUFFIX = ".MRD"
SPR_SUFFIX = ".SPR"


def _is_junk(name: str) -> bool:
    """
    macOS tars carry an AppleDouble ._<name> shadow member alongside every real file. Reading one
    as a .MRD yields garbage, so drop anything with a ._ path component.
    """
    return any(part.startswith("._") for part in PurePosixPath(name).parts)


def read_scan_tar(stream: BinaryIO) -> Tuple[str, int, List[Tuple[str, bytes]]]:
    """
    Drain a tar stream holding one scan directory and return its contents in memory.

    Everything is buffered because none of it can be acted on incrementally: the MR Solutions
    format appends its ASCII parameter block after the data at EOF, the .SPR sidecar may follow
    the .MRD members in tar order and streaming mode cannot rewind to it, and the MRD header needs
    the repetition count across every file in the scan before the first acquisition can be written.

    Args:
        - stream: readable binary stream positioned at the start of a tar archive. May be
                  non-seekable (a FIFO, a socket, sys.stdin.buffer)
    Returns:
        - meas_id: the scan directory name, i.e. the single top level path component shared by the
                   .MRD members. Empty if the archive has no single root directory, in which case
                   the caller must supply one
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
                parts = PurePosixPath(member.name).parts
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
