"""
Tyger entrypoint: convert one tarred spectral (1pul/fid) scan directory into an MRD2 stream.

Reads a tar archive of one scan from --input and writes the MRD2 binary stream to --output. Both
default to the Tyger buffer pipes, so the codespec needs no arguments:

    tar cf - -C /data PYR_HepG22.mrs | tyger run exec -f tyger_deploy/convert_spectral_codespec.yml > raw.mrd2

Both paths may be named FIFOs. Nothing here seeks, so a pipe works exactly like a file, which keeps
the script testable locally against plain tar files. EPSI scans go through tyger_convert_epsi.py
instead.

The difference that matters between the two: an EPSI scan spreads its repetitions across files, one
per file, whereas a spectral scan holds all of them on the nex axis inside a single file and
generate_acquisition numbers them from that axis. A spectral tar therefore carries exactly one .MRD.
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import List

from MRSreader import MRSdata
from MRStomrd2 import convert_group_to_mrd
from mrs_tar import read_scan_tar


def convert_spectral_tar(input_path: str, output_path: str, meas_id_override: str = "") -> None:
    """
    Convert a tar of one spectral scan directory to an MRD2 stream.
    Args:
        - input_path: path to a tar archive, or to a FIFO carrying one
        - output_path: path to write the MRD2 stream to, or a FIFO to write it into
        - meas_id_override: measurement id to record, overriding the tar's root directory name
    Returns:
        - None
    """
    # opening a FIFO for read blocks until the buffer sidecar opens the write end
    with open(input_path, "rb") as input_stream:
        meas_id, basefreq, members = read_scan_tar(input_stream)

    if meas_id_override:
        meas_id = meas_id_override
    if not meas_id:
        raise ValueError("could not derive a measurement id: tar the scan directory itself "
                         "(tar cf - -C /data PYR_HepG22.mrs) or pass --meas-id")
    if not basefreq:
        print("No .SPR sidecar in the tar, falling back to the default base frequency",
              file=sys.stderr)

    mrs_list: List[MRSdata] = []
    for name, filebytes in members:
        mrs = MRSdata()
        mrs.parse(filebytes, basefreq)
        # a spectral scan folder can contain a stray epsi file, drop it rather than failing the run
        if not ("1pul" in mrs.pplfile or "fid" in mrs.pplfile):
            print(f"Skipping non-spectral {name} (ppl={mrs.pplfile})", file=sys.stderr)
            continue
        mrs_list.append(mrs)

    if not mrs_list:
        raise ValueError(f"no spectral files among the {len(members)} .MRD members of the input tar")
    if len(mrs_list) > 1:
        # generate_acquisition numbers spectral repetitions from each file's own nex axis, so a
        # second file would restart at repetition 0 and silently overwrite the first in recon
        raise ValueError(f"expected one spectral .MRD in the tar but found {len(mrs_list)}; "
                         "tar a single scan directory")
    print(f"Converting spectral scan {meas_id} with {mrs_list[0].rawdata.shape[5]} repetitions",
          file=sys.stderr)

    # spectral repetitions come from the nex axis, not from the file count
    rep_count = mrs_list[0].rawdata.shape[5]
    with open(output_path, "wb") as output_stream:
        convert_group_to_mrd(mrs_list, meas_id, output_stream, rep_count)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Convert a tar of one spectral scan directory to an MRD2 stream")
    parser.add_argument("-i", "--input", default=os.environ.get("INPUT_PIPE"),
                        help="tar archive or FIFO to read (default: $INPUT_PIPE)")
    parser.add_argument("-o", "--output", default=os.environ.get("OUTPUT_PIPE"),
                        help="file or FIFO to write the MRD2 stream to (default: $OUTPUT_PIPE)")
    parser.add_argument("--meas-id", default="",
                        help="measurement id, overriding the tar's root directory name")
    args = parser.parse_args()

    if not args.input:
        parser.error("--input is required when $INPUT_PIPE is unset")
    if not args.output:
        parser.error("--output is required when $OUTPUT_PIPE is unset")

    convert_spectral_tar(args.input, args.output, args.meas_id)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
