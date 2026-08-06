"""
Tyger docker container entrypoint: convert one tarred EPSI scan directory into an MRD2 stream.

Reads a tar archive of scan folder from --input and writes the MRD2 binary stream to --output. Both
default to the Tyger buffer pipes, so the codespec needs no arguments:
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import List

from MRSreader import MRSdata
from MRStomrd2 import convert_group_to_mrd
from mrs_tar import read_scan_tar


def convert_epsi_tar(input_path: str, output_path: str, meas_id_override: str = "") -> None:
    """
    Convert a tar of one EPSI scan directory to an MRD2 stream.
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
                         "(tar cf - -C /data cirrhrat_0_1) or pass --meas-id")
    if not basefreq:
        print("No .SPR sidecar in the tar, falling back to the default base frequency",
              file=sys.stderr)

    mrs_list: List[MRSdata] = []
    for name, filebytes in members:
        mrs = MRSdata()
        mrs.parse(filebytes, basefreq)
        # an EPSI scan folder can contain a stray fid file, drop it rather than failing the run
        if "epsi" not in mrs.pplfile:
            print(f"Skipping non-EPSI {name} (ppl={mrs.pplfile})", file=sys.stderr)
            continue
        mrs_list.append(mrs)

    if not mrs_list:
        raise ValueError(f"no EPSI files among the {len(members)} .MRD members of the input tar")
    print(f"Converting {len(mrs_list)} EPSI files as {meas_id}", file=sys.stderr)

    # each EPSI file holds one repetition, so convert_group_to_mrd's default rep_count is correct
    with open(output_path, "wb") as output_stream:
        convert_group_to_mrd(mrs_list, meas_id, output_stream)


def main() -> int:
    parser = argparse.ArgumentParser(description="Convert a tar of single EPSI directory to an MRD2 stream")
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

    convert_epsi_tar(args.input, args.output, args.meas_id)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
