import os
from pathlib import Path
import sys
import argparse

from MRSreader import MRSdata

mrs = MRSdata()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Read MRS file')
    parser.add_argument('-i', '--input', type=str, required=False, help='Input file, defaults to stdin')
    args = parser.parse_args()
    mrs.mread3d(args.input)
    print(mrs.navg)