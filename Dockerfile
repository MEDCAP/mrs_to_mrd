# Container for the MRS to MRD2 conversion entrypoint: read a tar of one MR Solutions scan
# directory from the input buffer, write an MRD2 stream to the output buffer.
#
# One image serves both acquisition families. MRStomrd2.py reads the sequence name out of each
# .MRD and converts EPSI or spectral accordingly, so nothing outside has to know which it has:
#   docker build -t ghcr.io/medcap/mrs-convert:latest .
#
# Local run against plain files (no Tyger involved):
#   docker run --rm -v /path/to/data:/data ghcr.io/medcap/mrs-convert:latest \
#       --tar /data/scan.tar --output /data/raw.mrd2
#
# The same image also converts a mounted folder, which is how the *_recon.sh scripts drive it:
#   docker run --rm -v /path/to/data:/data ghcr.io/medcap/mrs-convert:latest \
#       --folder /data/cirrhrat_data -u 3
#
# CPU only: conversion is byte shuffling plus a reshape, there is no GPU work here.
FROM python:3.12-slim

LABEL org.opencontainers.image.source=https://github.com/MEDCAP/mrs_to_mrd

# MRD_VERSION_STRING is required at build time: the mrd fork's setup.py reads its version from it
ENV PYTHONUNBUFFERED=1 \
    MRD_VERSION_STRING=2.0.0

WORKDIR /app

COPY requirements-tyger.txt .

# git is only needed to fetch the mrd fork, so install, use, and drop it in one layer
RUN apt-get update \
    && apt-get install -y --no-install-recommends git \
    && pip install --no-cache-dir -r requirements-tyger.txt \
    && apt-get purge -y --auto-remove git \
    && rm -rf /var/lib/apt/lists/*

COPY MRSreader.py MRStomrd2.py mrs_tar.py ./

ENTRYPOINT ["python", "MRStomrd2.py"]
