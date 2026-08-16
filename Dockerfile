# Container for the Tyger conversion entrypoint: read a tar of one MR Solutions experiment folder
# from the input buffer, write one MRD2 stream to the output buffer.
#
# One image serves both acquisition families. MRStomrd2.py reads the sequence name out of each .MRD
# and converts EPSI or spectral accordingly, so nothing outside has to know which it has:
#   docker build -t ghcr.io/medcap/mrs-convert:latest .
#
# One experiment is one stream, which is what a job's single output buffer can carry. Files in the
# experiment that could not concatenate onto its acquisition ride in the same stream, flagged
# IS_NOISE_MEASUREMENT.
#
# Local run against plain files (no Tyger involved):
#   docker run --rm -v /path/to/data:/data ghcr.io/medcap/mrs-convert:latest \
#       --tar /data/cirrhrat_43_1.tar --output /data/cirrhrat_43_1_epsigre.mrd2
#
# The same image also converts a mounted folder, writing the stream beside the data:
#   docker run --rm -v /path/to/data:/data ghcr.io/medcap/mrs-convert:latest \
#       --folder /data/cirrhrat_43_1
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

# mrs_organize decides which files are one acquisition, so conversion does not run without it
COPY MRSreader.py mrs_organize.py mrs_tar.py MRStomrd2.py ./

ENTRYPOINT ["python", "MRStomrd2.py"]
