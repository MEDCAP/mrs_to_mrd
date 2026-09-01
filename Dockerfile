# Containers for the Tyger conversion entrypoints: read a tar of one MR Solutions scan directory
# from the input buffer, write an MRD2 stream to the output buffer.
#
# Two images from one file, because EPSI and spectral differ in how repetitions are encoded:
#   docker build --target epsi     -t ghcr.io/medcap/mrs-convert-epsi:latest .
#   docker build --target spectral -t ghcr.io/medcap/mrs-convert-spectral:latest .
#
# Local run against plain files (no Tyger involved):
#   docker run --rm -v /path/to/data:/data ghcr.io/medcap/mrs-convert-epsi:latest \
#       --input /data/scan.tar --output /data/raw.mrd2
#
# CPU only: conversion is byte shuffling plus a reshape, there is no GPU work here.
FROM python:3.12-slim AS base

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

COPY MRSreader.py MRStomrd2.py mrs_organize.py epsi_window.py ./

FROM base AS epsi
COPY tyger_convert_epsi.py ./
ENTRYPOINT ["python", "tyger_convert_epsi.py"]

FROM base AS spectral
COPY tyger_convert_spectral.py ./
ENTRYPOINT ["python", "tyger_convert_spectral.py"]
