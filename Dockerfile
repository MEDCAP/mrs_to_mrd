# Containers for the Tyger pipeline stages on this branch.
#
# Three images from one file, one per pipeline stage:
#   docker build --platform linux/amd64 --target convert -t ghcr.io/medcap/mrs-convert:latest .
#   docker build --platform linux/amd64 --target shift   -t ghcr.io/medcap/mrs-shift:latest .
#   docker build --platform linux/amd64 --target recon   -t ghcr.io/medcap/mrs-recon:latest .
#
# --platform linux/amd64 is not optional. The Tyger cluster runs amd64 nodes and a Mac builds
# arm64 by default, which lands as an exec format error inside the job rather than at build time.
#
# Local run against plain files (no Tyger involved):
#   docker run --rm -v /path/to/data:/data ghcr.io/medcap/mrs-convert:latest \
#       -t /data/scan.tar -o /data/raw.mrd2
#   docker run --rm -v /path/to/data:/data ghcr.io/medcap/mrs-shift:latest \
#       -i /data/raw.mrd2 -o /data/corrected.mrd2
#   docker run --rm -v /path/to/data:/data ghcr.io/medcap/mrs-recon:latest \
#       -i /data/corrected.mrd2 -o /data/recon.mrd2 -pyr_s 9.7 -lac_m 21.8
#
# CPU only. Conversion is byte shuffling plus a reshape, and the recon's fits are numpy and scipy.
FROM python:3.12-slim AS base

LABEL org.opencontainers.image.source=https://github.com/MEDCAP/mrs_to_mrd

# MRD_VERSION_STRING is required at build time: the mrd fork's setup.py reads its version from it
ENV PYTHONUNBUFFERED=1 \
    MRD_VERSION_STRING=2.0.0

WORKDIR /app

# requirements.txt is the whole dependency set, the same one used locally. An image-only
# requirements file is what let this Dockerfile go stale against a renamed module last time, so
# there is one file and the images carry matplotlib they never import.
COPY requirements.txt .

# git is only needed to fetch the mrd fork, so install, use, and drop it in one layer
RUN apt-get update \
    && apt-get install -y --no-install-recommends git \
    && pip install --no-cache-dir -r requirements.txt \
    && apt-get purge -y --auto-remove git \
    && rm -rf /var/lib/apt/lists/*

# MRStomrd2 and mrd2shift share a reader and a grouper, and MRStomrd2 imports mrd2shift for its
# peak check, so both entrypoints ship from one layer rather than from two near-identical copies.
FROM base AS pipeline
COPY MRSreader.py MRSorganize.py mrd2shift.py MRStomrd2.py ./

# Stage 1. Tar of one scan directory in, MRD2 stream out.
FROM pipeline AS convert
ENTRYPOINT ["python", "MRStomrd2.py"]

# Stage 2. Converted MRD2 stream in, the same stream with the EPSI echo drift taken out.
FROM pipeline AS shift
ENTRYPOINT ["python", "mrd2shift.py"]

# Stage 3. Corrected MRD2 stream in, the same stream plus the fitted NdArrays out.
FROM base AS recon
COPY lorentzian_fitter.py mrd2recon.py ./
ENTRYPOINT ["python", "mrd2recon.py"]
