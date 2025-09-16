FROM python:3.12-slim
WORKDIR /app
LABEL org.opencontainers.image.source=https://github.com/MEDCAP/mrs_to_mrd


# install dependencies
COPY requirements.txt ./
RUN pip install --user --no-cache-dir -r requirements.txt

# copy source codes
COPY lorn.py ./
COPY mrd2recon.py ./
COPY mrd2scrub.py ./
COPY MRSreader.py ./
COPY MRStomrd2.py ./
COPY stream_recon.py ./

# copy MRD library
COPY python ./