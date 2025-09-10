FROM python:3.12-slim
WORKDIR /app

# install dependencies
COPY requirements.txt ./
RUN pip install --user --no-cache-dir -r requirements.txt

# copy source codes
COPY lorn.py ./
COPY mrd2recon.py ./
COPY mrd2scrub.py ./
COPY MRSreader.py ./
COPY MRStomrd2.py ./

# copy MRD library
COPY python ./

# make output directory to save mrd file
RUN mkdir -p /app/output

ENTRYPOINT ["python", "./MRStomrd2.py"]