# copy python scripts 
# use conda environment to activate environment
FROM python:3.12-slim
WORKDIR /app
LABEL org.opencontainers.image.source="https://github.com/MEDCAP/mrs_to_mrd"

# install dependencies
COPY requirements.txt ./
RUN pip install --user --no-cache-dir -r requirements.txt

# pip install mrd-fork library first (less likely to change)
COPY mrd-fork/ ./mrd-fork/
RUN pip install --user ./mrd-fork/python

# copy source codes last (most likely to change)
COPY lorn.py ./
COPY MRSreader.py ./
COPY mrdplot.py ./
COPY readmrd.py ./
COPY mrd2recon.py ./

# Make script executable as the default command
# Override CMD to run different scripts
# -u flag runs Python in unbuffered mode so stderr appears immediately
ENTRYPOINT ["python"]
CMD ["mrd2recon.py"]