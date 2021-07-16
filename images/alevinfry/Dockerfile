FROM continuumio/miniconda3:4.9.2
LABEL maintainer="ccdl@alexslemonade.org"

# Install apt-getable packages needed for nextflow
RUN apt-get update && apt-get install -y procps

# Install alevin-fry using bioconda
RUN conda install -c conda-forge -c bioconda alevin-fry=0.4.0

WORKDIR /home
