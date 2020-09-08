FROM continuumio/miniconda3:4.8.2
LABEL maintainer="ccdl@alexslemonade.org"


#setup bioconda
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge


# Install kallisto-bustools with conda
#########################################
RUN conda install kallisto bustools