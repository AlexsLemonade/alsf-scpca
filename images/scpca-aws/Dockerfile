FROM continuumio/miniconda3:4.9.2
LABEL maintainer="ccdl@alexslemonade.org"

# Install apt-getable packages needed for nextflow and aws installation
RUN apt-get update && apt-get install -y curl procps zip

# Install boto3 using conda
RUN conda install boto3 

# Install awscli2
RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64-2.1.34.zip" -o "awscliv2.zip" \
    && unzip awscliv2.zip \
    && ./aws/install

WORKDIR /home
