FROM ubuntu:20.04
LABEL maintainer="ccdl@alexslemonade.org"

RUN apt update && apt install -y build-essential libzip-dev

# Install bcl2fastq
WORKDIR /tmp
COPY bcl2fastq2-v2.20.0.422-Source.tar.gz bcl2fastq2.tar.gz
COPY install-bcl2fastq.sh install-bcl2fastq.sh
RUN bash install-bcl2fastq.sh

# Install cellranger
WORKDIR /opt
COPY cellranger-6.0.1.tar.gz cellranger-6.0.1.tar.gz
RUN tar xzf cellranger-6.0.1.tar.gz && rm cellranger-6.0.1.tar.gz
ENV PATH="/opt/cellranger-6.0.1:${PATH}"

WORKDIR /home
