FROM centos:7
LABEL maintainer="ccdl@alexslemonade.org"

# Install bcl2fastq
WORKDIR /tmp
COPY bcl2fastq2-v2.20.0.422-Linux-x86_64.rpm bcl2fastq2.rpm
RUN rpm -i bcl2fastq2.rpm 

# Install cellranger
WORKDIR /opt
COPY cellranger-9.0.1.tar.gz cellranger-9.0.1.tar.gz
RUN tar xzf cellranger-9.0.1.tar.gz && rm cellranger-9.0.1.tar.gz
ENV PATH="/opt/cellranger-9.0.1:${PATH}"

WORKDIR /home
