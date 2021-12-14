FROM centos:7
LABEL maintainer="ccdl@alexslemonade.org"

# Install bcl2fastq
WORKDIR /tmp
COPY bcl2fastq2-v2.20.0.422-Linux-x86_64.rpm bcl2fastq2.rpm
RUN rpm -i bcl2fastq2.rpm 

# Install cellranger
WORKDIR /opt
COPY cellranger-6.1.2.tar.gz cellranger-6.1.2.tar.gz
RUN tar xzf cellranger-6.1.2.tar.gz && rm cellranger-6.1.2.tar.gz
ENV PATH="/opt/cellranger-6.1.2:${PATH}"

WORKDIR /home
