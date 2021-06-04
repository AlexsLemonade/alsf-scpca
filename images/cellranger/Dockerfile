FROM centos:centos8
LABEL maintainer="ccdl@alexslemonade.org"

RUN yum makecache && yum install -y unzip

# Install bcl2fastq
WORKDIR /tmp
COPY bcl2fastq2-v2-20-0-linux-x86-64.zip bcl2fastq2-v2-20-0-linux-x86-64.zip
RUN unzip bcl2fastq2-v2-20-0-linux-x86-64.zip
RUN rpm -i bcl2fastq2-v2.20.0.422-Linux-x86_64.rpm \
&& rm bcl2fastq2-v2.20.0.422-Linux-x86_64.rpm && \
 rm bcl2fastq2-v2-20-0-linux-x86-64.zip
# Install cellranger
WORKDIR /opt
COPY cellranger-6.0.1.tar.gz cellranger-6.0.1.tar.gz
RUN tar xzf cellranger-6.0.1.tar.gz && rm cellranger-6.0.1.tar.gz
ENV PATH="/opt/cellranger-6.0.1:${PATH}"

WORKDIR /home
