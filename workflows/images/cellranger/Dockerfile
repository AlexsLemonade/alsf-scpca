FROM ubuntu:20.04
LABEL maintainer="ccdl@alexslemonade.org"

WORKDIR /opt
COPY cellranger-4.0.0.tar.gz cellranger-4.0.0.tar.gz
RUN tar xzf cellranger-4.0.0.tar.gz && rm cellranger-4.0.0.tar.gz
ENV PATH="/opt/cellranger-4.0.0:${PATH}"
WORKDIR /home
