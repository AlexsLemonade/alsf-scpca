FROM rocker/tidyverse:4.0.2
LABEL maintainer="ccdl@alexslemonade.org"

### Install apt-getable packages to start
#########################################
RUN apt-get update && apt-get install -y --no-install-recommends apt-utils dialog

# Add curl, bzip2 and dev libs
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    bzip2 \
    curl \
    libbz2-dev \
    libgdal-dev \
    liblzma-dev \
    libreadline-dev \
    libudunits2-dev


#### R packages
###############

# Commonly used R packages
RUN install2.r --error --deps TRUE \
    optparse \
    rprojroot

##########################
# Install bioconductor packages
RUN R -e "BiocManager::install(c( \
    'AnnotationHub', \
    'scran', \
    'scater', \
    'tximport'), \
    update = FALSE)"

# set final workdir for commands
WORKDIR /home/rstudio
