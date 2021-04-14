FROM rocker/tidyverse:4.0.5
LABEL maintainer="ccdl@alexslemonade.org"

### Install apt-getable packages to start
#########################################
RUN apt-get update && apt-get install -y --no-install-recommends apt-utils dialog

# Add curl, bzip2 and dev libs
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    awscli \
    bzip2 \
    curl \
    libbz2-dev \
    libgdal-dev \
    libglpk-dev \
    liblzma-dev \
    libreadline-dev \
    libudunits2-dev


#### R packages
###############
# this comes first since bioconductor packages are 
# dependent on updated matrixStats
RUN Rscript -e "install.packages(c( \
    'here', \
    'optparse', \
    'matrixStats', \
    'rprojroot') \
    )"

##########################
# Install bioconductor packages
RUN Rscript -e "BiocManager::install(c( \
    'AnnotationHub', \
    'Biostrings', \
    'BSgenome', \
    'DropletUtils', \
    'eisaR', \
    'ensembldb', \
    'fishpond', \
    'GenomicFeatures', \
    'scran', \
    'scater', \
    'tximport'), \
    update = FALSE)"

# set final workdir for commands
WORKDIR /home/rstudio
