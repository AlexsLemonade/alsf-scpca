FROM rocker/tidyverse:4.0.2
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
##########################
# Install bioconductor packages
RUN Rscript -e "BiocManager::install(c( \
    'AnnotationHub', \
    'DropletUtils', \
    'ensembldb', \
    'fishpond', \
    'scran', \
    'scater', \
    'tximport'), \
    update = FALSE)"

###############
# Other R packages
RUN Rscript -e "install.packages(c( \
    'here', \
    'optparse', \
    'rprojroot') \
    )"

# set final workdir for commands
WORKDIR /home/rstudio
