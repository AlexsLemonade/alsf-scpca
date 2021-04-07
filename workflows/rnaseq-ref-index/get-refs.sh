#!/bin/bash
set -euo pipefail

# update to use most recent genome release
s3_base=s3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-103
# run R in the scpca_r dockerfile
# docker pull ghcr.io/alexslemonade/scpca_r
# rdocker="docker run --mount type=bind,target=/home/rstudio,source=$PWD ghcr.io/alexslemonade/scpca-r:4.0.5"

# Get reference fasta files and sync to S3
wget -N -P fasta -i fasta_ref_urls.txt
# combine cdna & ncrna fasta files
cat fasta/Homo_sapiens.GRCh38.cdna.all.fa.gz \
  fasta/Homo_sapiens.GRCh38.ncrna.fa.gz \
  > fasta/Homo_sapiens.GRCh38.txome.fa.gz

aws s3 sync fasta $s3_base/fasta

# get annotation files & sync
wget -N -P annotation -i annotation_ref_urls.txt

# running without docker for now due to docker error
# $rdocker Rscript make_pre_mrna_fasta.R

Rscript make_pre_mrna_fasta.R

aws s3 sync annotation $s3_base/annotation
