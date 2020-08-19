#!/bin/bash
set -euo pipefail

s3_base=s3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-100

# Get reference fasta files and sync to S3
wget -N -P fasta -i fasta_ref_urls.txt
aws s3 sync fasta $s3_base/fasta

# get annotation files & sync
wget -N -P annotation -i annotation_ref_urls.txt
aws s3 sync annotation $s3_base/annotation




