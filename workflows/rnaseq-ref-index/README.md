# Workflow for indexing RNA-seq references

The scripts in this directory implement a workflow for downloading and indexing genome adn transcriptome references for RNAseq and uploading them to S3.

## Download

All reference files are downloaded from Ensembl for consistency, using the current version as of this writing: v100.
Fasta files to be downloaded are listed in `fasta_ref_urls.txt` and other annotation files can be listed in `annotation_ref_urls.txt`.
These files are then downloaded by running the bash script `get-refs.sh` which will download the selected files to a local machine, then sync to
the following S3 location: `s3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-100`

## Indexing

Currently, only salmon indexing is implemented, without decoys.
The workflow expects to find the reference cdna file at `s3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-100/fasta/Homo_sapiens.GRCh38.ncrna.fa.gz` and will place the salmon index directories at `s3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-100/salmon_index` with subdirectories for each type of index.

To run the index script locally, navigate to this directory and run
```
nextflow -C ../nextflow.config run build-index.nf
```

This will use the provided `nextflow.config` file, which allows for running with docker images (the needed image will be pulled automatically)

To run on AWS Batch, use instead:

```
nextflow -C ../nextflow.config run build-index.nf -profile batch
```