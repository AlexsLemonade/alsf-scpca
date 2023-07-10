#!/usr/bin/env python3

import argparse
import json
import boto3
import pandas

# Parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument(
    '--library_file',
    default='s3://ccdl-scpca-data/sample_info/scpca-library-metadata.tsv',
    help='path or URI to library data file TSV'
)
parser.add_argument(
    '--bucket',
    default='nextflow-ccdl-results',
    help='S3 bucket where results files are located'
)
parser.add_argument(
    '--checkpoints_prefix',
    default='scpca/processed/checkpoints',
    help='directory containing checkpoint files to modify'
)
args = parser.parse_args()

# Read in library file
library_df = pandas.read_csv(args.library_file, sep="\t", keep_default_na=False)

# remove any extra / at the end
bucket = args.bucket.strip('/')
checkpoints_prefix = args.checkpoints_prefix.strip('/')

# define techs for setting up checkpionts directories for different modalities
sc_techs = ["10Xv2", "10Xv2_5prime", "10Xv3", "10Xv3.1"]
bulk_techs = ['single_end', 'paired_end']
spatial_techs = ['visium']
demux_techs = ['cellhash_10Xv2', 'cellhash_10Xv3', 'cellhash_10Xv3.1']

# go through every run id and modify scpca-meta.json file if present
for run in library_df.itertuples():

    print(f"Processing {run.scpca_run_id}")

    # depending on the technology for that run, define the directory where scpca-meta.json is located
    if run.technology in sc_techs:
        checkpoint_folder = f"{checkpoints_prefix}/rad/{run.scpca_library_id}/{run.scpca_run_id}-rna"
    elif run.technology in bulk_techs:
        checkpoint_folder = f"{checkpoints_prefix}/salmon/{run.scpca_library_id}"
    elif run.technology in spatial_techs:
        checkpoint_folder = f"{checkpoints_prefix}/spaceranger/{run.scpca_library_id}/{run.scpca_run_id}-spatial"
    elif run.technology in demux_techs:
        checkpoint_folder = f"{checkpoints_prefix}/vireo/{run.scpca_library_id}-vireo"
    else:
        continue

    # define key for scpca-meta.json file
    meta_json_key = f"{checkpoint_folder}/scpca-meta.json"

    # set up S3
    s3 = boto3.resource('s3')
    s3_bucket = s3.Bucket(bucket)

    # read scpca-meta.json, if present
    try:
        results_obj = s3_bucket.Object(meta_json_key).get()['Body'].read().decode('utf-8')
        results_meta = json.loads(results_obj)
    except s3.meta.client.exceptions.NoSuchKey:
        print(f"No scpca-meta.json file for {run.scpca_run_id}")
        continue

    # add mito file and ref fasta index if not already present
    if all(key in results_meta for key in ('mito_file', 'ref_fasta_index')):
        print(f"All keys are present, no updates to scpca-meta.json for {run.scpca_run_id}")
    else:
        results_meta['mito_file'] = "s3://scpca-references/homo_sapiens/ensembl-104/annotation/Homo_sapiens.GRCh38.104.mitogenes.txt"
        results_meta['ref_fasta_index'] = "homo_sapiens/ensembl-104/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai"

    # copy updated json file
    s3_bucket.put_object(
        Key=meta_json_key,
        Body=json.dumps(results_meta, indent=2)
    )
