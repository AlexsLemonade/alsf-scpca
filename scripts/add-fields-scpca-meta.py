#!/usr/bin/env python3

"""

The purpose of this script is to update the the 'scpca-meta.json' checkpoint files of already processed libraries to be up to date with changes in `scpca-nf`.
This includes addition of the following missing fields in the mapping checkpoints:
- 'ref_mito'
- 'ref_fasta_index'
- 'assay_ontology_term_id'
- 'submitter_cell_types_file'
- 'ref_assembly'
- 'star_index'

For all runs present in the `--library_file`, this script will check for an existing `scpca-meta.json` file in the provided `--checkpoints_prefix` on S3.
If the file is unavailable, the run will be skipped.
If the file exists, the JSON is loaded, and the missing fields are added if they are not already present.
NOTE: This script only updates the `scpca-meta.json` files for mapping results.
For cell type metadata changes, see `add-celltype-fields-scpca-meta.py`.

To run this script for modifying the `scpca-meta.json` files from runs that have already been processed for production do:

python add-fields-scpca-meta.py --checkpoints_prefix "scpca-prod/checkpoints"

"""

import argparse
import json
import boto3
import pandas

# Parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument(
    "--library_file",
    default="s3://ccdl-scpca-data/sample_info/scpca-library-metadata.tsv",
    help="path or URI to library data file TSV",
)
parser.add_argument(
    "--library_ids",
    help="comma separated list of library IDs to update checkpoints",
)
parser.add_argument(
    "--bucket",
    default="nextflow-ccdl-results",
    help="S3 bucket where results files are located",
)
parser.add_argument(
    "--checkpoints_prefix",
    default="scpca/processed/checkpoints",
    help="directory containing checkpoint files to modify",
)
args = parser.parse_args()

# Read in library file
library_df = pandas.read_csv(args.library_file, sep="\t", keep_default_na=False)

# if library ids are provided, filter to only include those
# otherwise keep the full library df
if args.library_ids:
    # get a list of library_ids
    library_ids = args.library_ids.split(",")

    # check that library ids are present in library metadata
    if not all(id in library_df["scpca_library_id"].tolist() for id in library_ids):
        raise ValueError(f"All {library_ids} not found in {args.library_file}")

    # filter library df to only include specified library ids
    library_df = library_df[library_df["scpca_library_id"].isin(library_ids)]

# remove any extra / at the end
bucket = args.bucket.strip("/")
checkpoints_prefix = args.checkpoints_prefix.strip("/")

# define techs for setting up checkpionts directories for different modalities
sc_techs = ["10Xv2", "10Xv2_5prime", "10Xv3", "10Xv3.1"]
bulk_techs = ["single_end", "paired_end"]
spatial_techs = ["visium"]
demux_techs = ["cellhash_10Xv2", "cellhash_10Xv3", "cellhash_10Xv3.1"]

# go through every run id and modify scpca-meta.json file if present
for run in library_df.itertuples():
    print(f"Processing {run.scpca_run_id}")

    # depending on the technology for that run, define the directory where scpca-meta.json is located
    if run.technology in sc_techs:
        checkpoint_folder = (
            f"{checkpoints_prefix}/rad/{run.scpca_library_id}/{run.scpca_run_id}-rna"
        )
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
    s3 = boto3.resource("s3")
    s3_bucket = s3.Bucket(bucket)

    # read scpca-meta.json, if present
    try:
        results_obj = (
            s3_bucket.Object(meta_json_key).get()["Body"].read().decode("utf-8")
        )
        results_meta = json.loads(results_obj)
    except s3.meta.client.exceptions.NoSuchKey:
        print(f"No scpca-meta.json file for {run.scpca_run_id}")
        continue

    # create a list of new fields to check for
    new_fields = {
        "mito_file": "s3://scpca-references/homo_sapiens/ensembl-104/annotation/Homo_sapiens.GRCh38.104.mitogenes.txt",
        "ref_fasta_index": "homo_sapiens/ensembl-104/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai",
        "star_index": "s3://scpca-references/homo_sapiens/ensembl-104/star_index/Homo_sapiens.GRCh38.104.star_idx",
        "assay_ontology_term_id": run.assay_ontology_term_id,
        "ref_assembly": run.sample_reference,
    }

    # check if any of the new fields are already present
    # if they are all present, make sure that submitter_cell_types_file is up to date
    if (
        all(key in results_meta for key in new_fields)
        and results_meta["submitter_cell_types_file"] == run.submitter_cell_types_file
    ):
        print(
            f"All fields are present, no updates to scpca-meta.json for {run.scpca_run_id}"
        )
    else:
        # update missing fields
        for key, default in new_fields.items():
            results_meta.setdefault(key, default)

        # update submitter cell types file if it's missing or if the value doesn't match the metadata file
        if (
            results_meta.get("submitter_cell_types_file")
            != run.submitter_cell_types_file
        ):
            results_meta["submitter_cell_types_file"] = run.submitter_cell_types_file

    # copy updated json file
    s3_bucket.put_object(Key=meta_json_key, Body=json.dumps(results_meta, indent=2))
