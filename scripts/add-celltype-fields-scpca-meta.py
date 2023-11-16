#!/usr/bin/env python3

"""

The purpose of this script is to update the the 'scpca-meta.json' checkpoint files for cell typing of already processed libraries to be up to date with changes in `scpca-nf`.
This includes addition of the following missing fields in the celltype checkpoints:
- `singler_ref_version`
- `cellassign_ref_version`

For all runs present in the `--library_file`, this script will check the `celltype` folder for an existing `scpca-meta.json` file in the provided `--checkpoints_prefix` on S3.
Both the `scpca-meta.json` for SingleR and CellAssign will be updated.
If the file is unavailable, the run will be skipped.
If the file exists, the JSON is loaded, and the missing fields are added if they are not already present.
To run this script for modifying the cell type `scpca-meta.json` files from runs that have already been processed for production do:

python add-celltype-fields-scpca-meta.py --checkpoints_prefix "scpca-prod/checkpoints"

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
    "--project_celltype_file",
    default="s3://ccdl-scpca-data/sample_info/scpca-project-celltype-metadata.tsv",
    help="path to URI to project cell type data file TSV",
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

# Read in project celltype file
project_df = pandas.read_csv(
    args.project_celltype_file, sep="\t", keep_default_na=False
)


# join library and project cell type metadata
all_metadata_df = library_df.join(
    project_df.set_index("scpca_project_id"), on="scpca_project_id", how="left"
)
# need to filter to only projects that have at least one reference (either SingleR or CellAssign)
all_metadata_df = all_metadata_df[
    all_metadata_df[["singler_ref_file", "cellassign_ref_file"]].notnull().all(1)
]

# remove any extra / at the end
bucket = args.bucket.strip("/")
checkpoints_prefix = args.checkpoints_prefix.strip("/")

# set up S3
s3 = boto3.resource("s3")
s3_bucket = s3.Bucket(bucket)

# go through every run id and modify scpca-meta.json files if present
for run in all_metadata_df.itertuples():
    print(f"Processing {run.scpca_run_id}")

    # build paths to singleR and CellAssign checkpoints
    celltype_checkpoints = {
        "singler_checkpoint": f"{checkpoints_prefix}/celltype/{run.scpca_library_id}/{run.scpca_library_id}_singler/scpca-meta.json",
        "cellassign_checkpoint": f"{checkpoints_prefix}/celltype/{run.scpca_library_id}/{run.scpca_library_id}_cellassign/scpca-meta.json",
    }

    # cell type fields that need to be updated
    new_fields = {
        "singler_model_file": f"s3://scpca-references/celltype/singler_models/{run.singler_ref_file}",
        "cellassign_reference_file": f"s3://scpca-references/celltype/cellassign_references/{run.cellassign_ref_file}",
    }

    # for each checkpoint file, add cellassign and singler file path
    for checkpoint_file in celltype_checkpoints.values():
        # make sure file exists, otherwise just skip
        try:
            checkpoint_obj = (
                s3_bucket.Object(checkpoint_file).get()["Body"].read().decode("utf-8")
            )
            checkpoint_meta = json.loads(checkpoint_obj)
        except s3.meta.client.exceptions.NoSuchKey:
            print(f"No cell type scpca-meta.json file for {run.scpca_run_id}")
            continue

        # only modify if existing fields don't match what's in the project metadata
        if (
            checkpoint_meta["singler_model_file"] == new_fields["singler_model_file"]
            and checkpoint_meta["cellassign_reference_file"]
            == new_fields["cellassign_reference_file"]
        ):
            print(
                f"All fields are present, no updates for {run.scpca_run_id} to {checkpoint_file}"
            )
        else:
            # before modifying, make sure that the entries in the project metadata aren't NA
            if run.singler_ref_file != "NA":
                checkpoint_meta["singler_model_file"] = new_fields["singler_model_file"]
            else:
                checkpoint_meta["singler_model_file"] = "NA"
            if run.cellassign_ref_file != "NA":
                checkpoint_meta["cellassign_reference_file"] = new_fields[
                    "cellassign_reference_file"
                ]
            else:
                checkpoint_meta["cellassign_reference_file"] = "NA"

        # copy updated json file
        s3_bucket.put_object(
            Key=checkpoint_file, Body=json.dumps(checkpoint_meta, indent=2)
        )
