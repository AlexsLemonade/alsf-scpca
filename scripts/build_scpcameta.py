#!/usr/bin/env python3

# script to build scpca-metadata.json files for checkpoint files

import argparse
import json

import boto3
import pandas
import smart_open


### set up arguments

parser = argparse.ArgumentParser()
parser.add_argument(
    '--library_file',
    default = 's3://ccdl-scpca-data/sample_info/scpca-library-metadata.tsv',
    help = 'path or URI to library data file TSV'
)
parser.add_argument(
    '--bucket',
    default = 'nextflow-ccdl-results',
    help = 'S3 bucket where results files are located'
)
parser.add_argument(
    '--prefix',
    default = 'scpca/processed',
    help = 'base location of scpca results files'
)
const = parser.add_argument_group('Common Values')
const.add_argument(
    "--ref_assembly",
    default = "Homo_sapiens.GRCh38.104",
    help = "reference assembly"
)
const.add_argument(
    "--ref_fasta",
    default = "s3://scpca-references/homo_sapiens/ensembl-104/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
    help = "path or uri to genome reference fasta"
)
const.add_argument(
    "--ref_gtf",
    default = "s3://scpca-references/homo_sapiens/ensembl-104/annotation/Homo_sapiens.GRCh38.104.gtf.gz",
    help = "path or uri to genome reference gtf"
)
const.add_argument(
    "--salmon_splici_index",
    default = "s3://scpca-references/homo_sapiens/ensembl-104/salmon_index/Homo_sapiens.GRCh38.104.spliced_intron.txome",
    help = "path or uri to splici index"
)
const.add_argument(
    "--salmon_bulk_index",
    default = "s3://scpca-references/homo_sapiens/ensembl-104/salmon_index/Homo_sapiens.GRCh38.104.spliced_cdna.txome",
    help = "path or uri to bulk index"
)
const.add_argument(
    "--t2g_3col_path",
    default = "s3://scpca-references/homo_sapiens/ensembl-104/annotation/Homo_sapiens.GRCh38.104.spliced_intron.tx2gene_3col.tsv",
    help = "path or uri to 3 column t2g file for splici index"
)
const.add_argument(
    "--t2g_bulk_path",
    default = "s3://scpca-references/homo_sapiens/ensembl-104/annotation/Homo_sapiens.GRCh38.104.spliced_cdna.tx2gene.tsv",
    help = "path or uri to bulk t3g file"
)
const.add_argument(
    "--cellranger_index",
    default = "s3://scpca-references/homo_sapiens/ensembl-104/cellranger_index/Homo_sapiens.GRCh38.104_cellranger_full",
    help = "path or uri to Cell Ranger index"
)
const.add_argument(
    "--star_index",
    default = "s3://scpca-references/homo_sapiens/ensembl-104/star_index/Homo_sapiens.GRCh38.104.star_idx",
    help = "path or uri to STAR index"
)
const.add_argument(
    "--nextflow_version",
    default = "22.10.0",
    help = "default nextflow version"
)
const.add_argument(
    "--scpca_version",
    default = "v0.3.3",
    help = "default scpca version"
)
const.add_argument(
    "--barcode_dir",
    default = "s3://scpca-references/barcodes/10X",
    help = "path or uri to directory with barcode files"
)
args = parser.parse_args()

def get_barcode(tech, barcode_dir):
    """Get the barcode file appropriate for a 10X technology"""
    if "10Xv3" in tech:
        barcode_file = f"{barcode_dir}/3M-february-2018.txt"
    elif "10Xv2" in tech:
        barcode_file = f"{barcode_dir}/737K-august-2016.txt"
    else:
        barcode_file = "NA"

    return barcode_file


def process_scrna(run, bucket, prefix):
    print(f"Processing {run.scpca_run_id}")

    publish_dir = f"s3://{bucket}/{prefix}/checkpoints/rad/{run.scpca_library_id}"
    rad_dir = f"s3://{bucket}/{prefix}/checkpoints/rad/{run.scpca_library_id}/{run.scpca_library_id}-rna"

    ### Build metadata object
    metadata = {
        "run_id": run.scpca_run_id,
        "library_id": run.scpca_library_id,
        "sample_id": run.scpca_sample_id ,
        "project_id": run.scpca_project_id ,
        "submitter": run.submitter,
        "technology": run.technology,
        "seq_unit": run.seq_unit,
        "feature_barcode_file": run.feature_barcode_file,
        "feature_barcode_geom": run.feature_barcode_geom,
        "files_directory": run.files_directory,
        "slide_serial_number": run.slide_serial_number,
        "slide_section": run.slide_section,
        "ref_assembly": args.ref_assembly,
        "ref_fasta": args.ref_fasta,
        "ref_gtf": args.ref_gtf,
        "salmon_splici_index": args.salmon_splici_index,
        "salmon_bulk_index": args.salmon_bulk_index,
        "t2g_3col_path": args.t2g_3col_path,
        "t2g_bulk_path": args.t2g_bulk_path,
        "cellranger_index": args.cellranger_index,
        "star_index": args.star_index,
        "scpca_version": args.scpca_version,
        "nextflow_version": args.nextflow_version,
        "rad_publish_dir": publish_dir,
        "rad_dir": rad_dir,
        "barcode_file": run.barcode_file
    }


    print(json.dumps(metadata, indent=2))

### Read in library file
library_df = pandas.read_csv(args.library_file, sep = "\t",  dtype = 'string', keep_default_na=False)
library_df['barcode_file'] = library_df['technology'].apply(get_barcode, barcode_dir = args.barcode_dir)


print(library_df.iloc[0])



# filter to scRNAseq runs
sc_techs = ["10Xv2", "10Xv2_5prime", "10Xv3", "10Xv3.1"]
scRNA_df = library_df.loc[library_df['technology'].isin(sc_techs)]

scRNA_df.apply(
    process_scrna, axis = 1,
    bucket = args.bucket,
    prefix = args.prefix
    )



# print json
# print(json.dumps(metadata, indent=2))
