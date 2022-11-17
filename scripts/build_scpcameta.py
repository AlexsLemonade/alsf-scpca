#!/usr/bin/env python3

# script to build scpca-metadata.json files for checkpoint files

import argparse
import json
import subprocess

import boto3
import pandas

### Functions ###
def get_barcode(tech, barcode_dir):
    """Get the barcode file appropriate for a 10X technology"""
    barcode_dir = barcode_dir.rstrip("/") # remove trailing slash if needed
    if "10Xv3" in tech:
        barcode_file = f"{barcode_dir}/3M-february-2018.txt"
    elif "10Xv2" in tech:
        barcode_file = f"{barcode_dir}/737K-august-2016.txt"
    else:
        barcode_file = "NA"

    return barcode_file


def process_scrna(run, consts, overwrite = False):
    print(f"Processing {run.scpca_run_id}")

    # remove any slashes that might be at ends of in bucket or prefix
    bucket = consts.bucket.strip('/')
    prefix = consts.prefix.strip('/')

    origin_prefix = f"{prefix}/internal/rad/{run.scpca_library_id}/{run.scpca_run_id}-rna"
    dest_prefix = f"{prefix}/checkpoints/rad/{run.scpca_library_id}/{run.scpca_run_id}-rna"

    publish_dir = f"s3://{bucket}/{prefix}/checkpoints/rad/{run.scpca_library_id}"
    rad_dir = f"{publish_dir}/{run.scpca_run_id}-rna"
    meta_file = f"{rad_dir}/scpca-meta.json"

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
        "ref_assembly": consts.ref_assembly,
        "ref_fasta": consts.ref_fasta,
        "ref_gtf": consts.ref_gtf,
        "salmon_splici_index": consts.salmon_splici_index,
        "salmon_bulk_index": consts.salmon_bulk_index,
        "t2g_3col_path": consts.t2g_3col_path,
        "t2g_bulk_path": consts.t2g_bulk_path,
        "cellranger_index": consts.cellranger_index,
        "star_index": consts.star_index,
        "scpca_version": consts.scpca_version,
        "nextflow_version": consts.nextflow_version,
        "rad_publish_dir": publish_dir,
        "rad_dir": rad_dir,
        "barcode_file": run.barcode_file
    }

    # print(json.dumps(metadata, indent=2))

    # set up S3
    s3 = boto3.resource('s3')
    s3_bucket = s3.Bucket(bucket)

    # list origin and destination to see if files exist
    origin_objs = list(s3_bucket.objects.filter(Prefix = origin_prefix))
    dest_objs = list(s3_bucket.objects.filter(Prefix = dest_prefix))
    if len(origin_objs) == 0:
        print(f"No files for {run.scpca_run_id}")
    elif len(dest_objs) > 0 and not overwrite:
        print(f"Files present at destination for {run.scpca_run_id} -- skipping")
    else:
        print(f"Copying files for {run.scpca_run_id}")
        ### S3 copying using awscli
        sync_command = [
            "aws", "s3", "sync",
            f"s3://{bucket}/{origin_prefix}",
            rad_dir,
            "--dryrun"
        ]
        subprocess.run(sync_command)




### Main code
def main():
    ### Parse command line arguments
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
    parser.add_argument(
        '--overwrite',
        action="store_true",
        help = "overwrite existing files at destination"
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

    ### Read in library file
    library_df = pandas.read_csv(args.library_file, sep = "\t",  dtype = 'string', keep_default_na=False)

    # add barcode files
    library_df['barcode_file'] = library_df['technology'].apply(get_barcode, barcode_dir = args.barcode_dir)

    # filter to scRNAseq runs
    sc_techs = ["10Xv2", "10Xv2_5prime", "10Xv3", "10Xv3.1"]
    scRNA_df = library_df.loc[library_df['technology'].isin(sc_techs)]

    # process scRNA runs
    scRNA_df = scRNA_df.iloc[0:3]
    scRNA_df.apply(process_scrna, axis = 1, consts = args, overwrite = args.overwrite)



# print json
# print(json.dumps(metadata, indent=2))

if __name__ == '__main__':
    main()
