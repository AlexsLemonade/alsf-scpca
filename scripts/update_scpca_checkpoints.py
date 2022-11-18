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


def process_scrna(run, consts, overwrite):
    """
    Process a scrna run, moving internal checkpoint files from source to destination locations
    and adding a scpca-meta.json file.
    """
    print(f"Processing {run.scpca_run_id}")

    # remove any slashes that might be at ends of bucket or prefix
    bucket = consts.bucket.strip('/')
    prefix = consts.prefix.strip('/')

    origin_prefix = f"{prefix}/{consts.source_dir}/rad/{run.scpca_library_id}/{run.scpca_run_id}-rna"
    dest_prefix = f"{prefix}/{consts.dest_dir}/rad/{run.scpca_library_id}/{run.scpca_run_id}-rna"
    metadata_key = f"{dest_prefix}/scpca-meta.json"

    publish_uri = f"s3://{bucket}/{prefix}/{consts.dest_dir}/rad/{run.scpca_library_id}"
    rad_uri = f"{publish_uri}/{run.scpca_run_id}-rna"


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
        "rad_publish_dir": publish_uri,
        "rad_dir": rad_uri,
        "barcode_file": run.barcode_file
    }

    # print(json.dumps(metadata, indent=2))

    ## copy files from source to destination
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
            rad_uri
        ]
        subprocess.run(sync_command)
        ### write json object
        s3_bucket.put_object(
            Key = metadata_key,
            Body = json.dumps(metadata, indent=2)
        )

def process_bulk(run, consts, overwrite):
    """
    Process a bulk run, moving internal checkpoint files from source to destination locations
    and adding a scpca-meta.json file.
    """
    print(f"Processing {run.scpca_run_id}")

    # remove any slashes that might be at ends of bucket or prefix
    bucket = consts.bucket.strip('/')
    prefix = consts.prefix.strip('/')

    origin_prefix = f"{prefix}/{consts.source_dir}/salmon/{run.scpca_library_id}"
    dest_prefix = f"{prefix}/{consts.dest_dir}/salmon/{run.scpca_library_id}"
    metadata_key = f"{dest_prefix}/scpca-meta.json"

    publish_uri = f"s3://{bucket}/{prefix}/{consts.dest_dir}/salmon"
    bulk_uri = f"{publish_uri}/{run.scpca_library_id}"

def process_spatial(run, consts, overwrite):
    """
    Process a spatial run, moving internal checkpoint files from source to destination locations
    and adding a scpca-meta.json file.
    """
    print(f"Processing {run.scpca_run_id}")

    # remove any slashes that might be at ends of bucket or prefix
    bucket = consts.bucket.strip('/')
    prefix = consts.prefix.strip('/')

    origin_prefix = f"{prefix}/{consts.source_dir}/spaceranger/{run.scpca_library_id}/{run.scpca_run_id}-spatial"
    dest_prefix = f"{prefix}/{consts.dest_dir}/spaceranger/{run.scpca_library_id}/{run.scpca_run_id}-spatial"
    metadata_key = f"{dest_prefix}/scpca-meta.json"

    publish_uri = f"s3://{bucket}/{prefix}/{consts.dest_dir}/spaceranger/{run.scpca_library_id}"
    spatial_uri = f"{publish_uri}/{run.scpca_run_id}-spatial"

def process_demux(run, consts, overwrite):
    """
    Process a spatial run, moving internal checkpoint files from source to destination locations
    Since the demux workflow has always included scpca-meta.json, we do not need to create it!
    """
    print(f"Processing {run.scpca_run_id}")

    # remove any slashes that might be at ends of bucket or prefix
    bucket = consts.bucket.strip('/')
    prefix = consts.prefix.strip('/')

    origin_prefix = f"{prefix}/{consts.source_dir}/vireo/{run.scpca_library_id}-vireo"
    dest_prefix = f"{prefix}/{consts.dest_dir}/vireo/{run.scpca_library_id}-vireo"

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
            f"s3://{bucket}/{dest_prefix}"
        ]
        subprocess.run(sync_command)



### Main
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
        '--source_dir',
        default = 'internal',
        help = 'current subdirectory for scpca checkpoint files'
    )
    parser.add_argument(
        '--dest_dir',
        default = 'checkpoints',
        help = 'destination subdirectory for scpca checkpoint files'
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

    ### scRNA processing
    # filter to scRNAseq runs
    sc_techs = ["10Xv2", "10Xv2_5prime", "10Xv3", "10Xv3.1"]
    scRNA_df = library_df.loc[library_df['technology'].isin(sc_techs)]

    # process scRNA runs
    scRNA_df = scRNA_df.iloc[0:2] # limit to 2 for testing
    scRNA_df.apply(process_scrna, axis = 1,
                   consts = args
                   overwrite = args.overwrite)

    ### bulk processing
    bulk_techs = ['single_end', 'paired_end']
    bulk_df = library_df.loc[library_df['technology'].isin(bulk_techs)]
    bulk_df.apply(process_bulk, axis = 1,
                  consts = args
                  overwrite = args.overwrite)

    ### spatial processing
    spatial_techs = ['visium']
    spatial_df = library_df.loc[library_df['technology'].isin(spatial_techs)]
    spatial_df.apply(process_spatial, axis = 1,
                     consts = args
                     overwrite = args.overwrite)

    ### demux processing
    demux_techs = ['cellhash_10Xv2', 'cellhash_10Xv3', 'cellhash_10Xv3.1']
    demux_df = library_df.loc[library_df['technology'].isin(demux)]
    demux_df.apply(process_demux, axis = 1,
                   consts = args
                   overwrite = args.overwrite)

if __name__ == '__main__':
    main()
