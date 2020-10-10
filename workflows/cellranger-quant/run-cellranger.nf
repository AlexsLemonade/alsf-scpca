#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// basic parameters
params.ref_dir = 's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-100'
params.index_dir = 'cellranger_index'
params.index_name = 'cdna'

params.run_metafile = 's3://ccdl-scpca-data/sample_info/scpca-library-metadata.tsv'
params.run_ids = "SCPCR000001,SCPCR000002" //comma separated list to be parsed into a list
params.outdir = 's3://nextflow-ccdl-results/scpca/cellranger-quant'

// build full paths
params.index_path = "${params.ref_dir}/${params.index_dir}/${params.index_name}"

process cellranger{
  container '589864003899.dkr.ecr.us-east-1.amazonaws.com/scpca-cellranger:4.0.0'
  publishDir "${params.outdir}", mode: 'copy'
  label 'cpus_8'
  input:
    tuple val(id), path(read1), path(read2)
    path index
  output:
    path run_dir
  script:
    run_dir = "${id}-${index}"
    """
    cellranger
    """
}

workflow{
  run_ids = params.run_ids?.tokenize(',') ?: []
  ch_reads = Channel.fromPath(params.run_metafile)
    .splitCsv(header: true, sep: '\t')
    .filter{it.scpca_run_id in run_ids} // use only the rows in the sample list
    // create tuple of [sample_id, [Read1 files], [Read2 files]]
    .map{row -> tuple(row.scpca_run_id,
                      file("s3://${row.s3_prefix}/*_R1_*.fastq.gz"),
                      file("s3://${row.s3_prefix}/*_R2_*.fastq.gz"),
                      )}
  // run cellranger
  cellranger(ch_reads, params.index_path)
}
