#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// run parameters
params.run_metafile = 's3://ccdl-scpca-data/sample_info/scpca-library-metadata.tsv'
params.outdir = 's3://nextflow-ccdl-results/scpca-checks/'

process check_md5{
  container 'ubuntu:20.04'
  publishDir "${params.outdir}"
  input:
    tuple val(id), val(md5_file), path(files)
  output:
    path "${id}-md5check.txt"
  script:
    """
    md5sum -c ${md5_file} > ${id}-md5check.txt
    """
}
workflow{
  ch_runs = Channel.fromPath(params.run_metafile)
    .splitCsv(header: true, sep: '\t')
    .map{row -> tuple(row.scpca_run_id,
                      row.md5_file,
                      file("s3://${row.s3_prefix}/*")
                      )}
  check_md5(ch_runs)
}