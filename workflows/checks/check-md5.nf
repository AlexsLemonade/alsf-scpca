#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// run parameters
params.run_metafile = 's3://ccdl-scpca-data/sample_info/scpca-library-metadata.tsv'
params.outdir = 's3://nextflow-ccdl-results/scpca-checks/'

// comma separated list of run ids
// For all samples, use "All"
params.run_ids = "SCPCR000001,SCPCR000002"


process check_md5{
  container 'ubuntu:20.04'
  label 'bigdisk'
  publishDir "${params.outdir}"
  input:
    tuple val(id), val(md5_file), path(files)
  output:
    path outfile
  script:
    outfile = "${id}-md5check.txt"
    """
    md5sum -c ${md5_file} > ${outfile}
    cat ${outfile}
    """
}
workflow{
  run_ids = params.run_ids?.tokenize(',') ?: []
  run_all = run_ids[0] == "All"
  ch_runs = Channel.fromPath(params.run_metafile)
    .splitCsv(header: true, sep: '\t')
    .filter{run_all || (it.scpca_run_id in run_ids)}
    .map{row -> tuple(row.scpca_run_id,
                      row.md5_file,
                      file("s3://${row.s3_prefix}/*")
                      )}
  check_md5(ch_runs)
}
