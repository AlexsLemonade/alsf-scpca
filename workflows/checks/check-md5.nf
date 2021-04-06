#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// run parameters
params.run_metafile = 's3://ccdl-scpca-data/sample_info/scpca-library-metadata.tsv'
params.outdir = 's3://nextflow-ccdl-results/scpca-checks/'

// comma separated list of run ids
// For all samples, use "All"
params.run_ids = "SCPCR000001,SCPCR000002"


process check_md5{
  container 'ghcr.io/alexslemonade/scpca-aws'
  publishDir "${params.outdir}"
  input:
    tuple val(id), val(prefix), path(md5_file) 
  output:
    path outfile
  script:
    outfile = "${id}-md5check.txt"
    """
    check_md5_s3.py -m ${md5_file} -p ${prefix} > ${outfile}
    """
}

process cat_md5{
  // combine md5 results files, removing blank lines
  container 'ghcr.io/alexslemonade/scpca-aws'
  publishDir "${params.outdir}"
  input:
    path(files)
  output:
    path outfile
  script:
    outfile = "00_all-md5check.txt"
    """
    cat ${files} | sed '/^\\s+\$/d' > ${outfile}
    """
}

workflow{
  run_ids = params.run_ids?.tokenize(',') ?: []
  run_all = run_ids[0] == "All"
  ch_runs = Channel.fromPath(params.run_metafile)
    .splitCsv(header: true, sep: '\t')
    .filter{run_all || (it.scpca_run_id in run_ids)}
    .map{row -> tuple(row.scpca_run_id,
                      row.s3_prefix,
                      file("s3://${row.s3_prefix}/${row.md5_file}")
                      )}
  check_md5(ch_runs)
  cat_md5(check_md5.out.collect())
  // print failures
  cat_md5.out
    .splitText()
    .map{it.toString().trim()} // remove whitespace
    .filter{!(it =~/OK$/)} // print if not OK
    .view()
}
