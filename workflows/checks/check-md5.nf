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
  // tasks sometimes fail due to disk space but we can't specify that directly
  errorStrategy 'retry' 
  maxRetries = 1
  cpus { 1 + 3 * (task.attempt - 1)  } //use more CPUs to take over a machine if failure

  publishDir "${params.outdir}"
  input:
    tuple val(id), val(md5_file), path(files)
  output:
    path outfile
  script:
    outfile = "${id}-md5check.txt"
    """
    md5sum -c ${md5_file} > ${outfile}
    """
}

process cat_md5{
  // combine md5 results files, removing blank lines
  container 'ubuntu:20.04'
  publishDir "${params.outdir}"
  input:
    path(files)
  output:
    path outfile
  script:
    outfile = "00_all_md5check.txt"
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
                      row.md5_file,
                      file("s3://${row.s3_prefix}/*")
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
