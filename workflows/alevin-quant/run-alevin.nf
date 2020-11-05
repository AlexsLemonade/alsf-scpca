#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// run parameters
params.ref_dir = 's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-100'
params.index_dir = 'salmon_index'
params.index_name = 'cdna_k31'
params.annotation_dir = 'annotation'
params.t2g = 'Homo_sapiens.ensembl.100.tx2gene.tsv'
params.mitolist = 'Homo_sapiens.ensembl.100.mitogenes.txt'

params.run_metafile = 's3://ccdl-scpca-data/sample_info/scpca-library-metadata.tsv'
params.run_ids = "SCPCR000001,SCPCR000002" //comma separated list to be parsed into a list
params.outdir = 's3://nextflow-ccdl-results/scpca/alevin-quant'

// build full paths
params.index_path = "${params.ref_dir}/${params.index_dir}/${params.index_name}"
params.t2g_path = "${params.ref_dir}/${params.annotation_dir}/${params.t2g}"
params.mito_path = "${params.ref_dir}/${params.annotation_dir}/${params.mitolist}"

process alevin{
  container 'quay.io/biocontainers/salmon:1.3.0--hf69c8f4_0'
  label 'cpus_8'
  tag "${id}-${index}"
  publishDir "${params.outdir}"
  input:
    tuple val(id), path(read1), path(read2)
    path index
    path tx2gene
  output:
    path "${id}-${index}"
  script:
    """
    salmon alevin \
      -l ISR \
      --chromium \
      -1 ${read1} \
      -2 ${read2} \
      -i ${index} \
      --tgMap ${tx2gene} \
      -o ${id}-${index} \
      -p ${task.cpus} \
      --dumpFeatures \
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
  // run Alevin
  alevin(ch_reads, params.index_path, params.t2g_path)
}