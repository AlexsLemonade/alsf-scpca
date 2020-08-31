#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// run parameters
params.ref_dir = 's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-100'
params.index_dir = 'salmon_index'
params.index_name = 'cdna_k31'
params.annotation_dir = 'annotation'
params.t2g = 'Homo_sapiens.ensembl.100.tx2gene.tsv'
params.mitolist = 'Homo_sapiens.ensembl.100.mitogenes.txt'
params.sample_dir = 's3://ccdl-scpca-data/raw/green_adam'
params.sample_id = ['834', '905_3']
params.outdir = 's3://nextflow-ccdl-results/scpca/alevin-quant'

// build full paths
params.index_path = "${params.ref_dir}/${params.index_dir}/${params.index_name}"
params.t2g_path = "${params.ref_dir}/${params.annotation_dir}/${params.t2g}"
params.mito_path = "${params.ref_dir}/${params.annotation_dir}/${params.mitolist}"

process alevin{
  container 'quay.io/biocontainers/salmon:1.3.0--hf69c8f4_0'
  cpus 8
  publishDir "${params.outdir}"
  input:
    tuple val(id), path(read1), path(read2)
    path index
    path tx2gene
  output:
    path "${id}_${index}"
  script:
    """
    salmon alevin \
      -l ISR \
      --chromium \
      -1 ${read1} \
      -2 ${read2} \
      -i ${index} \
      --tgMap ${tx2gene} \
      -o ${id}_${index} \
      -p ${task.cpus} \
      --dumpFeatures \
    """
}

workflow{
  ch_reads = Channel.fromList(params.sample_id)
    // create tuple of [sample_id, [Read1 files], [Read2 files]]
    .map{ id -> tuple("$id",
                      file("${params.sample_dir}/${id}/*_R1_*.fastq.gz"),
                      file("${params.sample_dir}/${id}/*_R2_*.fastq.gz"),
                      )}
  // run Alevin
  alevin(ch_reads, params.index_path, params.t2g_path)
}