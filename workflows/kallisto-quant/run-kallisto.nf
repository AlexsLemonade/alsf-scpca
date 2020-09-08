#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// run parameters
params.ref_dir = 's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-100'
params.index_dir = 'kallisto_index'
params.index_name = 'cdna_k31'
params.annotation_dir = 'annotation'
params.t2g = 'Homo_sapiens.ensembl.100.tx2gene.tsv'
params.mitolist = 'Homo_sapiens.ensembl.100.mitogenes.txt'

params.sample_dir = 's3://ccdl-scpca-data/raw/green_adam'
params.sample_ids = "834,905_3" //comma separated list to be parsed into a list
params.outdir = 's3://nextflow-ccdl-results/scpca/kallisto-quant'

// build full paths
params.index_path = "${params.ref_dir}/${params.index_dir}/${params.index_name}"
params.t2g_path = "${params.ref_dir}/${params.annotation_dir}/${params.t2g}"
params.mito_path = "${params.ref_dir}/${params.annotation_dir}/${params.mitolist}"


process kallisto{
  container 'quay.io/biocontainers/kallisto:0.46.2--h4f7b962_1'
  cpus 8
  // try dynamic memory (28.GB so 2x will fit in r4.2xlarge)
  memory { 28.GB * task.attempt}
  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
  maxRetries 1
  tag "${id}-${index}"
  publishDir "${params.outdir}"
  input:
    tuple val(id), path(read1), path(read2)
    path index
  output:
    path "${id}-${index}_bus"
  script:
    // interleave read1 & read2 files
    reads = [read1, read2].transpose().flatten().join(' ')
    """
    kallisto bus \
      -i ${index} \
      -o ${id}-${index}_bus \
      -x 10xv3 \
      -t ${task.cpus} \
      ${reads}
    """
}

workflow{
  sample_ids = params.sample_ids?.tokenize(',') ?: []
  ch_reads = Channel.fromList(sample_ids)
    // create tuple of [sample_id, [Read1 files], [Read2 files]]
    .map{ id -> tuple("$id",
                      file("${params.sample_dir}/${id}/*_R1_*.fastq.gz"),
                      file("${params.sample_dir}/${id}/*_R2_*.fastq.gz"),
                      )}
  // run Kallisto
  kallisto(ch_reads, params.index_path)
}