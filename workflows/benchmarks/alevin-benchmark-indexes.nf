#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// run parameters
params.run_metafile = 's3://ccdl-scpca-data/sample_info/scpca-library-metadata.tsv'
params.run_ids = "SCPCR000001,SCPCR000002" //comma separated list to be parsed into a list
params.outdir = 's3://nextflow-ccdl-results/scpca-benchmark/alevin-quant'

process alevin{
  container 'quay.io/biocontainers/salmon:1.3.0--hf69c8f4_0'
  cpus 8
  // try dynamic memory (28.GB so 2x will fit in r4.2xlarge)
  memory { 28.GB * task.attempt}
  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
  maxRetries 1
  tag "${sample_id}-${index_id}"
  publishDir "${params.outdir}"
  input:
    tuple val(sample_id), path(read1), path(read2), val(index_id), path(index), path(tx2gene)
  output:
    path "${sample_id}-${index_id}"
  script:
    """
    salmon alevin \
      -l ISR \
      --chromium \
      -1 ${read1} \
      -2 ${read2} \
      -i ${index} \
      --tgMap ${tx2gene} \
      -o ${sample_id}-${index_id} \
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
                      file("${row.files_directory}/*_R1_*.fastq.gz"),
                      file("${row.files_directory}/*_R2_*.fastq.gz"),
                      )}
  ch_indexes = Channel.fromList([
    ['cdna_k31_no_sa',
     's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-100/salmon_index/cdna_k31',
     's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-100/annotation/Homo_sapiens.ensembl.100.tx2gene.tsv'],
    ['cdna_k23_no_sa',
     's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-100/salmon_index/cdna_k23',
     's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-100/annotation/Homo_sapiens.ensembl.100.tx2gene.tsv'],
    ['txome_k31_no_sa',
     's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-100/salmon_index/txome_k31',
     's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-100/annotation/Homo_sapiens.ensembl.100.tx2gene.tsv'],
    ['txome_k23_no_sa',
     's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-100/salmon_index/txome_k23',
     's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-100/annotation/Homo_sapiens.ensembl.100.tx2gene.tsv'],
    ['cdna_k31_full_sa',
     's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-100/salmon_index/cdna_k31_full_sa',
     's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-100/annotation/Homo_sapiens.ensembl.100.tx2gene.tsv'],
    ['txome_k31_full_sa',
     's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-100/salmon_index/txome_k31_full_sa',
     's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-100/annotation/Homo_sapiens.ensembl.100.tx2gene.tsv'],
    ['cdna_k31_partial_sa',
     's3://nextflow-ccdl-data/reference/homo_sapiens/refgenomes-hg38/salmon_partial_sa_index',
     's3://nextflow-ccdl-data/reference/homo_sapiens/refgenomes-hg38/annotation/Homo_sapiens.ensembl.97.tx2gene.tsv'],
  ])
  ch_testset = ch_reads.combine(ch_indexes)

  // run Alevin
  alevin(ch_testset)
}
