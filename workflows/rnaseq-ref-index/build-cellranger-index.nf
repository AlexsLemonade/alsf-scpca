#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// basic parameters
params.ref_dir = 's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-100'
params.genome = 'fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz'
params.gtf = 'annotation/Homo_sapiens.GRCh38.100.gtf.gz'


process cellranger_index{
  container '589864003899.dkr.ecr.us-east-1.amazonaws.com/scpca-cellranger:4.0.0'
  publishDir "${params.ref_dir}/cellranger_index", mode: 'copy'
  label 'cpus_8'
  input:
    tuple val(index_base), path(fasta), path(gtf)
  output:
    path "${index_base}"
  script:
    """
    gunzip -c ${fasta} > genome.fasta
    gunzip -c ${gtf} > genome.gtf
    cellranger mkgtf \
      genome.gtf \
      filtered.gtf \
      --attribute=gene_biotype:protein_coding

    cellranger mkref \
      --genome=${index_base} \
      --fasta=genome.fasta \
      --genes=filtered.gtf \
      --nthreads=${task.cpus}
    """
}

workflow {
  // channel of the reference files and labels
  ch_ref = Channel
    .fromList([ //currently only one reference, for testing
      ["cdna", params.ref_dir + "/" + params.genome, params.ref_dir + "/" + params.gtf ],
    ])
  cellranger_index(ch_ref)
}
