#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// basic parameters
params.ref_dir = 's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-103'
params.spliced_txome = 'fasta/Homo_sapiens.GRCh38.103.spliced.txome.fa.gz'
params.spliced_intron_txome = 'fasta/Homo_sapiens.GRCh38.103.spliced_intron.txome.fa.gz'
params.ensembl_txome = 'fasta/Homo_sapiens.GRCh38.txome.fa.gz'
params.genome = 'fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz'
params.kmer = '31'


// remove .gz if it exists
def get_base(file){
  if (file.extension == "gz"){
    return file.baseName
  }else{
    return file.fileName
  }
}


process salmon_index_no_sa{
  container 'quay.io/biocontainers/salmon:1.4.0--hf69c8f4_0'
  publishDir "${params.ref_dir}/salmon_index", mode: 'copy'
  memory { 28.GB * task.attempt}
  cpus 8
  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
  maxRetries 1
  input:
    tuple val(index_base), path(reference), val(kmer)
  output:
    path "${index_base}_k${kmer}"
  script:
    """
    salmon index \
      -t ${reference} \
      -i ${index_base}_k${kmer} \
      -k ${kmer} \
      -p ${task.cpus} \
    """
}

process salmon_index_full_sa{
  container 'quay.io/biocontainers/salmon:1.4.0--hf69c8f4_0'
  publishDir "${params.ref_dir}/salmon_index", mode: 'copy'
  // try dynamic memory (28.GB so 2x will fit in r4.2xlarge)
  memory { 28.GB * task.attempt}
  cpus 8
  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
  maxRetries 1
  input:
    tuple val(index_base), path(reference), val(kmer)
    path genome
  output:
    path "${index_base}_k${kmer}_full_sa"
  script:
    """
    gunzip -c ${genome} \
      |grep "^>" | cut -d " " -f 1 \
      |sed -e 's/>//g' > decoys.txt
    cat ${reference} ${genome} > gentrome.fa.gz
    salmon index \
      -t gentrome.fa.gz \
      -d decoys.txt \
      -i ${index_base}_k${kmer}_full_sa \
      -k ${kmer} \
      -p ${task.cpus} \
    """
}

process kallisto_index{
  container 'quay.io/biocontainers/kallisto:0.46.2--h4f7b962_1'
  publishDir "${params.ref_dir}/kallisto_index", mode: 'copy'
  memory { 120.GB * task.attempt}
  cpus 8
  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
  maxRetries 1
  input:
    tuple val(index_base), path(reference), val(kmer)
  output:
    path "${index_base}_k${kmer}"
  script:
    """
    kallisto index \
      -i ${index_base}_k${kmer} \
      -k ${kmer} \
      ${reference}
    """
}

workflow {
  // channel of the reference files and labels
  ch_ref = Channel
    .fromList([["ensembl_txome", params.ref_dir + "/" + params.ensembl_txome, params.kmer],
               ["spliced_txome", params.ref_dir + "/" + params.spliced_txome, params.kmer],
               ["spliced_intron_txome", params.ref_dir + "/" + params.spliced_intron_txome, params.kmer]])

  salmon_index_no_sa(ch_ref)
  salmon_index_full_sa(ch_ref, params.ref_dir + "/" + params.genome)
  kallisto_index(ch_ref)
}
