#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// basic parameters
params.ref_dir = 's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-100'
params.cdna = 'fasta/Homo_sapiens.GRCh38.cdna.all.fa.gz'
params.txome = 'fasta/Homo_sapiens.GRCh38.txome.fa.gz'
params.genome = 'Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz'
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
  container 'quay.io/biocontainers/salmon:1.3.0--hf69c8f4_0'
  publishDir "${params.ref_dir}/salmon_index", mode: 'copy'
  input:
    tuple val(index_base), path(reference), val(kmer)
  output:
    path "${index_base}_k${kmer}"
  script:
    """
    salmon index \
      -t ${reference} \
      -i ${index_base}_k${kmer} \
      -k ${kmer}
    """
}

process salmon_index_full_sa{
  container 'quay.io/biocontainers/salmon:1.3.0--hf69c8f4_0'
  publishDir "${params.ref_dir}/salmon_index", mode: 'copy'
  input:
    tuple val(index_base), path(reference), val(kmer)
  output:
    path "${index_base}_k${kmer}_full_sa"
  script:
    """
    gunzip -c ${params.genome} \
      |grep "^>" | cut -d " " -f 1 \
      |sed -e 's/>//g' > decoys.txt
    cat ${reference} ${params.genome} > gentrome.fa.gz
    salmon index \
      -t gentrome.fa.gz \
      -d decoys.txt \
      -i ${index_base}_k${kmer} \
      -k ${kmer}
    """
}

process kallisto_index{
  container 'quay.io/biocontainers/kallisto:0.46.2--h4f7b962_1'
  publishDir "${params.ref_dir}/kallisto_index", mode: 'copy'
  memory 32.GB
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
    .fromList([["cdna", params.ref_dir + "/" + params.cdna, params.kmer],
               ["txome", params.ref_dir + "/" + params.txome, params.kmer]])

  salmon_index_no_sa(ch_ref)
  salmon_index_full_sa(ch_ref)
  kallisto_index(ch_ref)
}