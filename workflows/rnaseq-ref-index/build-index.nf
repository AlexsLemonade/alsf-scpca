#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// basic parameters
params.ref_dir = 's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-100'
params.cdna = 'fasta/Homo_sapiens.GRCh38.cdna.all.fa.gz'
params.txome = 'fasta/Homo_sapiens.GRCh38.txome.fa.gz'
params.kmer = [23, 31]


// remove .gz if it exists
def get_base(file){
  if (file.extension == "gz"){
    return file.baseName
  }else{
    return file.fileName
  }
}


process salmon_index{
  container 'quay.io/biocontainers/salmon:1.3.0--hf69c8f4_0'
  publishDir "${params.ref_dir}/salmon_index" , mode: 'copy'
  input:
    tuple path(reference), val(index_base), val(kmer)
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

process kallisto_index{
  container 'quay.io/biocontainers/kallisto:0.46.2--h4f7b962_1'
  publishDir "${params.ref_dir}/kallisto_index" , mode: 'copy'
  memory 32.GB
  input:
    tuple path(reference), val(index_base), val(kmer)
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
    .fromList([[params.ref_dir + "/" + params.cdna, "cdna"],
               [params.ref_dir + "/" + params.txome, "txome"]])
  // possible kmer values
  ch_kmer = Channel.fromList(params.kmer)
  // create a channel with all k values for each ref
  ch_index = ch_ref.combine(ch_kmer)

  salmon_index(ch_index)
  kallisto_index(ch_index)
}
