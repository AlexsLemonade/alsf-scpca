#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// basic parameters
params.ref_dir = 's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-100'
params.cdna = 'fasta/Homo_sapiens.GRCh38.cdna.all.fa.gz'
params.ncrna = 'fasta/Homo_sapiens.GRCh38.ncrna.fa.gz'
params.k = [23, 31]


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
    tuple path(reference), val(index_base), val(k)
  output:
    path "${index_base}_k${k}"
  script:
    """
    salmon index \
      -t ${reference} \
      -i ${index_base}_k${k} \
      -k ${k}
    """
}

workflow {
  // channel of the reference file(s) and a label
  ch_ref = Channel
    .fromList([[params.ref_dir + "/" + params.cdna, "cdna"]])
  // possible k values
  ch_k = Channel.fromList(params.k)
  // create a channel with all k values for each ref
  ch_salmon = ch_ref.combine(ch_k)

  salmon_index(ch_salmon)
}