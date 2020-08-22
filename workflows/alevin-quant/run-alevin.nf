#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// run parameters
params.index_dir = 's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-100/salmon_index'
params.index_name = 'cdna_k31'
params.t2g = 's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-100/annotation/Homo_sapiens.ensembl.100.tx2gene.txt'
params.mitolist = 's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-100/annotation/Homo_sapiens.ensembl.100.mitogenes.txt'
params.sample_dir = 's3://ccdl-scpca-data/raw/green_adam'
params.sample_id = ['834', '905_3']



// process alevin{
//   input:
//     infiles
//   output:

// }

workflow{
  ch_files = Channel.fromList(params.sample_id)
    .map{ id -> tuple("$id",
                      file("${params.sample_dir}/${id}/*_R1_*.fastq.gz"),
                      file("${params.sample_dir}/${id}/*_R2_*.fastq.gz"),
                      )}
    .view()
}