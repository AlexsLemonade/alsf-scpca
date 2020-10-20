#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// basic parameters
params.ref_dir = 's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-100'
params.index_dir = 'cellranger_index'
params.index_name = 'cdna'

params.run_metafile = 's3://ccdl-scpca-data/sample_info/scpca-library-metadata.tsv'
params.run_ids = "SCPCR000001" //comma separated list to be parsed into a list
params.outdir = 's3://nextflow-ccdl-results/scpca/cellranger-quant'

// build full paths
params.index_path = "${params.ref_dir}/${params.index_dir}/${params.index_name}"

process cellranger{
  container '589864003899.dkr.ecr.us-east-1.amazonaws.com/scpca-cellranger:4.0.0'
  publishDir "${params.outdir}", mode: 'copy'
  label 'cpus_8'
  label 'bigdisk'
  input:
    tuple val(id), val(samples), path(fastq_dir)
    path index
  output:
    path output_id
  script:
    output_id = "${id}-${index}"
    """
    cellranger count \
      --id=${output_id} \
      --transcriptome=${index} \
      --fastqs=${fastq_dir} \
      --sample=${samples} \
      --localcores=${task.cpus} \
      --localmem=${task.memory.toGiga()}
    """
}

def getCRsamples(filelist){
  // takes a string with semicolon separated file names
  // returns just the 'sample info' portion of the file names,
  // as cellranger would interpret them, comma separated
  fastq_files = filelist.tokenize(';')
  samples = []
  fastq_files.each{
    // append sample names to list, using regex to extract element before S001, etc.
    // [0] for the first match set, [1] for the first extracted element
    samples << (it =~ /^(.+)_S.+_L.+_R.+.fastq.gz$/)[0][1]
  }
  // set to remove dup
  return samples.toSet().join(',')
}


workflow{
  run_ids = params.run_ids?.tokenize(',') ?: []
  ch_reads = Channel.fromPath(params.run_metafile)
    .splitCsv(header: true, sep: '\t')
    .filter{it.scpca_run_id in run_ids} // use only the rows in the sample list
    // create tuple of [sample_id, fastq dir]
    .map{row -> tuple(row.scpca_run_id,
                      getCRsamples(row.files),
                      file("s3://${row.s3_prefix}")
                      )}
  // run cellranger
  cellranger(ch_reads, params.index_path)
}
