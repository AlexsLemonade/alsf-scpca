#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// run parameters
params.ref_dir = 's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-103'
params.index_dir = 'salmon_index'
params.index_name = 'spliced_txome_k31_full_sa'
params.annotation_dir = 'annotation'
params.t2g = 'Homo_sapiens.GRCh38.103.spliced.tx2gene.tsv'
params.mitolist = 'Homo_sapiens.GRCh38.103.mitogenes.txt'

params.run_metafile = 's3://ccdl-scpca-data/sample_info/scpca-library-metadata.tsv'
// run_ids are comma separated list to be parsed into a list of run ids,
// or "All" to process all samples in the metadata file
params.run_ids = "SCPCR000001,SCPCR000002"

params.outdir = 's3://nextflow-ccdl-results/scpca/alevin-quant'

// build full paths
params.index_path = "${params.ref_dir}/${params.index_dir}/${params.index_name}"
params.t2g_path = "${params.ref_dir}/${params.annotation_dir}/${params.t2g}"
params.mito_path = "${params.ref_dir}/${params.annotation_dir}/${params.mitolist}"

// supported single cell technologies
tech_list = ['10Xv2', '10Xv3', '10Xv3.1'] 

process alevin{
  container 'quay.io/biocontainers/salmon:1.5.2--h84f40af_0'
  label 'cpus_8'
  tag "${id}-${index}"
  publishDir "${params.outdir}"
  input:
    tuple val(id), val(tech), path(read1), path(read2)
    path index
    path tx2gene
  output:
    path "${id}-${index}"
  script:
    // choose flag by technology
    tech_flag = ['10Xv2': '--chromium',
                 '10Xv3': '--chromiumV3',
                 '10Xv3.1': '--chromiumV3']
    """
    salmon alevin \
      -l ISR \
      ${tech_flag[tech]} \
      -1 ${read1} \
      -2 ${read2} \
      -i ${index} \
      --tgMap ${tx2gene} \
      -o ${id}-${index} \
      -p ${task.cpus} \
      --dumpFeatures
    """
}

workflow{
  run_ids = params.run_ids?.tokenize(',') ?: []
  run_all = run_ids[0] == "All"
  samples_ch = Channel.fromPath(params.run_metafile)
    .splitCsv(header: true, sep: '\t')
    .filter{it.technology in tech_list} 
    // use only the rows in the sample list
    .filter{run_all || (it.scpca_run_id in run_ids)}
  // create tuple of [sample_id, technology, [Read1 files], [Read2 files]]
  reads_ch = samples_ch
    .map{row -> tuple(row.scpca_run_id,
                      row.technology,
                      file("s3://${row.s3_prefix}/*_R1_*.fastq.gz"),
                      file("s3://${row.s3_prefix}/*_R2_*.fastq.gz"),
                      )}
  // run Alevin
  alevin(reads_ch, params.index_path, params.t2g_path)
}
