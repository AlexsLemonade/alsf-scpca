#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// run parameters
params.resolution = 'cr-like' //default resolution is cr-like, can also use full, cr-like-em, parsimony, and trivial
params.barcode_dir = 's3://nextflow-ccdl-data/reference/10X/barcodes' 
params.run_metafile = 's3://ccdl-scpca-data/sample_info/scpca-library-metadata.tsv'
params.outdir = "s3://nextflow-ccdl-results/scpca/alevin-fry-features"

// 10X barcode files
barcodes = ['CITEseq_10Xv2': '737K-august-2016.txt',
            'CITEseq_10Xv3': '3M-february-2018.txt',
            'CITEseq_10Xv3.1': '3M-february-2018.txt']
// supported single cell technologies
tech_list = barcodes.keySet()

// run_ids are comma separated list to be parsed into a list of run ids,
// or "All" to process all samples in the metadata file
params.run_ids = "SCPCR000084,SCPCR000085,SCPCR000340"

// containers
SALMON_CONTAINER = 'quay.io/biocontainers/salmon:1.5.2--h84f40af_0'
ALEVINFRY_CONTAINER = 'quay.io/biocontainers/alevin-fry:0.4.1--h7d875b9_0'

process index_feature{
  container SALMON_CONTAINER
  
  input:
    tuple val(id), path(feature_file)
  output:
    tuple val(id), path("feature_index")
  script:
    """
    salmon index \
      -t ${feature_file} \
      -i feature_index \
      --features \
      -k 7 

    awk '{print \$1"\\t"\$1;}' ${feature_file} > feature_index/t2g.tsv
    """
}

// generates RAD file for alevin feature matrix using alevin
process alevin_feature{
  container SALMON_CONTAINER
  label 'cpus_8'
  tag "${id}-features"
  input:
    tuple val(id), val(tech), path(read1), path(read2), val(feature_geom), path(feature_index)
  output:
    tuple val(id), path(run_dir), path(feature_index)
  script:
    // label the run directory by id
    run_dir = "${id}-features"
    // Define umi geometries
    umi_geoms = ['CITEseq_10Xv2': '1[17-26]',
                 'CITEseq_10Xv3': '1[17-28]',
                 'CITEseq_10Xv3.1': '1[17-28]']
    """
    mkdir -p ${run_dir}
    salmon alevin \
      -l ISR \
      -1 ${read1} \
      -2 ${read2} \
      -i ${feature_index} \
      --read-geometry ${feature_geom} \
      --bc-geometry 1[1-16] \
      --umi-geometry ${umi_geoms[tech]} \
      --sketch \
      --rad \
      -o ${run_dir} \
      -p ${task.cpus} 
    """
}

// quantify from rad input
process permit_collate_quant_feature{
  container ALEVINFRY_CONTAINER
  label 'cpus_8'
  publishDir "${params.outdir}"

  input:
    tuple val(id), path(run_dir), path(feature_index)
    path barcode_file
  output:
    path run_dir
  
  script: 
    """
    alevin-fry generate-permit-list \
      -i ${run_dir} \
      --expected-ori fw \
      -o ${run_dir} \
      -u ${barcode_file}

    alevin-fry collate \
      --input-dir ${run_dir} \
      --rad-dir ${run_dir} \
      -t ${task.cpus}
    
    alevin-fry quant \
      --input-dir ${run_dir} \
      --tg-map ${feature_index}/t2g.tsv \
      -r ${params.resolution} \
      -o ${run_dir} \
      --use-mtx \
      -t ${task.cpus} \
      
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

  //get and map the feature barcode files
  features_ch = samples_ch
    .map{row -> tuple(row.feature_barcode_file,
                      file("s3://${row.feature_barcode_file}"))}
    .unique()
  index_feature(features_ch)

  // create tuple of [sample_id, technology, [Read1 files], [Read2 files], feature_geometry, feature_index]
  reads_ch = samples_ch
    .map{row -> tuple(row.feature_barcode_file,
                      row.scpca_run_id,
                      row.technology,
                      file("s3://${row.s3_prefix}/*_R1_*.fastq.gz"),
                      file("s3://${row.s3_prefix}/*_R2_*.fastq.gz"),
                      row.feature_barcode_geom
                      )}
    .combine(index_feature.out, by: 0) // combine by the feature_barcode_file
    .map{ it.subList(1, it.size())} // remove the first element

  // run Alevin
  alevin_feature(reads_ch)
  
  barcodes_ch = samples_ch
    .map{row -> file("${params.barcode_dir}/${barcodes[row.technology]}")}
  // generate permit list from alignment 
  permit_collate_quant_feature(alevin_feature.out, barcodes_ch)
}
