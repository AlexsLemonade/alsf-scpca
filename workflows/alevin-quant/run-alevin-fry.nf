#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// run parameters
params.ref_dir = 's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-103'
params.index_dir = 'salmon_index'
params.annotation_dir = 'annotation'
params.index_type = 'cdna' // default index type is cdna, can also use splici
params.sketch = false // use sketch mode for mapping with flag `--sketch`
params.resolution = 'full' //default resolution is full, can also use cr-like, cr-like-em, parsimony, and trivial
params.t2g_3col = 'Homo_sapiens.GRCh38.103.spliced_intron.tx2gene_3col.tsv'
params.barcode_dir = 's3://nextflow-ccdl-data/reference/10X/barcodes' 
params.filter = 'unfiltered' // default filtering strategy is to use unfiltered, other options include knee (--knee-distance)
params.run_metafile = 's3://ccdl-scpca-data/sample_info/scpca-library-metadata.tsv'
params.outdir = "s3://nextflow-ccdl-results/scpca/alevin-fry-${params.filter}-quant"

// run_ids are comma separated list to be parsed into a list of run ids,
// or "All" to process all samples in the metadata file
params.run_ids = "SCPCR000001,SCPCR000002"

// index files 
index_names_map = ['cdna': 'spliced_txome_k31',
                   'splici': 'spliced_intron_txome_k31']

// tx2gene files
t2g_map = ['cdna': 'Homo_sapiens.GRCh38.103.spliced.tx2gene.tsv',
           'splici': 'Homo_sapiens.GRCh38.103.spliced_intron.tx2gene.tsv']

// 10X barcode files
barcodes = ['10Xv2': '737K-august-2016.txt',
            '10Xv3': '3M-february-2018.txt',
            '10Xv3.1': '3M-february-2018.txt']

// supported single cell technologies
tech_list = barcodes.keySet()

// file paths
index_path = "${params.ref_dir}/${params.index_dir}/${index_names_map[params.index_type]}"
index_prefix = "${params.ref_dir}/${params.annotation_dir}"
t2g_path = "${index_prefix}/${t2g_map[params.index_type]}"
  
// if using splici and cr-like use the 3 column t2g file for alevin-fry quant
if(params.resolution == 'cr-like' && params.index_type == 'splici'){
    t2g_quant_path = "${index_prefix}/${params.t2g_3col}"
    use_mtx = true
} else{
    t2g_quant_path = t2g_path
    use_mtx = false
}

// generates RAD file using alevin
process alevin{
  container 'quay.io/biocontainers/salmon:1.4.0--hf69c8f4_0'
  label 'cpus_8'
  tag "${id}-${index}"
  publishDir "${params.outdir}"
  input:
    tuple val(id), val(tech), path(read1), path(read2)
    path index
    path tx2gene
  output:
    path run_dir
  script:
    // label the run directory by id, index, and mapping mode
    run_dir = "${id}-${index}-${params.sketch ? 'sketch' : 'salign'}-${params.resolution}"
        if( params.filter != 'unfiltered'){
          run_dir += "-${params.filter}"
        }
    // choose flag by technology
    tech_flag = ['10Xv2': '--chromium',
                 '10Xv3': '--chromiumV3',
                 '10Xv3.1': '--chromiumV3']
    // run alevin like normal with the --justAlign flag 
    // creates output directory with RAD file needed for alevin-fry
    // uses sketch mode if --sketch was included at invocation
    """
    mkdir -p ${run_dir}
    salmon alevin \
      -l ISR \
      ${tech_flag[tech]} \
      -1 ${read1} \
      -2 ${read2} \
      -i ${index} \
      --tgMap ${tx2gene} \
      -o ${run_dir} \
      -p ${task.cpus} \
      --dumpFeatures \
      --justAlign \
      ${params.sketch ? '--sketch' : ''}
    """
}

//generate permit list from RAD input 
process generate_permit{
  container 'ghcr.io/alexslemonade/scpca-alevin-fry:latest'
  publishDir "${params.outdir}"
  input:
    path run_dir
    path barcode_file
  output:
    path run_dir
  script: 
    filter_map = ['unfiltered' : "-u ${barcode_file}",
                  'knee' : '--knee-distance']
    """
    alevin-fry generate-permit-list \
      -i ${run_dir} \
      --expected-ori fw \
      -o ${run_dir} \
      ${filter_map[params.filter]}
    """
}

// given permit list and barcode mapping, collate RAD file 
process collate_fry{
  container 'ghcr.io/alexslemonade/scpca-alevin-fry:latest'
  label 'cpus_8'
  publishDir "${params.outdir}"
  input: 
    path run_dir
  output: 
    path run_dir
  script:
    """
    alevin-fry collate \
      --input-dir ${run_dir} \
      --rad-dir ${run_dir} \
      -t ${task.cpus}
    """
}

// then quantify collated RAD file
process quant_fry{
  container 'ghcr.io/alexslemonade/scpca-alevin-fry:latest'
  label 'cpus_8'
  publishDir "${params.outdir}"
  input: 
    path run_dir
    path tx2gene
    val use_mtx
  output: 
    path run_dir
  script:
    """
    alevin-fry quant \
     --input-dir ${run_dir} \
     --tg-map ${tx2gene} \
     --output-dir ${run_dir} \
     -r ${params.resolution} \
     -t ${task.cpus} \
     ${use_mtx ? '--use-mtx' : ''}
    """
}
// run quant command with default full resolution strategy 

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

  barcodes_ch = samples_ch
    .map{row -> file("${params.barcode_dir}/${barcodes[row.technology]}")}

  // run Alevin
  alevin(reads_ch, index_path, t2g_path)
  // generate permit list from alignment 
  generate_permit(alevin.out, barcodes_ch)
  // collate RAD files 
  collate_fry(generate_permit.out)
  // create gene x cell matrix
  quant_fry(collate_fry.out, t2g_quant_path, use_mtx)
}
