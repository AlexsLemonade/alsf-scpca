#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// run parameters
params.ref_dir = 's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-100'
params.index_dir = 'salmon_index'
params.index_name = 'txome_k31_full_sa'
params.annotation_dir = 'annotation'
params.t2g = 'Homo_sapiens.ensembl.100.tx2gene.tsv'
params.mitolist = 'Homo_sapiens.ensembl.100.mitogenes.txt'

params.run_metafile = 's3://ccdl-scpca-data/sample_info/scpca-library-metadata.tsv'
// run_ids are comma separated list to be parsed into a list of run ids,
// or "All" to process all samples in the metadata file
params.run_ids = "SCPCR000001,SCPCR000002"

params.outdir = 's3://nextflow-ccdl-results/scpca/alevin-fry-quant'

// build full paths
params.index_path = "${params.ref_dir}/${params.index_dir}/${params.index_name}"
params.t2g_path = "${params.ref_dir}/${params.annotation_dir}/${params.t2g}"
params.mito_path = "${params.ref_dir}/${params.annotation_dir}/${params.mitolist}"

// generates RAD file using alevin
process alevin{
  container 'quay.io/biocontainers/salmon:1.4.0--hf69c8f4_0'
  label 'cpus_8'
  tag "${id}-${index}"
  publishDir "${params.outdir}"
  input:
    tuple val(id), path(read1), path(read2)
    path index
    path tx2gene
  output:
    path run_dir
  script:
    run_dir = "${id}-${index}"
  // run alevin like normal with the --justAlign flag 
  // creates output directory with RAD file needed for alevin-fry
    """
    mkdir -p ${run_dir}/alevin
    salmon alevin \
      -l ISR \
      --chromium \
      -1 ${read1} \
      -2 ${read2} \
      -i ${index} \
      --tgMap ${tx2gene} \
      -o ${run_dir}/alevin \
      -p ${task.cpus} \
      --dumpFeatures \
      --justAlign
    """
}

//generate permit list from RAD input 
process generate_permit{
  container 'ghcr.io/alexslemonade/scpca-alevin-fry:latest'
  label 'cpus_8'
  publishDir "${params.outdir}"
  input:
    path run_dir
  output:
    path permit_list
  script: 
    permit_list = "${run_dir}/permit_list"
  // expected-ori either signifies no filtering of alignments 
  // based on orientation, not sure if this is what we want? 
    """
    alevin-fry generate-permit-list \
      -i ${run_dir}/alevin \
      --expected-ori either \
      -o ${permit_list} \
      -k
    """
}

// given permit list and barcode mapping, collate RAD file 
process collate{
  container 'ghcr.io/alexslemonade/scpca-alevin-fry:latest'
  label 'cpus_8'
  publishDir "${params.outdir}"
  input: 
    path permit_list
    path run_dir
  output: 
    file collate_file
  script:
    collate_file = "${run_dir}/permit_list/map.collated.rad"
    """
    alevin-fry collate \
      --input-dir ${permit_list} \
      --rad-dir ${run_dir}/alevin \
      -t ${task.cpus}
    """
}

// quantify collated RAD file
process quant{
  container 'ghcr.io/alexslemonade/scpca-alevin-fry:latest'
  label 'cpus_8'
  publishDir "${params.outdir}"
  input: 
    file collate_file
    path run_dir
    path tx2gene
  output:
    path counts
  script:
    collate_file = "${run_dir}/permit_list/map.collated.rad"
    counts = "${run_dir}/counts"
    // run quant command with default full resolution strategy 
    // gene-level, UMI-deduplicated, equivalence class counts file + gene counts matrix output
    """
      mkdir -p ${run_dir}/counts

      alevin-fry quant \
        -i ${collate_file} \
        --tg-map ${tx2gene} \
        -o ${counts} \
        -t ${task.cpus} \
        --dump-eqclasses
    """
}


workflow{
  run_ids = params.run_ids?.tokenize(',') ?: []
  run_all = run_ids[0] == "All"
  ch_reads = Channel.fromPath(params.run_metafile)
    .splitCsv(header: true, sep: '\t')
    .filter{it.technology == "10Xv3"} // only 10X data
    // use only the rows in the sample list
    .filter{run_all || (it.scpca_run_id in run_ids)}
    // create tuple of [sample_id, [Read1 files], [Read2 files]]
    .map{row -> tuple(row.scpca_run_id,
                      file("s3://${row.s3_prefix}/*_R1_*.fastq.gz"),
                      file("s3://${row.s3_prefix}/*_R2_*.fastq.gz"),
                      )}
  // run Alevin
  alevin(ch_reads, params.index_path, params.t2g_path)
  // generate permit list from alignment 
  generate_permit(alevin.out)
  // collate RAD files
  collate(generate_permit.out, alevin.out)
  // create gene x cell matrix
  quant(collate.out, alevin.out, params.t2g_path)
}