#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// basic parameters
params.ref_dir = 's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-100'
params.index_dir = 'cellranger_index'
params.index_name = 'cdna'

params.run_metafile = 's3://ccdl-scpca-data/sample_info/scpca-library-metadata.tsv'
// run_ids are comma separated list to be parsed into a list of run ids,
// or "All" to process all samples in the metadata file
params.run_ids = "SCPCR000001,SCPCR000002"
params.outdir = 's3://nextflow-ccdl-results/scpca/cellranger-quant'

// technology options
cellranger_tech_list = ["10Xv2", "10Xv3", "10Xv3.1", "10Xv2_5prime"]
all_tech_list = cellranger_tech_list + "spatial"

// build full paths
params.index_path = "${params.ref_dir}/${params.index_dir}/${params.index_name}"

process cellranger{
  container '589864003899.dkr.ecr.us-east-1.amazonaws.com/scpca-cellranger:6.1.2'
  publishDir "${params.outdir}", mode: 'copy'
  tag "${id}-${index}" // add tag for tracking sample names in trace file  
  label 'cpus_8'
  label 'bigdisk'
  input:
    tuple val(meta), path(fastq_dir)
    path index
  output:
    path output_id
  script:
    output_id = "${meta.library_id}-${index}-${meta.include_introns ? 'pre_mRNA' : 'mRNA'}"
    """
    cellranger count \
      --id=${output_id} \
      --transcriptome=${index} \
      --fastqs=${fastq_dir} \
      --sample=${meta.cr_samples} \
      --localcores=${task.cpus} \
      --localmem=${task.memory.toGiga()} \
      ${meta.include_introns ? '--include-introns' : ''}

    """
}

process spaceranger{
  container '589864003899.dkr.ecr.us-east-1.amazonaws.com/scpca-cellranger:6.0.1'
  publishDir "${params.outdir}", mode: 'copy'
  tag "${id}-${index}-spatial" 
  label 'cpus_8'
  label 'bigdisk'
  input:
    tuple val(meta), path(fastq_dir), file(image_file)
    path index
  output:
    path output_id
  script:
    output_id = "${meta.library_id}-${index}-spatial"
    """
    spaceranger count \
      --id=${output_id} \
      --transcriptome=${index} \
      --fastqs=${fastq_dir} \
      --sample=${meta.samples} \
      --localcores=${task.cpus} \
      --localmem=${task.memory.toGiga()} \
      --image=${image_file} \
      --slide=${meta.slide_serial_number} \
      --area=${meta.slide_section}

    """
}

def getCRsamples(filelist){
  // takes a string with semicolon separated file names
  // returns just the 'sample info' portion of the file names,
  // as cellranger would interpret them, comma separated
  fastq_files = filelist.tokenize(';').findAll{it.contains '.fastq.gz'}
  samples = []
  fastq_files.each{
    // append sample names to list, using regex to extract element before S001, etc.
    // [0] for the first match set, [1] for the first extracted element
    samples << (it =~ /^(.+)_S.+_L.+_[R|I].+.fastq.gz$/)[0][1]
  }
  // convert samples list to a `set` to remove duplicate entries,
  // then join to a comma separated string.
  return samples.toSet().join(',')
}


workflow{
  run_ids = params.run_ids?.tokenize(',') ?: []
  run_all = run_ids[0] == "All"
  ch_reads = Channel.fromPath(params.run_metafile)
    .splitCsv(header: true, sep: '\t')
    .filter{it.technology in all_tech_list}
    // use only the rows in the sample list
    .filter{run_all || (it.scpca_run_id in run_ids)}
    .map{it.cr_samples =  getCRsamples(it.files); it}
  
  cellranger_reads = ch_reads
    .filter{it.technology in cellranger_tech_list} // only cellranger 10X data
    // create tuple of [metadata, fastq dir]
    //.map{it.cr_samples =  getCRsamples(it.files); it}
    .map{meta -> tuple(meta,
                       file("s3://${meta.s3_prefix}")
                       )}

  spaceranger_reads = ch_reads
    .filter{it.technology == "spatial"}
    // create tuple of [metadata, fastq dir, and image filename]
    .map{meta -> tuple(meta,
                       file("s3://${meta.s3_prefix}"),
                       //getCRsamples(${meta.files}.findAll{it.contains '.fastq.gz'}),
                       file("s3://${meta.s3_prefix}/*.jpg")
                       )}

  // run cellranger
  cellranger(cellranger_reads, params.index_path)

  // run spaceranger 
  spaceranger(spaceranger_reads, params.index_path)
  
}
