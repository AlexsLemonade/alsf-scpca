#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// run parameters
params.ref_dir = 's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-100'
params.index_dir = 'kallisto_index'
params.index_name = 'cdna_k31'
params.annotation_dir = 'annotation'
params.t2g = 'Homo_sapiens.ensembl.100.tx2gene.tsv'
params.mitolist = 'Homo_sapiens.ensembl.100.mitogenes.txt'
params.barcodes = 's3://nextflow-ccdl-data/reference/10X/barcodes/3M-february-2018.txt' // 10X v3 barcodes

params.sample_dir = 's3://ccdl-scpca-data/raw/green_adam'
params.sample_ids = "834,905_3" //comma separated list to be parsed into a list
params.outdir = 's3://nextflow-ccdl-results/scpca/kallisto-quant'

// build full paths
params.index_path = "${params.ref_dir}/${params.index_dir}/${params.index_name}"
params.t2g_path = "${params.ref_dir}/${params.annotation_dir}/${params.t2g}"
params.mito_path = "${params.ref_dir}/${params.annotation_dir}/${params.mitolist}"


process kallisto_bus{
  container 'quay.io/biocontainers/kallisto:0.46.2--h4f7b962_1'
  label 'cpus_8'
  tag "${id}-${index}"
  publishDir "${params.outdir}"
  input:
    tuple val(id), path(read1), path(read2)
    path index
  output:
    path bus_dir
  script:
    // interleave read1 & read2 files
    bus_dir = "${id}-${index}"
    reads = [read1, read2].transpose().flatten().join(' ')
    """
    kallisto bus \
      -i ${index} \
      -o ${bus_dir} \
      -x 10xv3 \
      -t ${task.cpus} \
      ${reads}
    """
}

// process to create a barcode whitelist if not provided
// not currently used
process bustools_whitelist{
  container 'quay.io/biocontainers/bustools:0.40.0--h4f7b962_0'
  label 'cpus_8'
  publishDir "${params.outdir}"
  input:
    path bus_dir
  output:
    path whitelist
  script:
    whitelist = "${bus_dir}/whitelist.txt"
    """
    bustools sort \
      -o output.sorted.bus \
      -t ${task.cpus} \
      ${bus_dir}/output.bus
    bustools whitelist \
    -o ${whitelist} \
    output.sorted.bus
    """
}

process bustools_correct{
  container 'quay.io/biocontainers/bustools:0.40.0--h4f7b962_0'
  label 'cpus_8'
  input:
    path bus_dir
    path whitelist
  output:
    path bus_dir
  script:
    """
    bustools correct \
      -p \
      -w ${whitelist} \
      ${bus_dir}/output.bus \
    | bustools sort \
      -o ${bus_dir}/output.corrected.bus \
      -t ${task.cpus} \
      -
    """
}

process bustools_count{
  container 'quay.io/biocontainers/bustools:0.40.0--h4f7b962_0'
  label 'cpus_8'
  publishDir "${params.outdir}"
  input:
    path bus_dir
    path tx2gene
  output:
    path bus_dir
  script:
    """
    bustools count \
      -o ${bus_dir}/gene.count \
      -e ${bus_dir}/matrix.ec \
      -t ${bus_dir}/transcripts.txt \
      -g ${tx2gene} \
      --genecounts \
      ${bus_dir}/output.corrected.bus
    """
}


workflow{
  sample_ids = params.sample_ids?.tokenize(',') ?: []
  ch_reads = Channel.fromList(sample_ids)
    // create tuple of [sample_id, [Read1 files], [Read2 files]]
    .map{ id -> tuple("$id",
                      file("${params.sample_dir}/${id}/*_R1_*.fastq.gz"),
                      file("${params.sample_dir}/${id}/*_R2_*.fastq.gz"),
                      )}
  // run kallisto
  kallisto_bus(ch_reads, params.index_path)
  // correct busfiles
  bustools_correct(kallisto_bus.out, params.barcodes)
  // count genes
  bustools_count(bustools_correct.out, params.t2g_path)
}