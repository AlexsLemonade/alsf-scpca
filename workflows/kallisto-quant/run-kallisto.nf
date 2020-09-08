#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// run parameters
params.ref_dir = 's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-100'
params.index_dir = 'kallisto_index'
params.index_name = 'cdna_k31'
params.annotation_dir = 'annotation'
params.t2g = 'Homo_sapiens.ensembl.100.tx2gene.tsv'
params.mitolist = 'Homo_sapiens.ensembl.100.mitogenes.txt'

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
    path bus
  script:
    // interleave read1 & read2 files
    bus = "${id}-${index}_bus"
    reads = [read1, read2].transpose().flatten().join(' ')
    """
    kallisto bus \
      -i ${index} \
      -o ${bus} \
      -x 10xv3 \
      -t ${task.cpus} \
      ${reads}
    """
}

process bustools_sort{
  container 'quay.io/biocontainers/bustools:0.40.0--h4f7b962_0'
  label 'cpus_8'
  input:
    path bus
  output:
    path outfile
  script:
    outfile = "${bus}/output.sorted.bus"
    """
    bustools sort \
      -o  ${outfile} \
      -t ${task.cpus} \
      -m ${task.memory.toGiga()}G \
      ${bus}/output.bus
    """
}

process bustools_whitelist{
  container 'quay.io/biocontainers/bustools:0.40.0--h4f7b962_0'
  input:
    path busfile_sorted
  output:
    path whitelist_file
  script:
    whitelist_file = "${busfile_sorted.simpleName}_whitelist.txt"
    """
    bustools whitelist \
      -o ${whitelist_file} \
      ${busfile_sorted}
    """
}

process bustools_correct{
  container 'quay.io/biocontainers/bustools:0.40.0--h4f7b962_0'
  label 'cpus_8'
  input:
    path bus
    path whitelist
  output:
    path outfile
  script:
    outfile = "${bus}/output.corrected.bus"
    """
    bustools correct \
      -p \
      -w ${whitelist} \
      ${bus}/output.bus \
    | bustools sort \
      -o ${outfile} \
      -t ${task.cpus} \
      -m ${task.memory.toGiga()}G \
      -
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
  // run Kallisto
  kallisto_bus(ch_reads, params.index_path)
  // get busfiles
  // generate whitelist
  bustools_sort(kallisto_bus.out) | bustools_whitelist
  // correct busfiles
  bustools_correct(kallisto_bus.out, bustools_whitelist.out)
}