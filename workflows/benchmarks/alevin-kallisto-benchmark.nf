#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// run parameters
params.sample_dir = 's3://ccdl-scpca-data/raw/green_adam'
params.sample_ids = "834,905_3" //comma separated values to be parsed into a list
params.outdir = 's3://nextflow-ccdl-results/scpca-benchmark/quants'
params.t2g = 's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-100/annotation/Homo_sapiens.ensembl.100.tx2gene.tsv'
params.barcodes = 's3://nextflow-ccdl-data/reference/10X/barcodes/3M-february-2018.txt'

process alevin{
  container 'quay.io/biocontainers/salmon:1.3.0--hf69c8f4_0'
  label 'cpus_8'
  tag "${sample_id}-${index_id}"
  publishDir "${params.outdir}/alevin"
  input:
    tuple val(sample_id), path(read1), path(read2), val(index_id), path(index)
    path(tx2gene)
  output:
    path "${sample_id}-${index_id}"
  script:
    """
    salmon alevin \
      -l ISR \
      --chromium \
      -1 ${read1} \
      -2 ${read2} \
      -i ${index} \
      --tgMap ${tx2gene} \
      -o ${sample_id}-${index_id} \
      -p ${task.cpus} \
      --dumpFeatures \
    """
}

process kallisto_bus{
  container 'quay.io/biocontainers/kallisto:0.46.2--h4f7b962_1'
  label 'cpus_8'
  tag "${sample_id}-${index_id}"
  publishDir "${params.outdir}/kallisto"
  input:
    tuple val(sample_id), path(read1), path(read2), val(index_id), path(index)
  output:
    path run_dir
  script:
    run_dir = "${sample_id}-${index_id}"
    // interleave read1 & read2 file paths
    reads = [read1, read2].transpose().flatten().join(' ')
    """
    mkdir -p ${run_dir}/bus
    kallisto bus \
      -i ${index} \
      -o ${run_dir}/bus \
      -x 10xv3 \
      -t ${task.cpus} \
      ${reads}
    """
}

process bustools_correct{
  container 'quay.io/biocontainers/bustools:0.40.0--h4f7b962_0'
  label 'cpus_8'
  tag "${run_dir}-buscorrect"
  publishDir "${params.outdir}/kallisto"
  input:
    path run_dir
    path whitelist
  output:
    path run_dir
  script:
    """
    bustools correct \
      -p \
      -w ${whitelist} \
      ${run_dir}/bus/output.bus \
    | bustools sort \
      -o ${run_dir}/bus/output.corrected.bus \
      -t ${task.cpus} \
      -
    """
}

process bustools_count{
  container 'quay.io/biocontainers/bustools:0.40.0--h4f7b962_0'
  label 'cpus_8'
  tag "${run_dir}-buscount"
  publishDir "${params.outdir}/kallisto"
  input:
    path run_dir
    path tx2gene
  output:
    path count_dir
  script:
    count_dir ="${run_dir}/counts"
    """
    mkdir -p ${count_dir}
    bustools count \
      -o ${count_dir}/gene_count \
      -e ${run_dir}/bus/matrix.ec \
      -t ${run_dir}/bus/transcripts.txt \
      -g ${tx2gene} \
      --genecounts \
      ${run_dir}/bus/output.corrected.bus
    """
}

workflow{
  sample_ids = params.sample_ids?.tokenize(',') ?: []
  ch_reads = Channel.fromList(sample_ids)
    // create tuple of [sample_id, [Read1 files], [Read2 files]]
    .map{ id -> tuple(id,
                      file("${params.sample_dir}/${id}/*_R1_*.fastq.gz"),
                      file("${params.sample_dir}/${id}/*_R2_*.fastq.gz"),
                      )}
  ch_indexes_alevin = Channel.fromList([
    ['alevin_cdna_no_sa',
     's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-100/salmon_index/cdna_k31'],
    ['alevin_txome_no_sa',
     's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-100/salmon_index/txome_k31'],
    ['alevin_cdna_full_sa',
     's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-100/salmon_index/cdna_k31_full_sa'],
    ['alevin_txome_full_sa',
     's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-100/salmon_index/txome_k31_full_sa']
  ])

  ch_testset_alevin = ch_reads.combine(ch_indexes_alevin)

  ch_indexes_kallisto = Channel.fromList([
    ['kallisto_cdna',
     's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-100/kallisto_index/cdna_k31'],
    ['kallisto_txome',
     's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-100/kallisto_index/txome_k31']
  ])

  ch_testset_kallisto = ch_reads.combine(ch_indexes_kallisto)

  // run Alevin
  alevin(ch_testset_alevin, params.t2g)

  // run Kallisto
  kallisto_bus(ch_testset_kallisto)
  bustools_correct(kallisto_bus.out, params.barcodes)
  bustools_count(bustools_correct.out, params.t2g)
}