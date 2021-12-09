#!/usr/bin/env nextflow
nextflow.enable.dsl=2


params.run_metafile = 's3://ccdl-scpca-data/sample_info/scpca-library-metadata.tsv'
params.run_ids = 'SCPCR000533'
params.outdir = 's3://nextflow-ccdl-results/scpca/demux'
params.barcode_dir = 's3://nextflow-ccdl-data/reference/10X/barcodes'

params.ref_fasta = 's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-104/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz'
params.ref_fasta_index = 's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-104/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai'


process cellsnp{
  container params.CELLSNP_CONTAINER
  publishDir "${params.outdir}/cellsnp/${meta.library_id}"
  label "cpus_8"
  input:
    tuple val(meta_star), path(star_bam), path(star_bai), path(star_quant),
          val(meta_mpileup), path(vcf_file)
  output:
    tuple val(meta), path(outdir), path(vcf_file)
  script:
    meta = meta_star
    meta.sample_ids = meta_mpileup.sample_ids
    meta.bulk_run_ids = meta_mpileup.bulk_run_ids
    quant_dir = meta_star.seq_unit == "nucleus" ? "GeneFull" : "Gene"
    barcodes = "${star_quant}/Solo.out/${quant_dir}/filtered/barcodes.tsv"
    outdir = "${meta.library_id}_cellSNP"
    """
    cellsnp-lite \
      --samFile ${star_bam} \
      --barcodeFile ${barcodes} \
      --regionsVCF <(gunzip -c ${vcf_file}) \
      --nproc ${task.cpus} \
      --outDir ${outdir} \
      --minMAF=0.1 \
      --minCOUNT=20 \
      --gzip
    """
}

process vireo{
  container params.CONDA_CONTAINER
  publishDir "${params.outdir}/cellsnp/${meta.library_id}"
  label "cpus_8"
  input:
    tuple val(meta), path(cellsnp_dir), path(vcf_file)
  output:
    tuple val(meta), path(outdir)
  script:
    outdir = "${meta.library_id}_vireo"
    """
    pip install vireoSNP==0.5.6
    vireo \
      --cellData ${cellsnp_dir} \
      --donorFile ${vcf_file}  \
      --outDir ${outdir} \
      --nproc ${task.cpus}
    """
}

workflow cellsnp_vireo {
  take: 
    starsolo_bam_ch //channel of [meta, bamfile, bam.bai]
    starsolo_quant_ch //channel of [meta, starsolo_dir]
    mpileup_vcf_ch // channel of [meta, vcf_file]
  main:
    mpileup_ch = mpileup_vcf_ch
      .map{[it[0].multiplex_library_id, it[0], it[1]]} // pull out library id for combining
    star_mpileup_ch = starsolo_bam_ch
      .combine(starsolo_quant_ch, by: 0)
      .map{[it[0].library_id, it[0], it[1], it[2]]} // add library id at start
      .combine(mpileup_ch, by: 0)
      .map{[it[1], it[2], it[3], it[4], it[5], it[6]]} // drop library id
    
    cellsnp(star_mpileup_ch)
    vireo(cellsnp.out)
  emit:
    vireo.out
}
