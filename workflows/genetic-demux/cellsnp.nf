#!/usr/bin/env nextflow
nextflow.enable.dsl=2

CELLSNPCONTAINER = 'quay.io/biocontainers/cellsnp-lite:1.2.2--h22771d5_0'

params.run_metafile = 's3://ccdl-scpca-data/sample_info/scpca-library-metadata.tsv'
params.run_ids = 'SCPCR000533'
params.outdir = 's3://nextflow-ccdl-results/scpca/demux/cellsnp'
params.barcode_dir = 's3://nextflow-ccdl-data/reference/10X/barcodes' 

params.ref_fasta = 's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-104/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz'
params.ref_fasta_index = 's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-104/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai'

// 10X barcode files
cell_barcodes = ['10Xv2': '737K-august-2016.txt',
                 '10Xv3': '3M-february-2018.txt',
                 '10Xv3.1': '3M-february-2018.txt',
                 '10Xv2_5prime': '737K-august-2016.txt']

// bulk mpileup results as they would come from a channel
mpileup_results = [
  [
    [ //meta
      sample_ids: ['SCPCS000133', 'SCPCS000134', 'SCPCS000135', 'SCPCS000136'],
      multiplex_run_id: 'SCPCR000533',
      multiplex_library_id: 'SCPCL000533',
      bulk_run_ids: ['SCPCR000170', 'SCPCR000171', 'SCPCR000172', 'SCPCR000173']
    ],
    file('s3://nextflow-ccdl-results/scpca/demux/mpileup/SCPCL000533/SCPCL000533.vcf.gz')
  ]
]
// starsolo results as they would come
starsolo_results = [
  [
    [
      run_id: 'SCPCR000533', 
      sample_id: 'SCPCS000133_SCPCS000134_SCPCS000135_SCPCS000136', 
      technology: '10Xv3.1', 
      seq_unit: 'nucleus', 
      s3_prefix: 'sccdl-scpca-data/runs/SCPCR000533'
    ],
    file('s3://nextflow-ccdl-results/scpca/demux/starsolo/SCPCL000533/SCPCR000533.sorted.bam'),
    file('s3://nextflow-ccdl-results/scpca/demux/starsolo/SCPCL000533/SCPCR000533.sorted.bam.bai')
  ]
]


process cellsnp{
  container CELLSNPCONTAINER
  publishDir "${params.outdir}/${meta.library_id}"
  memory "32.GB"
  cpus "8"
  input:
    tuple val(meta_star), path(star_bam), path(star_bai)
    tuple val(meta_mpileup), path(vcf_file)
    path barcode_file
  output:
    tuple val(meta), path(outdir)
  script:
    meta = meta_star 
    meta.sample_ids = meta_mpileup.sample_ids
    meta.bulk_run_ids = meta_mpileup.bulk_run_ids
    """
    cellsnp-lite \
      --samFile ${star_bam} \
      --barcodeFile ${barcode_file} \
      --regionsVCF <(gunzip -c ${vcf_file}) \
      --nproc ${task.cpus} \
      --outDir ${meta.library_id}_cellSNP \
      --minMAF=0.1 \
      --minCOUNT=20 \
      --gzip
    """
}

workflow{
  // only multiplexed samples as [sample_id, run_id] pairs
  multiplexed_ch = Channel.from(starsolo_results)
  mpileup_ch = Channel.from(mpileup_results)
  cellbarcodes_ch = multiplexed_ch
      .map{file("${params.barcode_dir}/${cell_barcodes[it[0].technology]}")}
  cellsnp(multiplexed_ch, mpileup_ch, cellbarcodes_ch)
}
