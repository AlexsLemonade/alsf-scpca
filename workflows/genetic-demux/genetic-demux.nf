#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.run_metafile = 's3://ccdl-scpca-data/sample_info/scpca-library-metadata.tsv'
params.run_ids = 'SCPCR000533'
params.outdir = 's3://nextflow-ccdl-results/scpca/demux/cellsnp'
params.barcode_dir = 's3://nextflow-ccdl-data/reference/10X/barcodes'

params.ref_fasta = 's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-104/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz'
params.ref_fasta_index = 's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-104/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai'
params.star_index = 's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-104/star_index/Homo_sapiens.GRCh38.104.star_idx'


// 10X barcode files
cell_barcodes = [
  '10Xv2': '737K-august-2016.txt',
  '10Xv2_5prime': '737K-august-2016.txt',
  '10Xv3': '3M-february-2018.txt',
  '10Xv3.1': '3M-february-2018.txt'
  ]

// supported technologies
single_cell_techs= cell_barcodes.keySet()
bulk_techs = ['single_end', 'paired_end']
all_techs = single_cell_techs + bulk_techs


// include processes
include { bulkmap_star; index_bam } from './map-bulk-star.nf'
include { starsolo } from './map-sc-star.nf'
include { mpileup } from './mpileup.nf'
include { cellsnp; vireo } from './cellsnp.nf'

workflow{
  run_ids = params.run_ids?.tokenize(',') ?: []
  run_all = run_ids[0] == "All"

  runs_ch = Channel.fromPath(params.run_metafile)
    .splitCsv(header: true, sep: '\t')
    // convert row data to a metadata map, keeping only columns we will need (& some renaming)
    // sample_id can have multiple `;`-separated values, change to `_`
    .map{[
      run_id: it.scpca_run_id,
      library_id: it.scpca_library_id,
      sample_id: it.scpca_sample_id.replaceAll(";", "_"),
      project_id: it.scpca_project_id?: "no_project",
      submitter: it.submitter,
      technology: it.technology,
      seq_unit: it.seq_unit,
      s3_prefix: it.s3_prefix
    ]}   
    // only technologies we know how to process
    .filter{it.technology in all_techs} 

  multiplex_ch = runs_ch
    .filter{it.technology in single_cell_techs}
    .filter{run_all 
             || (it.run_id in run_ids) 
             || (it.library_id in run_ids)
            }
    .filter{it.sample_id.contains("_")} 
    .map{[it.sample_id.tokenize("_"), it]} // split out sample ids into a tuple
    .transpose() // one element per sample (meta objects repeated)
  
  bulk_ch = runs_ch
    .filter {it.technology in bulk_techs}
    .map{[it.sample_id, it]}
  
  multiplex_bulk_ch = multiplex_ch
    .combine(bulk_ch, by:0)
    .view()
}
