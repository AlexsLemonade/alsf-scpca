#!/usr/bin/env nextflow
nextflow.enable.dsl=2

STARCONTAINER = 'quay.io/biocontainers/star:2.7.9a--h9ee0642_0'
SAMTOOLSCONTAINER = 'quay.io/biocontainers/samtools:1.14--hb421002_0'

// parameters
params.ref  = 'Homo_sapiens.GRCh38.104'
params.star_index = 's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-104/star_index/Homo_sapiens.GRCh38.104.star_idx'
params.barcode_dir = 's3://nextflow-ccdl-data/reference/10X/barcodes' 

params.run_metafile = 's3://ccdl-scpca-data/sample_info/scpca-library-metadata.tsv'
params.run_ids = 'SCPCR000533'

params.outdir = 's3://nextflow-ccdl-results/scpca/starsolo/'

// 10X barcode files
cell_barcodes = ['10Xv2': '737K-august-2016.txt',
                 '10Xv3': '3M-february-2018.txt',
                 '10Xv3.1': '3M-february-2018.txt',
                 '10Xv2_5prime': '737K-august-2016.txt']

// supported single cell technologies
single_cell_techs = cell_barcodes.keySet()

process starsolo{
  container STARCONTAINER
  label 'bigdisk'
  memory "32.GB"
  cpus "8"
  publishDir "${params.outdir}/${meta.library_id}"
  input:
    tuple val(meta), path(read1), path(read2)
    path star_index
    path barcode_file
  output:
    tuple val(meta), path(output_dir)
  script:
    tech_flag = ['10Xv2': '',
                 '10Xv2_5prime': '',
                 '10Xv3': '--soloUMIlen 12',
                 '10Xv3.1': '--soloUMIlen 12']
    output_dir = "${meta.run_id}"
    """
    mkdir -p ${output_dir}/Solo.out/Gene/raw
    STAR \
      --soloType CB_UMI_Simple \
      --genomeDir ${star_index} \
      --runThreadN ${task.cpus} \
      --readFilesIn ${read2.join(',')} ${read1.join(',')} \
      --readFilesCommand gunzip -c \
      --soloCBwhitelist ${barcode_file} \
      ${tech_flag[meta.technology]} \
      --soloCellFilter EmptyDrops_CR \
      --outSAMtype BAM Unsorted \
      --outSAMattributes NH HI nM AS CR UR CB UB CY UY GX GN \
      --runDirPerm All_RWX \
      --outFileNamePrefix ${output_dir}/ 
    """
}

process index_bam{
  container SAMTOOLSCONTAINER
  
  input:
    tuple val(meta), path(bamfile)
  output:
    tuple val(meta), path(bamfile), path(bamfile_index)
  script:
    bamfile_index = "${bamfile}.bai"
    """
    samtools index ${bamfile} ${bamfile_index}
    """
}

workflow{
  sc_techs = ['single_end', 'paired_end']
  run_ids = params.run_ids?.tokenize(',') ?: []
  run_all = run_ids[0] == "All"

  singlecell_ch = Channel.fromPath(params.run_metafile)
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
      feature_barcode_file: it.feature_barcode_file,
      feature_barcode_geom: it.feature_barcode_geom,
      s3_prefix: it.s3_prefix,
    ]}
    // only bulk samples
    .filter{it.technology in single_cell_techs} 
    // use only the rows in the run_id list (run, library, or sample can match)
    // or run by project or submitter if the project parameter is set
    .filter{run_all 
             || (it.run_id in run_ids) 
            }
    
    sc_reads_ch = singlecell_ch
      .map{meta -> tuple(meta,
                         file("s3://${meta.s3_prefix}/*_R1_*.fastq.gz"),
                         file("s3://${meta.s3_prefix}/*_R2_*.fastq.gz"))}
    
    cellbarcodes_ch = singlecell_ch
      .map{file("${params.barcode_dir}/${cell_barcodes[it.technology]}")}
    
    starsolo(sc_reads_ch, params.star_index, cellbarcodes_ch)
}
