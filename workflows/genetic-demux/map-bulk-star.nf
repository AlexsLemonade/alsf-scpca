#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// parameters
params.ref  = 'Homo_sapiens.GRCh38.104'
params.star_index = 's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-104/star_index/Homo_sapiens.GRCh38.104.star_idx'

params.run_metafile = 's3://ccdl-scpca-data/sample_info/scpca-library-metadata.tsv'
params.run_ids = 'SCPCR000166'

params.outdir = 's3://nextflow-ccdl-results/scpca/star-bulk/'

STARCONTAINER = 'quay.io/biocontainers/star:2.7.9a--h9ee0642_0'

process bulkmap_star{
  container STARCONTAINER
  memory "32.GB"
  cpus "8"
  publishDir "${params.outdir}"
  input:
    tuple val(meta), path(read1), path(read2)
    path star_index
  output:
    tuple val(meta), path(output_dir)
  script:
    output_dir = "${meta.library_id}"
    """
    STAR \
      --genomeDir ${star_index} \
      --runThreadN ${task.cpus} \
      --readFilesIn ${read1.join(',')} \
      ${meta.technology == 'paired_end' ? read2.join(',') : ""} \
      --readFilesCommand gunzip -c \
      --outFileNamePrefix ${output_dir}/ \
      --outSAMtype BAM SortedByCoordinate 
    """
}

workflow{
  bulk_techs = ['single_end', 'paired_end']
  run_ids = params.run_ids?.tokenize(',') ?: []
  run_all = run_ids[0] == "All"

  bulk_ch = Channel.fromPath(params.run_metafile)
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
    .filter{it.technology in bulk_techs} 
    // use only the rows in the run_id list (run, library, or sample can match)
    // or run by project or submitter if the project parameter is set
    .filter{run_all 
             || (it.run_id in run_ids) 
            }
    
    bulk_reads_ch = bulk_ch
      .map{meta -> tuple(meta,
                         file("s3://${meta.s3_prefix}/*_R1_*.fastq.gz"),
                         file("s3://${meta.s3_prefix}/*_R2_*.fastq.gz"))}
    
    bulkmap_star(bulk_reads_ch, params.star_index)
}
