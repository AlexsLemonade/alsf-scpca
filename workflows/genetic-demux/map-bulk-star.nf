#!/usr/bin/env nextflow
nextflow.enable.dsl=2

STARCONTAINER = 'quay.io/biocontainers/star:2.7.9a--h9ee0642_0'
SAMTOOLSCONTAINER = 'quay.io/biocontainers/samtools:1.14--hb421002_0'

// parameters
params.ref  = 'Homo_sapiens.GRCh38.104'
params.star_index = 's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-104/star_index/Homo_sapiens.GRCh38.104.star_idx'

params.run_metafile = 's3://ccdl-scpca-data/sample_info/scpca-library-metadata.tsv'
params.run_ids = 'SCPCR000170,SCPCR000171,SCPCR000172,SCPCR000173'

params.outdir = 's3://nextflow-ccdl-results/scpca/demux/star-bulk/'



process bulkmap_star{
  container STARCONTAINER
  memory "32.GB"
  cpus "8"
  input:
    tuple val(meta), path(read1), path(read2)
    path star_index
  output:
    tuple val(meta), path(output_bam)
  script:
    output_bam = "${meta.run_id}.sorted.bam"
    """
    STAR \
      --genomeDir ${star_index} \
      --runThreadN ${task.cpus} \
      --readFilesIn ${read1.join(',')} \
      ${meta.technology == 'paired_end' ? read2.join(',') : ""} \
      --readFilesCommand gunzip -c \
      --outSAMtype BAM SortedByCoordinate

    mv Aligned.sortedByCoord.out.bam ${output_bam}
    """
}

process index_bam{
  container SAMTOOLSCONTAINER
  publishDir "${params.outdir}/${meta.sample_id}"
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
      s3_prefix: it.s3_prefix,
    ]}
    // only bulk samples
    .filter{it.technology in bulk_techs} 
    // use only the rows in the run_id list (run, library, or sample can match)
    .filter{run_all 
             || (it.run_id in run_ids) 
            }
    
    bulk_reads_ch = bulk_ch
      .map{meta -> tuple(meta,
                         file("s3://${meta.s3_prefix}/*_R1_*.fastq.gz"),
                         file("s3://${meta.s3_prefix}/*_R2_*.fastq.gz"))}
    
    bulkmap_star(bulk_reads_ch, params.star_index) \
      | index_bam
}
