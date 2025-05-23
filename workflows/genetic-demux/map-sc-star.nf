#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process starsolo{
  container params.STAR_CONTAINER
  publishDir "${params.outdir}/starsolo/${meta.library_id}"
  tag "${meta.run_id}"
  label 'bigdisk'
  memory "32.GB"
  cpus "8"
  input:
    tuple val(meta), path(read1), path(read2)
    path star_index
    path barcode_file
  output:
    tuple val(meta), path(output_dir), emit: starsolo_dir
    tuple val(meta), path(output_bam), emit: star_bam
  script:
    tech_flag = ['10Xv2': '',
                 '10Xv2_5prime': '',
                 '10Xv3': '--soloUMIlen 12',
                 '10Xv3.1': '--soloUMIlen 12']
    features_flag = meta.seq_unit == "nucleus" ? "--soloFeatures Gene GeneFull" : "--soloFeatures Gene"
    output_dir = "${meta.run_id}_star"
    output_bam = "${meta.run_id}.sorted.bam"
    
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
      ${features_flag} \
      --soloCellFilter EmptyDrops_CR \
      --outSAMtype BAM SortedByCoordinate \
      --outSAMattributes NH HI nM AS CR UR CB UB CY UY GX GN \
      --outBAMsortingThreadN 2 \
      --limitBAMsortRAM 20000000000 \
      --runDirPerm All_RWX \
      --outFileNamePrefix ${output_dir}/

    mv ${output_dir}/Aligned.sortedByCoord.out.bam ${output_bam}
    """
}

process index_bam{
  container params.SAMTOOLS_CONTAINER
  label 'bigdisk'
  publishDir "${params.outdir}/starsolo/${meta.library_id}"
  tag "${meta.run_id}"
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

workflow starsolo_sc{
  take:
    singlecell_ch

  main: 
    sc_reads_ch = singlecell_ch
      .map{meta -> tuple(meta,
                         file("${meta.files_directory}/*_R1_*.fastq.gz"),
                         file("${meta.files_directory}/*_R2_*.fastq.gz"))}

    cellbarcodes_ch = singlecell_ch
      .map{file("${params.barcode_dir}/${params.cell_barcodes[it.technology]}")}

    starsolo(sc_reads_ch, params.star_index, cellbarcodes_ch)
    index_bam(starsolo.out.star_bam)
  
  emit:
    bam = index_bam.out
    quant = starsolo.out.starsolo_dir
}
