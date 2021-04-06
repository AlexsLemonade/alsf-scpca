#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// run parameters
params.ref_dir = 's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-100'
params.index_dir = 'kallisto_index'
params.index_name = 'txome_k31'
params.annotation_dir = 'annotation'
params.t2g = 'Homo_sapiens.ensembl.100.tx2gene.tsv'
params.mitolist = 'Homo_sapiens.ensembl.100.mitogenes.txt'
params.barcode_dir = 's3://nextflow-ccdl-data/reference/10X/barcodes' 
// 10X barcode files
barcodes = ['10Xv2': '737K-august-2016.txt',
            '10Xv3': '3M-february-2018.txt',
            '10Xv3.1': '3M-february-2018.txt']

params.run_metafile = 's3://ccdl-scpca-data/sample_info/scpca-library-metadata.tsv'
params.run_ids = "SCPCR000003" //comma separated list to be parsed into a list
params.outdir = 's3://nextflow-ccdl-results/scpca/kallisto-quant'

// build full paths
params.index_path = "${params.ref_dir}/${params.index_dir}/${params.index_name}"
params.t2g_path = "${params.ref_dir}/${params.annotation_dir}/${params.t2g}"
params.mito_path = "${params.ref_dir}/${params.annotation_dir}/${params.mitolist}"


process kallisto_bus{
  container 'quay.io/biocontainers/kallisto:0.46.2--h4f7b962_1'
  label 'cpus_8'
  tag "${id}-${index}"
  input:
    tuple val(id), val(tech), path(read1), path(read2)
    path index
  output:
    path run_dir
  script:
    run_dir = "${id}-${index}"
    // interleave read1 & read2 file paths
    reads = [read1, read2].transpose().flatten().join(' ')
    // translate for kallisto technology flag
    tech_x = ['10Xv2': '10xv2',
              '10Xv3': '10xv3',
              '10Xv3.1': '10xv3']
    """
    mkdir -p ${run_dir}/bus
    kallisto bus \
      -i ${index} \
      -o ${run_dir}/bus \
      -x ${tech_x[tech]} \
      -t ${task.cpus} \
      ${reads}
    """
}

// process to create a barcode whitelist if not provided
// not currently used
process bustools_whitelist{
  container 'quay.io/biocontainers/bustools:0.40.0--h4f7b962_0'
  label 'cpus_8'
  input:
    path run_dir
  output:
    path barcode_file
  script:
    barcode_file = "${run_dir}/bus/barcodes.txt"
    """
    bustools sort \
      -o output.sorted.bus \
      -t ${task.cpus} \
      ${run_dir}/bus/output.bus
    bustools whitelist \
      -o ${barcode_file} \
      output.sorted.bus
    """
}

process bustools_correct{
  container 'quay.io/biocontainers/bustools:0.40.0--h4f7b962_0'
  label 'cpus_8'
  publishDir "${params.outdir}"
  input:
    path run_dir
    path barcode_file
  output:
    path run_dir
  script:
    """
    bustools correct \
      -p \
      -w ${barcode_file} \
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
  publishDir "${params.outdir}"
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
  run_ids = params.run_ids?.tokenize(',') ?: []
  run_all = run_ids[0] == "All"
  samples_ch = Channel.fromPath(params.run_metafile)
    .splitCsv(header: true, sep: '\t')
    .filter{it.technology in barcodes} // only 10X data
    // use only the rows in the sample list
    .filter{run_all || (it.scpca_run_id in run_ids)}
  // create tuple of [sample_id, technology, [Read1 files], [Read2 files]]
  reads_ch = samples_ch
    .map{row -> tuple(row.scpca_run_id,
                      row.technology,
                      file("s3://${row.s3_prefix}/*_R1_*.fastq.gz"),
                      file("s3://${row.s3_prefix}/*_R2_*.fastq.gz"),
                      )}
  barcodes_ch = samples_ch
    .map{row -> file("${params.barcode_dir}/${barcodes[row.technology]}")}
  // run kallisto
  kallisto_bus(reads_ch, params.index_path)
  // correct busfiles using expected barcodes
  bustools_correct(kallisto_bus.out, barcodes_ch)
  // count genes
  bustools_count(bustools_correct.out, params.t2g_path)
}
