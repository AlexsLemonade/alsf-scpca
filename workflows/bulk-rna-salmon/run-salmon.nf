#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// run parameters
params.ref_dir = 's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-103'
params.index_dir = 'salmon_index'
params.index_name = 'spliced_txome_k31'
params.annotation_dir = 'annotation'
params.t2g = 'Homo_sapiens.GRCh38.103.spliced.tx2gene.tsv'
params.mitolist = 'Homo_sapiens.GRCh38.103.mitogenes.txt'

params.run_ids = "SARC0067"

params.outdir = 's3://soragni-sarcomas'

// build full paths
params.index_path = "${params.ref_dir}/${params.index_dir}/${params.index_name}"
params.t2g_path = "${params.ref_dir}/${params.annotation_dir}/${params.t2g}"
params.mito_path = "${params.ref_dir}/${params.annotation_dir}/${params.mitolist}"
params.fastq = "${params.outdir}/fastq_raw"

process fastqc{
    container
    label 'cpus_8'
    tag "${id}"
    publishDir "${params.outdir}"
    input:
        path fastq
    output: 
        path fastqc_reports 
    script: 
        fastqc_reports = "${params.outdir}/reports/fastqc"
        """
        mkdir -p reports/fastqc 
        fastqc ${fastq} \
        -o ${fastqc_reports}
        """

}

process fastp{
    container
    label 'cpus_8'
    tag "${id}"
    publishDir "${params.outdir}"
    input: 
        path fastq
        val id
    output: 
        path trimmed_data
        path fastp_reports
    script: 
        trimmed_data = "${params.outdir}/data/trimmed"
        fastp_reports = "${params.outdir}/reports/fastp"
        """
        mkdir -p data/trimmed
        mkdir -p reports/fastp

        fastp --in1 {input.r1} --in2 {input.r2} /
        --out1 "${trimmed_data}/${id}_R1_001_fastq.gz" /
        --out2 "${trimmed_data}/${id}_R2_001_fastq.gz" /
        --html "${fastp_reports}/${id}_fastp.html" /
        --json "${fastp_reports}/${id}_fastp.json" /
        --trim_poly_g /
        --report_title '${id} report'
        """

}

process salmon{
    container 'quay.io/biocontainers/salmon:1.4.0--hf69c8f4_0'
    label 'cpus_8'
    tag "${id}-${index}"
    publishDir "${params.outdir}"
    input: 
        path index
        path trimmed_data
    output: 
        path salmon_results
    script:
        salmon_results = "${params.outdir}/salmon-quant"
        """
        salmon quant -i ${index} /
        -l A /
        -1 "${trimmed_data}/${id}_R1_001_fastq.gz" /
        -2 "${trimmed_data}/${id}_R2_001_fastq.gz" /
        -o ${salmon_results} /
        --threads ${task.cpus}
        """

}

process multiqc{
    container
    label 'cpus_8'
    tag "${id}"
    publishDir "${params.outdir}"
    input:
        path fastqc_reports
        path fastp_reports
        path salmon_results
    output: 
        path multiqc_reports
    script: 
        multiqc_reports = "${params.outdir}/reports/multiqc"
        """
        mkdir -p reports/multiqc
        multiqc -n ${multiqc_reports}/multiqc_report.html /
        "${fastqc_reports}/${id}_R1_001_fastqc.zip" /
        "${fastqc_reports}/${id}_R2_001_fastqc.zip" /
        "${fastp_reports}/${id}_fastp.json" /
        "${salmon_results}/${id}/aux_info/meta_info.json" /
        "${salmon_results}/${id}/libParams/flenDist.txt"
        """

    
}
