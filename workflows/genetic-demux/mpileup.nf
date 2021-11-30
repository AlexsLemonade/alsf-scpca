#!/usr/bin/env nextflow
nextflow.enable.dsl=2

SAMTOOLSCONTAINER = 'quay.io/biocontainers/samtools:1.14--hb421002_0'
BCFTOOLSCONTAINER = 'quay.io/biocontainers/bcftools:1.14--h88f3f91_0'

params.run_metafile = 's3://ccdl-scpca-data/sample_info/scpca-library-metadata.tsv'
params.run_ids = 'SCPCR000533'
params.outdir = 's3://nextflow-ccdl-results/scpca/demux/mpileup'

params.ref_fasta = 's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-104/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz'
params.ref_fasta_index = 's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-104/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai'

// supported single cell technologies
single_cell_techs = ['10Xv2', '10Xv3', '10Xv3.1', '10Xv2_5prime']

process mpileup{
  container BCFTOOLSCONTAINER
  publishDir "${params.outdir}/${meta.multiplex_library_id}"
  cpus "2"
  memory "2.G"
  input:
    tuple val(meta), path(bamfiles), path(bamfiles_index)
    tuple path(reference), path(reference_index)
  output:
    tuple val(meta), path(mpileup_file)
  script:
    mpileup_file = "${meta.multiplex_library_id}.vcf.gz"
    """
    # create sample file to use for header
    echo "${meta.sample_ids.join('\n')}" > samples.txt
    
    # call genotypes against reference, filter sites with missing genotypes 
    # & use sample names for header (replacing file names)
    bcftools mpileup -Ou \
      --fasta-ref ${reference} \
      ${bamfiles} \
    | bcftools call -Ou --multiallelic-caller --variants-only \
    | bcftools view -Oz --genotype ^miss \
    | bcftools reheader -s samples.txt \
    > ${mpileup_file}
    """
}


// bulk mapping results as they would come from a channel
bulk_results = [
  [ [ // meta 
      run_id: 'SCPCR000170', 
      sample_id: 'SCPCS000133', 
      technology: 'single_end', 
      seq_unit: 'bulk', 
      s3_prefix: 'sccdl-scpca-data/runs/SCPCR000170'
    ],
    file('s3://nextflow-ccdl-results/scpca/star-bulk/SCPCS000133/SCPCR000170.sorted.bam'),
    file('s3://nextflow-ccdl-results/scpca/star-bulk/SCPCS000133/SCPCR000170.sorted.bam.bai')
  ],
  [ [ 
      run_id: 'SCPCR000171', 
      sample_id: 'SCPCS000134', 
      technology: 'single_end', 
      seq_unit: 'bulk', 
      s3_prefix: 'sccdl-scpca-data/runs/SCPCR000171'
    ],
    file('s3://nextflow-ccdl-results/scpca/star-bulk/SCPCS000134/SCPCR000171.sorted.bam'),
    file('s3://nextflow-ccdl-results/scpca/star-bulk/SCPCS000134/SCPCR000171.sorted.bam.bai')
  ],
  [ [
      run_id: 'SCPCR000172', 
      sample_id: 'SCPCS000135', 
      technology: 'single_end', 
      seq_unit: 'bulk', 
      s3_prefix: 'sccdl-scpca-data/runs/SCPCR000172'
    ],
    file('s3://nextflow-ccdl-results/scpca/star-bulk/SCPCS000135/SCPCR000172.sorted.bam'),
    file('s3://nextflow-ccdl-results/scpca/star-bulk/SCPCS000135/SCPCR000172.sorted.bam.bai')
  ],
  [ [
      run_id: 'SCPCR000173', 
      sample_id: 'SCPCS000136', 
      technology: 'single_end', 
      seq_unit: 'bulk', 
      s3_prefix: 'sccdl-scpca-data/runs/SCPCR000173'
    ],
    file('s3://nextflow-ccdl-results/scpca/star-bulk/SCPCS000136/SCPCR000173.sorted.bam'),
    file('s3://nextflow-ccdl-results/scpca/star-bulk/SCPCS000136/SCPCR000173.sorted.bam.bai')
  ]
]

workflow{
  run_ids = params.run_ids?.tokenize(',') ?: []
  run_all = run_ids[0] == "All"

  all_ch = Channel.fromPath(params.run_metafile)
    .splitCsv(header: true, sep: '\t')
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
    .filter{it.run_id in run_ids}

  // only multiplexed samples as [sample_id, run_id] pairs
  multiplexed_ch = all_ch.filter{it.technology in single_cell_techs}
    .filter{it.sample_id.contains("_")} 
    .map{[tuple(it.sample_id.split("_")), it]} // split out sample ids into a tuple
    .transpose() // one element per sample (meta objects repeated)

  bulk_ch = Channel.from(bulk_results)
    // pull sample id out, leave other three fields as is
    .map{[it[0].sample_id, it[0], it[1], it[2]]}
  
  pileup_ch = multiplexed_ch.combine(bulk_ch, by: 0) // combine by sample id
    .groupTuple(by: 1) // group by the multiplex run meta object
    .map{[
      [ // create a meta object for each group of files
        sample_ids: it[0],
        multiplex_run_id: it[1].run_id,
        multiplex_library_id: it[1].library_id,
        bulk_run_ids: it[2].collect{it.run_id},
        bulk_run_prefixes: it[2].collect{it.s3_prefix}
      ],
      it[3], // bamfiles
      it[4]  // bamfile indexes
     ]}
  
  mpileup(pileup_ch, [params.ref_fasta, params.ref_fasta_index])
}
