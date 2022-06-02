#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process mpileup{
  container params.BCFTOOLS_CONTAINER
  publishDir "${params.outdir}/mpileup/${meta.multiplex_library_id}"
  tag "${meta.multiplex_run_id}"
  cpus "2"
  memory "8.GB"
  input:
    tuple val(meta), path(bamfiles), path(bamfiles_index)
    tuple path(ref_fasta), path(ref_index)
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
      --fasta-ref ${ref_fasta} \
      ${bamfiles} \
    | bcftools call -Ou --multiallelic-caller --variants-only \
    | bcftools view -Oz --genotype ^miss \
    | bcftools reheader -s samples.txt \
    > ${mpileup_file}
    """
}

// not currently used: faster than mpileup but may be less accurate
process mpileup_cellsnp{
  container params.CELLSNP_CONTAINER
  publishDir "${params.outdir}/mpileup_cellsnp/${meta.multiplex_library_id}"
  label "cpus_8"
  input:
    tuple val(meta), path(bamfiles), path(bamfiles_index)
    tuple path(ref_fasta), path(ref_index)
  output:
    tuple val(meta), path(mpileup_file)
  script:
    mpileup_file = "${meta.multiplex_library_id}.vcf.gz"
    outdir = "${meta.multiplex_library_id}_mpileup_cellsnp"
    """
    # create sample and bam file lists
    echo "${meta.sample_ids.join('\n')}" > samples.txt
    echo "${bamfiles.join('\n')}" > bamfiles.txt
    # decompress reference
    gunzip -c ${ref_fasta} > reference.fa
    mv ${ref_index} reference.fa.fai

    # call genotypes against reference with cellsnp
    cellsnp-lite \
      --samFileList bamfiles.txt \
      --sampleList samples.txt \
      --refseq reference.fa \
      --UMItag None \
      --cellTAG None \
      --nproc ${task.cpus} \
      --outDir cellsnp \
      --minMAF=0.1 \
      --minCOUNT=20 \
      --maxDEPTH=100000 \
      --genotype \
      --gzip
    
    mv cellsnp/cellSNP.cells.vcf.gz ${mpileup_file}
    """
}

workflow pileup_multibulk{
  take:
    multiplex_ch // a channel of multiplex meta objects, with sample_ids joined by `_` in their ids
    bulk_mapped_ch // output of bulk mapping: [meta, bamfile, bamfile_index]
  
  main:
    //pull sample to front of bulk_mapped_ch
    sample_bulk_ch = bulk_mapped_ch
      .map{[it[0].sample_id] + it} 

    pileup_ch = multiplex_ch 
      .map{[it.sample_id.tokenize("_"), it]} // split out sample ids into a tuple
      .transpose() // one element per sample (meta objects repeated)
      .combine(sample_bulk_ch, by: 0) // combine by sample id
      .groupTuple(by: 1) // group by the multiplex run meta object
      .map{[
        [ // create a meta object for each group of files
          sample_ids: it[0],
          multiplex_run_id: it[1].run_id,
          multiplex_library_id: it[1].library_id,
          bulk_run_ids: it[2].collect{it.run_id},
          bulk_run_prefixes: it[2].collect{it.files_directory}
        ],
        it[3], // bamfiles
        it[4]  // bamfile indexes
      ]}
    mpileup(pileup_ch, [params.ref_fasta, params.ref_fasta_index])
  
  emit:
    mpileup.out

}
