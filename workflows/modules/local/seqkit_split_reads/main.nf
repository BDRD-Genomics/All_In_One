#!/usr/bin/ nextflow

nextflow.enable.dsl=2

/*
========================================================================================
   Seqkit Split Reads
========================================================================================
   Github   : 
   Contact  :     
----------------------------------------------------------------------------------------

*/
process Split_Short_Reads {

    tag { sample_id }
    publishDir { "${params.outdir}/${params.project_id}/${sample_id}/trim/quality_control/" }, mode: 'copy'
    label 'seqkit_split_reads'
    input:
    tuple val(sample_id), path(fastq)

    output:
    tuple val(sample_id), path("*.fastq.gz"), emit: split_shorts_reads_ch

    script:
    """
    seqkit split2 --threads ${task.cpus} \
	${fastq} \
	-l 435M \
        -O . \
	-f
    """
}

process Split_Long_Reads {

    tag { sample_id }
    publishDir { "${params.outdir}/${params.project_id}/${sample_id}/trim/quality_control/" }, mode: 'copy'
    label 'seqkit_split_reads'
    input:
    tuple val(sample_id), path(fastq)

    output:
    tuple val(sample_id), path("*.fastq.gz"), emit: split_long_reads_ch

    script:
    """
    seqkit split2 --threads ${task.cpus} \
        ${fastq} \
        -l 435M \
        -O . \
        -f
    """
}
