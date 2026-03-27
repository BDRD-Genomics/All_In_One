#!/usr/bin/ nextflow

nextflow.enable.dsl=2

/*
========================================================================================
                                    QUAST
========================================================================================
   Github   : 
   Contact  :     
----------------------------------------------------------------------------------------

*/

process Quast {

    tag { sample_id }

    errorStrategy 'ignore'
    publishDir "${params.outdir}", mode: 'copy'
    label 'low'

    input:
    tuple val(sample_id), val(qc_SR_Read1), val(qc_SR_Read2), val(qc_long_read), val(circos_fasta)

    output:
    stdout emit: stage_7a_quast_ch

    when:
    params.nanopore_assembly

    script:
    """
    eval "\$(command conda 'shell.bash' 'hook' 2> /dev/null)"
    conda activate quast

    mkdir -p ${params.outdir}/${params.project_id}/${sample_id}/quast

    quast \\
        ${params.outdir}/${params.project_id}/${sample_id}/fastp/${sample_id}_fastp_ONT.fastq.gz \\
        ${params.outdir}/${params.project_id}/${sample_id}/fastp/${sample_id}_fastp_R1.fastq.gz \\
        ${params.outdir}/${params.project_id}/${sample_id}/fastp/${sample_id}_fastp_R2.fastq.gz \\
        --output-dir ${params.outdir}/${params.project_id}/${sample_id}/quast \\
        --threads ${params.threads} \\
        --circos ${circos_fasta}
    """
}
