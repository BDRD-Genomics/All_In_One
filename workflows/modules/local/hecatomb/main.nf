#!/usr/bin/ nextflow

nextflow.enable.dsl=2

/*
========================================================================================
                                    Hecatomb
========================================================================================
   Github   : 
   Contact  :     
----------------------------------------------------------------------------------------

*/

process Hecatomb {

    tag { sample_id }

    errorStrategy 'ignore'
    publishDir "${params.outdir}", mode: 'copy'
    label 'high'

    input:
    tuple val(sample_id), val(qc_SR_Read1), val(qc_SR_Read2), val(qc_long_read)

    output:
    tuple val(sample_id), val("merged_assembly.fasta"), emit: hecatomb_merged_fasta
    stdout emit: hecatomb_ch

    when:
    params.nanopore_assembly

    script:
    // Dynamically define rpath
    def rpath = null
    if (qc_SR_Read1 && qc_SR_Read2 && qc_SR_Read1 != '-' && qc_SR_Read2 != '-') {
        rpath = qc_SR_Read1.toString().replaceAll(/\/[^\/]+$/, '')  // Get directory from Read1
    } else if (qc_long_read && qc_long_read != '-') {
        rpath = qc_long_read.toString().replaceAll(/\/[^\/]+$/, '') // Get directory from long read
    }

    """
    eval "\$(command conda 'shell.bash' 'hook' 2> /dev/null)"
    conda activate hecatomb

    echo "Run hecatomb"
    mkdir -p ${params.outdir}/${params.project_id}/${sample_id}/hecatomb

    hecatomb run all --profile slurm \
        --threads ${task.cpus} \
        --reads ${rpath} \
        --output ${params.outdir}/${params.project_id}/${sample_id}/hecatomb \
        ${params.hecatomb_config}

    echo "hecatomb complete"
    """
}
