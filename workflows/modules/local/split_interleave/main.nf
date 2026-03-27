#!/usr/bin/ nextflow

nextflow.enable.dsl=2

process Split_Interleave {
    tag { sample_id }
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/split_interleaved_host_removed/", mode: 'copy'
    label 'optimized_split_interleaved'

    input:
    tuple val(sample_id), path(interleaved_fastq), val(mode)

    output:
    tuple val(sample_id), path("${sample_id}_host_removed_sr_R1.fastq"), path("${sample_id}_host_removed_sr_R2.fastq"), emit: split_interleave_ch

    when:
    mode in ['short', 'hybrid']

    script:
    """
    eval "\$(command conda 'shell.bash' 'hook' 2>/dev/null)"
    conda activate bbmap
    reformat.sh in=${interleaved_fastq} out1=${sample_id}_host_removed_sr_R1.fastq out2=${sample_id}_host_removed_sr_R2.fastq
    """
} 
