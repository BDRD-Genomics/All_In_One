#!/usr/bin/ nextflow

nextflow.enable.dsl=2

/*
========================================================================================
   CheckM Workflow
========================================================================================
   Github   : 
   Contact  :     
----------------------------------------------------------------------------------------

*/

process CheckM_Assemblies {
    tag { "${assembler} | ${sample_id}" }
    publishDir path: { "${params.outdir}/${params.project_id}/${sample_id}/checkm/assemblies/${assembler}/" }, mode: 'copy'
    label 'normal'
    errorStrategy 'ignore'
    cpus { 32 }
    memory { '64 GB'}
    time '36h'
    beforeScript 'export CHECKM_DATA_DIR=/export/database/checkm'
    conda "${baseDir}/env/checkm-genome.yml"

    when:
    params.run_checkm
    input:
    tuple val(sample_id), val(assembler), path(contigs_fasta)

    output:
    tuple val(sample_id), val(assembler), file("checkm_summary.tsv"), emit: checkm_summary_ch
    tuple val(sample_id), val(assembler), file("checkm_full_stats.tsv"), emit: checkm_full_stats_ch
    tuple val(sample_id), val(assembler), file("lineage.ms"), emit: checkm_lineage_ch
    tuple val(sample_id), val(assembler), file("*"), emit: checkm_all_files_ch

    script:
    """
    set -euo pipefail
    #eval "\$(command conda 'shell.bash' 'hook' 2> /dev/null)"
    #conda activate checkm-genome
    # Normalize to .fasta so we can use one -x flag regardless of source
    ln -sf "${contigs_fasta}" "${sample_id}.fasta"
    checkm lineage_wf -t ${task.cpus} -x fasta --tab_table -f checkm_summary.tsv . . --tmpdir /tmp
    checkm qa -t ${task.cpus} --tab_table -f checkm_full_stats.tsv -o 2 lineage.ms .
    """
}
