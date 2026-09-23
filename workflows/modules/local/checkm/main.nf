#!/usr/bin/ nextflow

nextflow.enable.dsl=2

/*
========================================================================================
   CheckM Workflow
========================================================================================
*/

process CheckM_Assemblies {
    tag { "${assembler} | ${sample_id}" }
    publishDir path: { "${params.outdir}/${params.project_id}/${sample_id}/assembly_verification/${assembler}/checkm" }, mode: 'copy'
    label 'normal'
    errorStrategy 'ignore'
    cpus { 32 }
    memory { '64 GB'}
    time '36h'

    input:
    tuple val(sample_id), val(assembler), path(contigs_fasta)

    output:
    tuple val(sample_id), val(assembler), file("checkm_summary.tsv"), emit: checkm_summary_ch
    tuple val(sample_id), val(assembler), file("checkm_full_stats.tsv"), emit: checkm_full_stats_ch
    tuple val(sample_id), val(assembler), file("lineage.ms"), emit: checkm_lineage_ch
    tuple val(sample_id), val(assembler), file("*"), emit: checkm_all_files_ch

    when:
    params.run_checkm

    script:
    """
    set -euo pipefail
    export CHECKM_DATA_PATH="${params.checkm_data_dir}"
    if [[ "${contigs_fasta}" != "${sample_id}.fasta" ]]; then
        cp "${contigs_fasta}" "${sample_id}.fasta"
    fi
    checkm lineage_wf -t ${task.cpus} -x fasta --tab_table -f checkm_summary.tsv . . --tmpdir "${params.tmp_dir}"
    checkm qa -t ${task.cpus} --tab_table -f checkm_full_stats.tsv -o 2 lineage.ms .
    """
}

process CheckM2_Assemblies {
    tag { "${assembler} | ${sample_id}" }
    publishDir path: { "${params.outdir}/${params.project_id}/${sample_id}/assembly_verification/${assembler}/" }, mode: 'copy'
    label 'normal'
    errorStrategy 'ignore'
    cpus { 32 }
    memory { '64 GB'}
    time '36h'

    input:
    tuple val(sample_id), val(assembler), path(contigs_fasta)

    output:
    tuple val(sample_id), val(assembler), file("*"), emit: checkm2_all_files_ch

    when:
    params.run_checkm

    script:
    """
    set -euo pipefail
    if [[ "${contigs_fasta}" != "${sample_id}.fasta" ]]; then
        cp "${contigs_fasta}" "${sample_id}.fasta"
    fi
    checkm2 predict -t ${task.cpus} -x fasta --specific --force --input ${sample_id}.fasta --output-directory ./checkm2 --database_path ${params.checkm2_db}
    """
}

process CheckV_Assemblies {
    tag { "${assembler} | ${sample_id}" }
    publishDir path: { "${params.outdir}/${params.project_id}/${sample_id}/assembly_verification/${assembler}/checkv" }, mode: 'copy'
    label 'normal'
    errorStrategy 'ignore'
    cpus { 32 }
    memory { '64 GB'}
    time '36h'

    input:
    tuple val(sample_id), val(assembler), path(contigs_fasta)

    output:
    tuple val(sample_id), val(assembler), file("*"), emit: checkv_all_files_ch

    when:
    params.run_checkm

    script:
    """
    set -euo pipefail
    if [[ "${contigs_fasta}" != "${sample_id}.fasta" ]]; then
        cp "${contigs_fasta}" "${sample_id}.fasta"
    fi
    checkv end_to_end ${sample_id}.fasta ./ --remove_tmp -t ${task.cpus} -d ${params.checkv_db}
    """
}
