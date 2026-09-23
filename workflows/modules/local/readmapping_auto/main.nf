#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
   Read Mapping Auto

Pulls a local reference FASTA for every kraken2/bracken/mash/BLAST(contigs) hit, optionally
dereplicates near-identical references by ANI, then runs CoverM (wraps minimap2).
*/

process GET_REFERENCES {
    tag "${sample_id}"
    label 'readmapping_auto'
    conda "${baseDir}/env/readmapping_auto.yml"
    publishDir { "${params.outdir}/${params.project_id}/${sample_id}/readmapping/references" }, mode: 'copy'
    errorStrategy 'ignore'
    cpus 64

    input:
    tuple val(sample_id), path(kraken_report), path(bracken_table), path(mash_screen), path(blast_tsvs)

    output:
    tuple val(sample_id), path("${sample_id}.references.fasta"), path("${sample_id}.references.tsv"), path("genomes"), emit: references_ch
    path "${sample_id}.not_found.tsv", emit: not_found_ch

    script:
    """
    export BLASTDB_LMDB_MAP_SIZE=100000000

    python3 ${params.scripts}/readmapping_auto_get_references.py \\
        --sample-id ${sample_id} \\
        --kraken-report ${kraken_report} \\
        --bracken-table ${bracken_table} \\
        --mash-screen ${mash_screen} \\
        --blast-tsv ${blast_tsvs} \\
        --blast-db ${params.blast_core_nt} \\
        --assembly-summary ${params.assembly_summary_refseq} \\
        --assembly-summary-historical ${params.assembly_summary_refseq_historical} \\
        --gcf-genomes-dir ${params.readmapping_auto_gcf_genomes_dir} \\
        --refseq-prok-dir ${params.readmapping_auto_refseq_prok_dir} \\
        --refseq-viruses-dir ${params.readmapping_auto_refseq_viruses_dir} \\
        --threads 12 \\
        --outdir .
    """
}

process DEDUP_ANI {
    tag "${sample_id}"
    label 'readmapping_auto'
    conda "${baseDir}/env/readmapping_auto.yml"
    publishDir { "${params.outdir}/${params.project_id}/${sample_id}/readmapping/references" }, mode: 'copy'
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), path(references_fasta), path(references_tsv), path(genomes)

    output:
    tuple val(sample_id), path("${sample_id}.references_dedup.fasta"), emit: dedup_fasta_ch
    path "${sample_id}.cluster_manifest.tsv", emit: cluster_manifest_ch

    script:
    """
    python3 ${params.scripts}/readmapping_auto_dedup_ani.py \\
        --sample-id ${sample_id} \\
        --references-tsv ${references_tsv} \\
        --ani-threshold ${params.readmapping_auto_ani_threshold} \\
        --threads ${task.cpus} \\
        --outdir .
    """
}

process COVERM {
    tag "${sample_id}"
    label 'readmapping_auto'
    conda "${baseDir}/env/readmapping_auto.yml"
    publishDir { "${params.outdir}/${params.project_id}/${sample_id}/readmapping" }, mode: 'copy'
    errorStrategy 'ignore'
    cpus 64

    input:
    tuple val(sample_id), path(reference_fasta), path(long_read)

    output:
    tuple val(sample_id), path("bams/*.bam"), emit: bam_ch
    tuple val(sample_id), path("stats/${sample_id}.coverm.tsv"), emit: stats_ch

    script:
    """
    mkdir -p bams stats
    minimap2 -ax map-ont -t ${task.cpus} -I 50g "${reference_fasta}" "${long_read}" | samtools sort -@ ${task.cpus} -o bams/${sample_id}_sorted.bam
    samtools index bams/${sample_id}_sorted.bam

    coverm contig --bam-files bams/${sample_id}_sorted.bam -m mean covered_bases covered_fraction length count -t ${task.cpus} -o stats/${sample_id}.coverm.tsv
    """
}
