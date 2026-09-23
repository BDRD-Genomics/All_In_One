#!/usr/bin/ nextflow

nextflow.enable.dsl=2

/*
========================================================================================
   Dorado Polish Workflow
========================================================================================
   Github   : 
   Contact  :     
----------------------------------------------------------------------------------------

*/

process Dorado_Split_BamFiles {

    tag "bamfile"
    publishDir "${params.outdir}/${params.project_id}/bam_files/", mode: 'copy'
    label 'normal'

    input:
    path(bamfile)

    output:
    path("*.bam")

    script:
    """
    ${params.scripts}/split_rename_bams.sh $bamfile
    """
}

process Dorado_Aligner {
    tag { sample_id }
    label 'normal'

    // publish ONLY the bam_pass/<sample_id>/ files
    publishDir { "${params.outdir}/${params.project_id}/${sample_id}/dorado_aligner" },
        mode: 'copy',
        pattern: { "unknown/**/bam_pass/${sample_id}/*" }

    // publish the summary too (it lives in --output-dir)
    publishDir { "${params.outdir}/${params.project_id}/${sample_id}/dorado_aligner" },
        mode: 'copy',
        pattern: "alignment_summary.txt"

    input:
    tuple val(sample_id), path(bam_file), path(contigs_fasta)

    output:
    // keep your old 3-tuple stable
    tuple val(sample_id),
          path("unknown/**/bam_pass/${sample_id}/*"), emit: dorado_aligner_ch

    // emit summary separately to avoid tuple-shape breakage
    tuple val(sample_id),
          path("alignment_summary.txt"), emit: dorado_aligner_summary_ch

    script:
    """
    ${params.dorado_software}/dorado aligner ${contigs_fasta} ${bam_file} \
      --output-dir . --emit-summary -t 32
    """
}

process Dorado_Polisher {
    tag { sample_id }
    publishDir { "${params.outdir}/${params.project_id}/${sample_id}/dorado_polish/" }, mode: 'copy'
    clusterOptions '--partition="normal" --gpus=2'
    errorStrategy 'ignore'
    input:
    tuple val(sample_id), path(bam_file), path(bam_index), path(contigs_fasta)

    output:
    tuple val(sample_id), file("*"), emit: dorado_polisher_ch

    script:
    """
    ${params.dorado_software}/dorado polish \\
        --device auto \\
        --bacteria \\
        --models-directory ${params.dorado_polish_model} \\
        -o . \\
        ${bam_file} \\
        ${contigs_fasta}
    """
}
