#!/usr/bin/ nextflow

nextflow.enable.dsl=2

/*
========================================================================================
   MMSEQ
========================================================================================
   Github   : 
   Contact  :     
----------------------------------------------------------------------------------------

*/

process Parse_MMSEQ_dragonflye_contigs {
    tag {sample_id}
    errorStrategy 'ignore'
    publishDir { "${params.outdir}/${params.project_id}/${sample_id}/blast/" }, mode: 'copy'
    label 'lowmem'

    input:
    tuple val(sample_id), file(mmseq_out)

    output:
    tuple val(sample_id),file("*.parsed") ,emit: mmseqs_dragonflye_contigs_ch
    tuple val(sample_id),file("*.log") ,emit: mmseqs_dragonflye_contigs_log_ch

    script:
    """
    tail -n +2 ${mmseq_out} > tmp && mv tmp ${mmseq_out}
    python3 ${params.scripts}/VS_MD_diamond_parser_linFilt_working.py -i ${mmseq_out} -t mmseqs -r allRanks -v ${params.vhunter} -n ${params.ncbi_taxa} >> Parse_MMSEQ_dragonflye_contigs.log 2>&1
    """
}

process Parse_MMSEQ_dragonflye_medaka_contigs {
    tag {sample_id}
    errorStrategy 'ignore'
    publishDir { "${params.outdir}/${params.project_id}/${sample_id}/blast/" }, mode: 'copy'
    label 'lowmem'

    input:
    tuple val(sample_id), file(mmseq_out)

    output:
    tuple val(sample_id),file("*.parsed") ,emit: mmseqs_dragonflye_medaka_contigs_ch
    tuple val(sample_id),file("*.log") ,emit: mmseqs_dragonflye_medaka_contigs_log_ch

    script:
    """
    tail -n +2 ${mmseq_out} > tmp && mv tmp ${mmseq_out}
    python3 ${params.scripts}/VS_MD_diamond_parser_linFilt_working.py -i ${mmseq_out} -t mmseqs -r allRanks -v ${params.vhunter} -n ${params.ncbi_taxa} >> Parse_MMSEQ_dragonflye_medaka_contigs.log 2>&1
    """
}

process Parse_MMSEQ_metaspades_contigs {
    tag {sample_id}
    errorStrategy 'ignore'
    publishDir { "${params.outdir}/${params.project_id}/${sample_id}/blast/" }, mode: 'copy'
    label 'lowmem'

    input:
    tuple val(sample_id), file(mmseq_out)

    output:
    tuple val(sample_id),file("*.parsed") ,emit: mmseqs_metaspades_contigs_ch
    tuple val(sample_id),file("*.log") ,emit: mmseqs_metaspades_contigs_log_ch

    script:
    """
    tail -n +2 ${mmseq_out} > tmp && mv tmp ${mmseq_out}
    python3 ${params.scripts}/VS_MD_diamond_parser_linFilt_working.py -i ${mmseq_out} -t mmseqs -r allRanks -v ${params.vhunter} -n ${params.ncbi_taxa} >> Parse_MMSEQ_metaspades_contigs.log 2>&1
    """
}

process Parse_MMSEQ_unicycler_contigs {
    tag {sample_id}
    errorStrategy 'ignore'
    publishDir { "${params.outdir}/${params.project_id}/${sample_id}/blast/" }, mode: 'copy'
    label 'lowmem'

    input:
    tuple val(sample_id), file(mmseq_out)

    output:
    tuple val(sample_id),file("*.parsed") ,emit: mmseqs_unicycler_contigs_ch
    tuple val(sample_id),file("*.log") ,emit: mmseqs_unicycler_contigs_log_ch

    script:
    """
    tail -n +2 ${mmseq_out} > tmp && mv tmp ${mmseq_out}
    python3 ${params.scripts}/VS_MD_diamond_parser_linFilt_working.py -i ${mmseq_out} -t mmseqs -r allRanks -v ${params.vhunter} -n ${params.ncbi_taxa} >> Parse_MMSEQ_unicycler_contigs.log 2>&1
    """
}
