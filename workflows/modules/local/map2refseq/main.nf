#!/usr/bin/ nextflow

nextflow.enable.dsl=2

/*
========================================================================================
   Map 2 Reference
========================================================================================

*/

/* map reads to reference sequence and use for downstream analysis (Optional)*/

process Map_Reads_2_RefSeq {
    tag { sample_id }
    publishDir { "${params.outdir}/${params.run_id}/${sample_id}/trim/quality_control/targeted_read_mapping/" }, mode: 'copy'
    errorStrategy 'ignore'
    cpus { cpus }
    memory { mem }

    input:
    tuple val(sample_id), file(fastq_1), file(fastq_2), file(long_read), val(mode), val(cpus), val(mem)

    output:
    tuple val(sample_id), file("${sample_id}_host_removed_sr.fastq.gz"),    optional: true, emit: host_removed_fq_ch_sr
    tuple val(sample_id), file("${sample_id}_host_removed_LR.fastq.gz"),    optional: true, emit: host_removed_fq_ch_lr
    tuple val(sample_id), file("${sample_id}_host_LR.fastq.gz"),            optional: true, emit: host_fq_ch_lr
    tuple val(sample_id), file("${sample_id}_host_sr.fastq.gz"),            optional: true, emit: host_fq_ch_sr
    tuple val(sample_id), file("${sample_id}_host_removed_sr_R1.fastq.gz"), optional: true, emit: host_removed_fq_ch_r1
    tuple val(sample_id), file("${sample_id}_host_removed_sr_R2.fastq.gz"), optional: true, emit: host_removed_fq_ch_r2

    when:
    params.map2reads

    script:
    """


    if [[ "${mode}" == "long" || "${mode}" == "hybrid" ]]; then
        minimap2 -ax map-ont --secondary=no -t $task.cpus --split-prefix tmpfile ${params.host_db} \\
            ${long_read} | samtools sort -@ $task.cpus -o ${sample_id}_host_removed_LR.bam --write-index -
        samtools fastq -f 4 ${sample_id}_host_removed_LR.bam | pigz > ${sample_id}_host_removed_LR.fastq.gz
        samtools fastq -F 4 ${sample_id}_host_removed_LR.bam | pigz > ${sample_id}_host_LR.fastq.gz
    fi

    if [[ "${mode}" == "short" || "${mode}" == "hybrid" ]]; then
        reformat.sh in1=${fastq_1} in2=${fastq_2} out=${sample_id}_IR.fastq.gz ow=t
        minimap2 -ax sr --secondary=no -t $task.cpus --split-prefix tmpfile ${params.host_db} \\
            ${sample_id}_IR.fastq.gz | samtools sort -@ $task.cpus -o ${sample_id}_host_removed_sr.bam --write-index -
        samtools fastq -f 4 ${sample_id}_host_removed_sr.bam | pigz > ${sample_id}_host_removed_sr_before_repair.fastq.gz
        samtools fastq -F 4 ${sample_id}_host_removed_sr.bam | pigz > ${sample_id}_host_sr.fastq.gz
      
        repair.sh in=${sample_id}_host_removed_sr_before_repair.fastq.gz out=${sample_id}_host_removed_sr.fastq.gz
      
        reformat.sh in=${sample_id}_host_removed_sr.fastq.gz out1=${sample_id}_host_removed_sr_R1.fastq out2=${sample_id}_host_removed_sr_R2.fastq
        pigz ${sample_id}_host_removed_sr_R1.fastq
        pigz ${sample_id}_host_removed_sr_R2.fastq
    fi
    """
}


process Post_host_remove_fastqc {
    tag { sample_id }
    errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.run_id}/fastqc/targeted_read_mapping", mode: 'copy'
    label 'qc'

    input:
    // qc_files will be a list (grouped per sample)
    tuple val(sample_id), path(qc_files)

    output:
    tuple val(sample_id), file("*"), emit: post_host_remove_fastqc_ch

    when:
    params.run_qc_stats

    script:
    """
    fastqc --memory 2000 --outdir . ${qc_files}
    """
}

process Multiqc_QC_host_removal {
   
    errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.run_id}/", mode: 'copy'
    label 'qc'

    input:
    file(post_map2refseq_files)
    when:
    params.run_qc_stats
    script:
    """
    mkdir -p ${params.outdir}/${params.run_id}/multiqc/targeted_read_mapping/
    mkdir -p ${params.outdir}/${params.run_id}/qc_stats/targeted_read_mapping/
    multiqc ${params.outdir}/${params.run_id}/fastqc/targeted_read_mapping/ --data-format csv --outdir ${params.outdir}/${params.run_id}/multiqc/targeted_read_mapping//
    Rscript ${params.scripts}/create_host_removed_qc_stats.R -p ${params.outdir}/${params.run_id}/multiqc/targeted_read_mapping/multiqc_data/multiqc_general_stats.csv \
                                                -o ${params.outdir}/${params.run_id}/qc_stats/targeted_read_mapping/
    """
}

process Map_Reads_2_Contigs {
    tag { "${sample_id}_${assembler}" }
    publishDir { "${params.outdir}/${params.run_id}/${sample_id}/assembly_verification/map2assembly/" },
                mode: 'copy',
                overwrite: true,
                saveAs: { fn -> "${assembler}/${fn}" }
    label 'lowmem'
    errorStrategy 'ignore'

    input:
    tuple val(sample_id),
          file(fastq_1),
          file(fastq_2),
          file(long_read),
          file(contigs),
          val(mode),
          val(assembler)

    output:
    tuple val(sample_id), val(assembler), path("*_assembly_mapped_*.fastq.gz"),
        optional: true, emit: map2assembly_ch

    tuple val(sample_id),
          path("${sample_id}_${assembler}_unmapped_LR.fastq.gz"),
          val(assembler),
          val(mode),
          emit: unmap2assembly_ch
    when:
    params.run_assembly && mode in ['long', 'hybrid']

    script:
    """

    # Empty SR placeholders preserve the six-field Assembly_Workflow input.
    # The emitted mode is forced to 'long', so these files are never assembled.
    printf '' | gzip -c > ${sample_id}_assembly_unmapped_SR_R1.fastq.gz
    printf '' | gzip -c > ${sample_id}_assembly_unmapped_SR_R2.fastq.gz
    printf '' | gzip -c > ${sample_id}_assembly_unmapped_LR.fastq.gz

    if [[ "${mode}" == "long" || "${mode}" == "hybrid" ]]; then
        minimap2 -ax map-ont -t ${task.cpus} ${contigs} ${long_read} \
            | samtools view -b -o ${sample_id}_assembly_mapping_LR.bam

        # -F 4 = mapped; -f 4 = unmapped
        samtools fastq -F 4 ${sample_id}_assembly_mapping_LR.bam \
            | pigz > ${sample_id}_assembly_mapped_LR.fastq.gz
        samtools fastq -f 4 ${sample_id}_assembly_mapping_LR.bam \
            > ${sample_id}_assembly_unmapped_LR.unfiltered.fastq

        # Strictly longer than 700 bp means a minimum accepted length of 701.
        reformat.sh \
            in=${sample_id}_assembly_unmapped_LR.unfiltered.fastq \
            out=${sample_id}_${assembler}_unmapped_LR.fastq.gz \
            minlength=700 ow=t

        rm ${sample_id}_assembly_unmapped_LR.unfiltered.fastq

        samtools sort -o ${sample_id}_assembly_mapping_LR_sorted.bam \
            ${sample_id}_assembly_mapping_LR.bam
        samtools index ${sample_id}_assembly_mapping_LR_sorted.bam
        samtools idxstats ${sample_id}_assembly_mapping_LR_sorted.bam \
            > ${sample_id}_assembly_mapping_LR_covstats.txt

    fi

    """
}
