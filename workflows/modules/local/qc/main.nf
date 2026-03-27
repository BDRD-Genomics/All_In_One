#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// QC MODULE: Mode-aware processes

process Interleave {
    tag { sample_id }
    publishDir "${params.outdir}/${params.project_id}/interleave_fastq/", mode: 'copy'
    label 'optimized_qc_workflow'
    conda "${baseDir}/env/aio_qc.yml"

    input:
    tuple val(sample_id), file(fastq_1), file(fastq_2), val(long_read), val(mode), val(cpus), val(mem)

    cpus {cpus} // setting slurm allocation dynamically 
    memory {mem} // setting slurm allocation dynamically 

    output:
    tuple val(sample_id), file("${sample_id}_IR.fastq.gz"), optional: true, emit: interleave_fastq_file_ch
    tuple val(sample_id), file("*"), optional: true, emit: interleave_ch
 
    when:
    params.run_qc_stats
    script:
    """
    if [ "${params.longreads}" == "true" ]; then
	echo "Skip this step" > step_skipped.txt
    
    elif [ "${params.hybrid}" == "true" ]; then
    	mkdir -p ${params.outdir}/${params.project_id}/interleave_fastq/
    	reformat.sh in1=${fastq_1} in2=${fastq_2} out=${sample_id}_IR.fastq.gz ow=t tossbrokenreads=t
    elif [ "${params.shortreads}" == "true" ]; then
        mkdir -p ${params.outdir}/${params.project_id}/interleave_fastq/
        reformat.sh in1=${fastq_1} in2=${fastq_2} out=${sample_id}_IR.fastq.gz ow=t tossbrokenreads=t
    else 
	echo "skipped" > interleave.skipped
    fi
    """
}

process Pretrim_fastqc_merged {
    tag {sample_id}
    errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.project_id}/fastqc/pretrim/", mode: 'copy'
    label 'optimized_qc_workflow'
    conda "${baseDir}/env/aio_qc.yml"

    cpus {cpus} // setting slurm allocation dynamically
    memory {mem} // setting slurm allocation dynamically

    input:
    tuple val(sample_id), val(fastq_1),val(fastq_2), file(long_read), val(mode), val(cpus), val(mem)
    tuple val(sample_id), file(interleave_files)

    output:

    tuple val(sample_id), file("*"), emit: pretrim_fastqc_ch

    when: 
    params.run_qc_stats
    script:
    """
    if [ "${params.longreads}" == "true" ] || [ "${params.hybrid}" == "true" ]; then
    	mkdir -p ${params.outdir}/${params.project_id}/fastqc/pretrim/
    	fastqc --memory 2000 --outdir . \
               ${long_read} 
    fi
    if [ "${params.shortreads}" == "true" ] || [ "${params.hybrid}" == "true" ]; then
    	mkdir -p ${params.outdir}/${params.project_id}/fastqc/pretrim/
    	fastqc --memory 2000 --outdir . ${interleave_files}
    fi
    """
}

process Pretrim_fastqc {
    tag {sample_id}
    errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.project_id}/fastqc/pretrim/", mode: 'copy'
    label 'optimized_qc_workflow'
    conda "${baseDir}/env/aio_qc.yml"

    cpus {cpus} // setting slurm allocation dynamically
    memory {mem} // setting slurm allocation dynamically

    input: 
    tuple val(sample_id), path(fastq_1),path(fastq_2),path(data_dir)

    output:
    
    tuple val(sample_id), file("*"), emit: pretrim_fastqc_ch

    when:
    params.run_qc_stats

    script:
    """
    mkdir -p ${params.outdir}/${params.project_id}/fastqc/pretrim/
    fastqc --memory 2000 --outdir . \
           ${data_dir}/${sample_id}*.fastq.gz
    """
}

process Quality_Control {
    tag { sample_id }
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/trim/quality_control/", mode: 'copy'
    //label 'optimized_qc_workflow'
    conda "${baseDir}/env/aio_qc.yml"
    input:
    tuple val(sample_id), file(fastq_1), file(fastq_2), file(long_read), val(mode), val(cpus), val(mem)

    cpus {cpus} // setting slurm allocation dynamically 
    memory {mem} // setting slurm allocation dynamically 


    output:
    tuple val(sample_id), file("${sample_id}_fastp_R1.fastq.gz"), file("${sample_id}_fastp_R2.fastq.gz"), optional: true, emit: trimmed_short_ch
    tuple val(sample_id), file("${sample_id}_fastp_merged.fastq.gz"), optional: true, emit: interleaved_trimmed_short_ch

    when:
    params.shortreads || params.hybrid

    script:
    """
    if [ "${params.trim}" == "true" ]; then
    mkdir -p ${params.outdir}/${params.project_id}/${sample_id}/status_log/

    fastp --verbose -xyp --dedup --thread=${task.cpus} --report_title "${sample_id}_fastp_report" \
		${params.fastp_opts} \
		-j ${sample_id}_fastp_ILLUMINA.json \
		-h ${sample_id}_fastp_ILLUMINA.html \
		--in1 ${fastq_1} --in2 ${fastq_2} \
		--out1 ${sample_id}_fastp_before_repair_R1.fastq.gz \
		--out2 ${sample_id}_fastp_before_repair_R2.fastq.gz >> fastp.log 2>&1
    reformat.sh in1=${sample_id}_fastp_before_repair_R1.fastq.gz in2=${sample_id}_fastp_before_repair_R2.fastq.gz out=${sample_id}_fastp_merged_before_repair.fastq.gz ow=t
    repair.sh in=${sample_id}_fastp_merged_before_repair.fastq.gz out=${sample_id}_fastp_merged.fastq.gz
   
    reformat.sh in=${sample_id}_fastp_merged.fastq.gz \
	out1=${sample_id}_fastp_R1.fastq \
        out2=${sample_id}_fastp_R2.fastq

    gzip ${sample_id}_fastp_R1.fastq
    gzip ${sample_id}_fastp_R2.fastq 
    echo "Stage 2 Quality Control Finished" > ${params.outdir}/${params.project_id}/${sample_id}/status_log/stage2_qc.finished
    fi 
    if [ "${params.trim}" == "false" ]; then
    cp ${fastq_1} ${sample_id}_fastp_R1.fastq.gz
    cp ${fastq_2} ${sample_id}_fastp_R2.fastq.gz
    reformat.sh in1=${sample_id}_fastp_R1.fastq.gz in2=${sample_id}_fastp_R2.fastq.gz out=${sample_id}_fastp_merged.fastq.gz ow=t
    fi
    """
}


process PoreChop {
    tag {sample_id} 
    errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/trim/quality_control/", mode: 'copy'
    //label 'optimized_qc_workflow'
    conda "${baseDir}/env/aio_qc.yml"

    input:
    tuple val(sample_id), file(fastq_1), file(fastq_2), file(long_read), val(mode), val(cpus), val(mem)

    cpus {cpus} // setting slurm allocation dynamically 
    memory {mem} // setting slurm allocation dynamically 

    output: 
    tuple val(sample_id), file("${sample_id}.LR.trimmed.fastq.gz"), emit: porechop_ch
    
    when:
    params.longreads || params.hybrid
    script:
    """
    porechop -i ${long_read} \
             -o ${sample_id}.LR.trimmed.fastq.gz  \
             --format fastq.gz \
             -t ${task.cpus} \
             --require_two_barcodes \
             --discard_unassigned --discard_middle >> PoreChop.log 2>&1
    """

}

process Post_trim_fastqc {
    tag {sample_id} 
    errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.project_id}/fastqc/post_trim/", mode: 'copy'
    label "qc"
    conda "$baseDir/env/aio_qc.yml"
    //cpus {cpus} // setting slurm allocation dynamically
    //memory {mem} // setting slurm allocation dynamically

    input:
    tuple val(sample_id), path(qc_files)

    output:

    tuple val(sample_id), file("*"), emit: posttrim_fastqc_ch

    when:
    params.run_qc_stats

    script:
    """
    mkdir -p ${params.outdir}/${params.project_id}/fastqc/post_trim/
    fastqc --memory 2000 --outdir . \
           ${params.outdir}/${params.project_id}/${sample_id}/trim/quality_control/${sample_id}*.fastq.gz

    """
}

// Multiqc - QC Stats

process Multiqc_QC_Stats {
   
    errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.project_id}/", mode: 'copy'
    label 'qc'
    conda "$baseDir/env/aio_qc.yml"

    //cpus {cpus} // setting slurm allocation dynamically
    //memory {mem} // setting slurm allocation dynamically

    input:
    file fastqc_pre_trim_fles
    file fastqc_post_trim_files 
    
    when:
    params.run_qc_stats
    script:
    """
    mkdir -p ${params.outdir}/${params.project_id}/multiqc/pretrim/
    multiqc ${params.outdir}/${params.project_id}/fastqc/pretrim/ --data-format csv --outdir ${params.outdir}/${params.project_id}/multiqc/pretrim/
    mkdir -p ${params.outdir}/${params.project_id}/multiqc/post_trim/
    multiqc ${params.outdir}/${params.project_id}/fastqc/post_trim/ --data-format csv --outdir ${params.outdir}/${params.project_id}/multiqc/post_trim/
    mkdir -p ${params.outdir}/${params.project_id}/qc_stats/
    Rscript ${params.scripts}/create_qc_stats.R -i ${params.outdir}/${params.project_id}/multiqc/pretrim/multiqc_data/multiqc_general_stats.csv \
                                                                                -p ${params.outdir}/${params.project_id}/multiqc/post_trim/multiqc_data/multiqc_general_stats.csv \
                                                                                -o ${params.outdir}/${params.project_id}/qc_stats/
    echo "Post Trim Multiqc Complete" > post_trim_multiqc.finished
    Rscript ${params.scripts}/Plot_QC_Counts.R \
        --input ${params.outdir}/${params.project_id}/qc_stats/qc_stats_final.xlsx \
        --outdir ${params.outdir}/${params.project_id}/qc_plots/ \
        --outfile ${params.project_id}_raw_reads_vs_trimmed_reads.jpeg
    """
}





process Read_Distribution {
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/read_distribution", mode: 'copy'
    label 'optimized_qc_workflow'
    conda "$baseDir/env/aio_qc.yml"
    errorStrategy 'ignore'
    tag { sample_id }   
    input:
    tuple val(sample_id), val(qc_long_read)

    when:
    params.longreads || params.hybrid
    
    output: 
    tuple val(sample_id), file("${sample_id}_read_length_distribution.png"), file("${sample_id}_read_length_summary.csv"), file("${sample_id}_Q1.csv"), emit: read_distribution_ch
    script:
    """
    # TMP Solution (pun intended)
    export TMPDIR=/export/tmp
    export TMP=/export/tmp
    export TEMP=/export/tmp
    Rscript ${params.scripts}/read_distribution.R -i ${qc_long_read} -p ${sample_id} -o .
    """
}
