#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// QC MODULE: Mode-aware processes

process Interleave {
    tag { sample_id }
    publishDir "${params.outdir}/${params.run_id}/interleave_fastq/", mode: 'copy'
    label 'optimized_qc_workflow'
    cpus { cpus }
    memory { mem }

    input:
    tuple val(sample_id), file(fastq_1), file(fastq_2), val(long_read), val(mode), val(cpus), val(mem)

    output:
    tuple val(sample_id), file("${sample_id}_IR.fastq.gz"), optional: true, emit: interleave_fastq_file_ch
    tuple val(sample_id), file("*"), optional: true, emit: interleave_ch
 
    when:
    params.run_qc_stats
    script:
    """
    if [ "${params.longreads}" == "true" ]; then
	echo "Skip this step" > lr_skipped.txt
    fi
    if [ "${params.hybrid}" == "true" ]; then
    	#mkdir -p ${params.outdir}/${params.run_id}/interleave_fastq/
    	reformat.sh in1=${fastq_1} in2=${fastq_2} out=${sample_id}_IR.fastq.gz ow=t tossbrokenreads=t
    elif [ "${params.shortreads}" == "true" ]; then
        #mkdir -p ${params.outdir}/${params.run_id}/interleave_fastq/
        reformat.sh in1=${fastq_1} in2=${fastq_2} out=${sample_id}_IR.fastq.gz ow=t tossbrokenreads=t
    else 
	echo "skipped" > interleave.skipped
    fi
    """
}

process Pretrim_fastqc_merged {
    tag {sample_id}
    errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.run_id}/fastqc/pretrim/", pattern: "*.{zip,html}", mode: 'copy'
    publishDir { "${params.outdir}/${params.run_id}/${sample_id}/q_stats/" }, pattern: "*.txt", mode: 'copy'
    label 'optimized_qc_workflow'

    cpus {cpus} // setting slurm allocation dynamically
    memory {mem} // setting slurm allocation dynamically

    input:
    tuple val(sample_id), val(fastq_1),val(fastq_2), file(long_read), val(mode), val(cpus), val(mem)
    tuple val(interleave_sample_id), file(interleave_files)

    output:
    val (sample_id), emit:pretrim_fastqc_sample_id
    tuple val(sample_id), file("*"), emit: pretrim_fastqc_ch

    when: 
    params.run_qc_stats
    script:
    """
    if [ "${params.longreads}" == "true" ] || [ "${params.hybrid}" == "true" ]; then
    	#mkdir -p ${params.outdir}/${params.run_id}/fastqc/pretrim/
    	fastqc --memory 2000 --outdir . \
               ${long_read} 
        seqkit stats -j ${task.cpus} -a -T ${long_read} > ${sample_id}_seqkit_rawReads_LR.txt
    fi
    if [ "${params.shortreads}" == "true" ] || [ "${params.hybrid}" == "true" ]; then
    	#mkdir -p ${params.outdir}/${params.run_id}/fastqc/pretrim/
    	fastqc --memory 2000 --outdir . ${interleave_files}
        seqkit stats -j ${task.cpus} -a -T ${interleave_files} > ${sample_id}_seqkit_rawReads_SR.txt
    fi
    """
}

process Pretrim_fastqc {
    tag {sample_id}
    errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.run_id}/fastqc/pretrim/", mode: 'copy'
    label 'optimized_qc_workflow'

    input: 
    tuple val(sample_id), path(fastq_1),path(fastq_2),path(data_dir)

    output:
    val(sample_id), emit: pretrim_fastqc_sample_id    
    tuple val(sample_id), file("*"), emit: pretrim_fastqc_ch

    when:
    params.run_qc_stats

    script:
    """
    fastqc --memory 2000 --outdir . \
           ${data_dir}/${sample_id}*.fastq.gz
    """
}

process Pretrim_NanoPlot {
    tag { sample_id }
    errorStrategy 'ignore'

    publishDir { "${params.outdir}/${params.run_id}/nanoplot/pretrim/${sample_id}/" }, mode: 'copy'

    label 'optimized_qc_workflow'
    cpus { cpus }
    memory { mem }

    input:
    tuple val(sample_id), file(fastq_1), file(fastq_2), file(long_read), val(mode), val(cpus), val(mem)

    output:
    tuple val(sample_id), path("${sample_id}_pretrim_nanoplot"), emit: pretrim_nanoplot_ch

    when:
    params.run_qc_stats && (params.longreads || params.hybrid)

    script:
    """
    mkdir -p ${sample_id}_pretrim_nanoplot

    NanoPlot \\
        --fastq ${long_read} \\
        --threads ${task.cpus} \\
        --outdir ${sample_id}_pretrim_nanoplot \\
        --prefix ${sample_id}_pretrim_
    """
}

process Quality_Control {
    tag { sample_id }
    publishDir { "${params.outdir}/${params.run_id}/${sample_id}/trim/quality_control/" }, mode: 'copy'
    cpus { cpus }
    memory { mem }

    input:
    tuple val(sample_id), file(fastq_1), file(fastq_2), file(long_read), val(mode), val(cpus), val(mem)


    output:
    tuple val(sample_id), file("${sample_id}_fastp_R1.fastq.gz"), file("${sample_id}_fastp_R2.fastq.gz"), optional: true, emit: trimmed_short_ch
    tuple val(sample_id), file("${sample_id}_fastp_merged.fastq.gz"), optional: true, emit: interleaved_trimmed_short_ch

    when:
    params.shortreads || params.hybrid

    script:
    """
    if [ "${params.trim}" == "true" ]; then

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
	out1=${sample_id}_fastp_R1.fastq.gz \
        out2=${sample_id}_fastp_R2.fastq.gz

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
    publishDir { "${params.outdir}/${params.run_id}/${sample_id}/trim/quality_control/" }, mode: 'copy'
    cpus { cpus }
    memory { mem }

    input:
    tuple val(sample_id), file(fastq_1), file(fastq_2), file(long_read), val(mode), val(cpus), val(mem)

    output: 
    tuple val(sample_id), file("${sample_id}.LR.porechop.fastq.gz"), emit: porechop_ch
    
    when:
    params.longreads || params.hybrid
    script:
    """
    porechop -i ${long_read} \
             -o ${sample_id}.LR.porechop.fastq.gz  \
             --format fastq.gz \
             -t ${task.cpus} \
             --require_two_barcodes \
             --discard_unassigned --discard_middle 
    """

}

process Fastp_LongReads {
    tag { sample_id }
    errorStrategy 'ignore'

    publishDir { "${params.outdir}/${params.run_id}/${sample_id}/trim/quality_control/" }, mode: 'copy'

    label 'optimized_qc_workflow'
    conda "${baseDir}/env/aio_qc.yml"

    input:
    tuple val(sample_id), file(porechop_fastq)

    output:
    tuple val(sample_id), file("${sample_id}.LR.trimmed.fastq.gz"), emit: fastp_long_ch
    tuple val(sample_id), file("*"), emit: fastp_long_all_files_ch

    when:
    params.longreads || params.hybrid

    script:
    """
    if [ "${params.trim}" == "true" ]; then

        fastplong \\
            -i ${porechop_fastq} \\
            -o ${sample_id}.LR.trimmed.fastq.gz \\
            --thread ${task.cpus} \\
            --disable_adapter_trimming \\
            -y \\
            -x \\
            --qualified_quality_phred ${params.quality_phread_fastp_long} \\
            --length_required 30 \\
            --json ${sample_id}_fastp_LONG.json \\
            --html ${sample_id}_fastp_LONG.html \\
            --report_title "${sample_id}_fastp_long_report" \\
            ${params.fastp_long_opts ?: ''} 

    else

        cp ${porechop_fastq} ${sample_id}.LR.fastp_long.fastq.gz

    fi
    """
}

process Post_trim_fastqc {
    tag { sample_id }
    errorStrategy 'ignore'

    publishDir "${params.outdir}/${params.run_id}/fastqc/post_trim/", pattern: "*.{zip,html}", mode: 'copy'
    publishDir { "${params.outdir}/${params.run_id}/${sample_id}/q_stats/" }, pattern: "*.txt", mode: 'copy'

    label "qc"
    conda "$baseDir/env/aio_qc.yml"

    input:
    tuple val(sample_id), path(qc_files)

    output:
    tuple val(sample_id), file("*"), emit: posttrim_fastqc_ch

    when:
    params.run_qc_stats

    script:
    """
    #mkdir -p ${params.outdir}/${params.run_id}/fastqc/post_trim/

    fastqc --memory 2000 --outdir . ${qc_files}

    seqkit stats -j ${task.cpus} -a -T ${qc_files} > ${sample_id}_seqkit_trimmedReads.txt
    """
}


process Posttrim_NanoPlot {
    tag { sample_id }
    errorStrategy 'ignore'

    publishDir { "${params.outdir}/${params.run_id}/nanoplot/post_trim/${sample_id}/" }, mode: 'copy'

    label 'optimized_qc_workflow'
    conda "${baseDir}/env/aio_qc.yml"

    input:
    tuple val(sample_id), file(trimmed_long_read)

    output:
    tuple val(sample_id), path("${sample_id}_posttrim_nanoplot"), emit: posttrim_nanoplot_ch

    when:
    params.run_qc_stats && (params.longreads || params.hybrid)

    script:
    """
    mkdir -p ${sample_id}_posttrim_nanoplot

    NanoPlot \\
        --fastq ${trimmed_long_read} \\
        --threads ${task.cpus} \\
        --outdir ${sample_id}_posttrim_nanoplot \\
        --prefix ${sample_id}_posttrim_
    """
}

// Multiqc - QC Stats

process Multiqc_QC_Stats {
   
    errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.run_id}/", mode: 'copy'
    label 'qc'
    conda "$baseDir/env/aio_qc.yml"

    //cpus {cpus} // setting slurm allocation dynamically
    //memory {mem} // setting slurm allocation dynamically

    input:
    file fastqc_pre_trim_files
    file fastqc_post_trim_files 
    
    output:
    val params.run_id, emit: multiqc_complete_ch
    //path("multiqc/pretrim/"), emit: multiqc_pretrim_ch
    //path("multiqc/post_trim/"), emit: multiqc_posttrim_ch
    //path("qc_stats/"), emit: qc_stats_ch
    when:
    params.run_qc_stats
    script:
    """
    rm -rf "${params.outdir}/${params.run_id}//multiqc/pretrim" "${params.outdir}/${params.run_id}//multiqc/post_trim" "${params.outdir}/${params.run_id}/qc_stats/"
    mkdir -p ${params.outdir}/${params.run_id}/multiqc/pretrim/
    multiqc ${params.outdir}/${params.run_id}/fastqc/pretrim/ --data-format csv --outdir ${params.outdir}/${params.run_id}/multiqc/pretrim/ --export
    mkdir -p ${params.outdir}/${params.run_id}/ multiqc/post_trim/
    multiqc ${params.outdir}/${params.run_id}/fastqc/post_trim/ --data-format csv --outdir ${params.outdir}/${params.run_id}/multiqc/post_trim/ --export
    mkdir -p ${params.outdir}/${params.run_id}/qc_stats/
    Rscript ${params.scripts}/create_qc_stats.R -i ${params.outdir}/${params.run_id}/multiqc/pretrim/multiqc_data/multiqc_general_stats.csv \
                                                                                -p ${params.outdir}/${params.run_id}/multiqc/post_trim/multiqc_data/multiqc_general_stats.csv \
                                                                                -o ${params.outdir}/${params.run_id}/qc_stats/
    echo "Post Trim Multiqc Complete" > post_trim_multiqc.finished
    Rscript ${params.scripts}/Plot_QC_Counts.R \
        --input ${params.outdir}/${params.run_id}/qc_stats/qc_stats_final.xlsx \
        --outdir ${params.outdir}/${params.run_id}/qc_plots/ \
        --outfile ${params.run_id}_raw_reads_vs_trimmed_reads.jpeg
    """
}

process Read_Distribution {
    publishDir { "${params.outdir}/${params.run_id}/${sample_id}/read_distribution" }, mode: 'copy'
    label 'optimized_qc_workflow'
    conda "$baseDir/env/aio_qc.yml"
    errorStrategy 'ignore'
    tag { sample_id }   
    input:
    tuple val(sample_id), val(qc_long_read)
    
    output: 
    tuple val(sample_id), 
	  file("${sample_id}_read_length_distribution.png"), 
	  file("${sample_id}_read_length_summary.csv"), 
          file("${sample_id}_Q1.csv"), 
          file("${sample_id}_trimmed_q1_removed.fastq.gz"),
          //path("fastcat"), 
          emit: read_distribution_ch
    
    tuple val(sample_id), path("${sample_id}_trimmed_q1_removed.fastq.gz"), emit: read_distribution_fastq_ch

    when:
    params.longreads || params.hybrid

    script:
    """
    mkdir -p tmp
    export TMPDIR="\$PWD/tmp"
    export TMP="\$TMPDIR"
    export TEMP="\$TMPDIR"
    Rscript ${params.scripts}/read_distribution.R \\
	-i ${qc_long_read} \\
	-p ${sample_id} \\
	-o .

    q1=\$(tail -n +2 "${sample_id}_Q1.csv" | head -n 1 | tr -d '\\r')
    
    echo "Using Q1 minimum read length for ${sample_id}: \${q1}"
  
    seqkit seq \\
    	-m "\${q1}" \\
        "${qc_long_read}" \\
        -o "${sample_id}_trimmed_q1_removed.fastq.gz"
    #fastcat fastq \\
    #    --force-error \\
    #    --min-length="\${q1}" \\
    #    --min-qscore=${params.dragonflye_min_quality} \\
    #    --output=./fastcat -v \\
    #    "${qc_long_read}" | pigz >  ${sample_id}_trimmed_q1_removed.fastq.gz
    """
}

process Fastcat_trim {
    publishDir { "${params.outdir}/${params.run_id}/${sample_id}/read_distribution" }, mode: 'copy'
    label 'optimized_qc_workflow'
    conda "$baseDir/env/fastcat.yml"
    errorStrategy 'ignore'
    tag { sample_id }
    input:
    tuple val(sample_id), val(min_length), val(min_quality), val(trim_file) 

    output:
    tuple val(sample_id),
          path("fastcat"),
          emit: read_distribution_ch

    when:
    params.longreads || params.hybrid

    script:
    """
    fastcat fastq \\
        --force-error \\
        --min-length=${min_length} \\
        --min-qscore=${min_quality} \\
        --output=./fastcat -v \\
        ${trim_file} | pigz >  ${sample_id}.LR.trimmedFiltered.fastq.gz
    """
}




process LongRead_NanoPlot_QC_Stats {

    errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.run_id}/", mode: 'copy'
    label 'qc'
    conda "$baseDir/env/aio_qc.yml"

    input:
    val nanoplot_pre_files
    val nanoplot_post_files

    output:
    val params.run_id, emit: longread_nanoplot_qc_complete_ch

    when:
    params.run_qc_stats && (params.longreads || params.hybrid)

    script:
    """
    set -euo pipefail

    PRE_NANOPLOT_DIR="${params.outdir}/${params.run_id}/nanoplot/pretrim/"
    POST_NANOPLOT_DIR="${params.outdir}/${params.run_id}/nanoplot/post_trim/"

    mkdir -p ${params.outdir}/${params.run_id}/qc_stats_post_Q1/
    mkdir -p ${params.outdir}/${params.run_id}/qc_plots/

    if [[ ! -d "\${PRE_NANOPLOT_DIR}" ]]; then
        echo "ERROR: Missing pretrim NanoPlot directory: \${PRE_NANOPLOT_DIR}"
        exit 1
    fi

    if [[ ! -d "\${POST_NANOPLOT_DIR}" ]]; then
        echo "ERROR: Missing posttrim NanoPlot directory: \${POST_NANOPLOT_DIR}"
        exit 1
    fi

    if ! find "\${PRE_NANOPLOT_DIR}" -name "*NanoStats.txt" -type f | grep -q . ; then
        echo "ERROR: No pretrim NanoStats.txt files found in \${PRE_NANOPLOT_DIR}"
        exit 1
    fi

    if ! find "\${POST_NANOPLOT_DIR}" -name "*NanoStats.txt" -type f | grep -q . ; then
        echo "ERROR: No posttrim NanoStats.txt files found in \${POST_NANOPLOT_DIR}"
        exit 1
    fi

    Rscript ${params.scripts}/create_nanoplot_qc_stats.R \
        -i \${PRE_NANOPLOT_DIR} \
        -p \${POST_NANOPLOT_DIR} \
        -o ${params.outdir}/${params.run_id}/qc_stats_post_Q1/

    Rscript ${params.scripts}/Plot_NanoPlot_QC_Counts.R \
        --input ${params.outdir}/${params.run_id}/qc_stats_post_Q1/qc_stats_final.xlsx \
        --outdir ${params.outdir}/${params.run_id}/qc_plots/ \
        --outfile ${params.run_id}_long_raw_reads_vs_Q1_trimmed_reads.jpeg

    echo "Long-read NanoPlot QC stats complete" > longread_nanoplot_qc_stats.finished
    """
}
