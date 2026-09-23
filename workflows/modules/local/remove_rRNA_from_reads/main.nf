#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process Remove_Common_Flora_rRNA_reads {
    tag { sample_id }
    publishDir { "${params.outdir}/${params.project_id}/${sample_id}/trim/quality_control/rRNA_removed/" }, mode: 'copy'
    label 'qc'
    cpus { cpus }
    memory { mem }

    input:
    tuple val(sample_id), val(sr_file), val(lr_file), val(cpus), val(mem)

    output:
    tuple val(sample_id), file("${sample_id}_host_contaminant_rRNA_removed_pe.fastq.gz"), optional: true, emit: rRNA_host_remove_short
    tuple val(sample_id), file("${sample_id}_host_contaminant_rRNA_removed_LR.fastq.gz"), optional: true, emit: rRNA_host_remove_long
    tuple val(sample_id), file("${sample_id}_rRNA.fastq.gz"), optional: true, emit: rRNA_reads_short
    tuple val(sample_id), file("${sample_id}_rRNA_LR.fastq.gz"), optional: true, emit: rRNA_reads_long
    tuple val(sample_id), file("${sample_id}_host_contaminant_rRNA_removed_R1.fastq.gz"), optional: true, emit: rRNA_host_remove_short_r1
    tuple val(sample_id), file("${sample_id}_host_contaminant_rRNA_removed_R2.fastq.gz"), optional: true, emit: rRNA_host_remove_short_r2
    tuple val(sample_id), file("*"), optional: true, emit: rRNA_all_files_ch 
    
    when:
    params.remove_rRNA_reads

    script:
    def java_mem = mem.tokenize()[0] + 'g'
    """
    echo "sr_file: ${sr_file}"
    echo "lr_file: ${lr_file}"

    if [[ -s "$lr_file" ]]; then
       minimap2 -ax map-ont --secondary=no -t ${task.cpus} ${params.host_db_silva} "${lr_file}" > ${sample_id}_host_contaminant_rRNA_removed_LR.sam
       samtools fastq -f 4 ${sample_id}_host_contaminant_rRNA_removed_LR.sam | pigz > ${sample_id}_host_contaminant_rRNA_removed_LR.fastq.gz
       samtools fastq -F 4 ${sample_id}_host_contaminant_rRNA_removed_LR.sam | pigz > ${sample_id}_rRNA_LR.fastq.gz
    fi


    if [[ -s "$sr_file" ]]; then
	bbmap.sh ${params.bbmap_args} \\
                    in="${sr_file}" \\
                    path=${params.bbmap_ref_silva} tossbrokenreads printunmappedcount=t \\
                    covstats=${sample_id}.rRNA_covstats.txt \\
                    outm=${sample_id}_rRNA.fastq.gz usejni=f \\
                    -Xmx${java_mem} \\
                    outu=${sample_id}_host_contaminant_rRNA_removed_pe_before_repair.fastq.gz overwrite=true


        repair.sh in=${sample_id}_host_contaminant_rRNA_removed_pe_before_repair.fastq.gz out=${sample_id}_host_contaminant_rRNA_removed_pe.fastq.gz

        reformat.sh in=${sample_id}_host_contaminant_rRNA_removed_pe.fastq.gz \
        out1=${sample_id}_host_contaminant_rRNA_removed_R1.fastq \
        out2=${sample_id}_host_contaminant_rRNA_removed_R2.fastq

        gzip ${sample_id}_host_contaminant_rRNA_removed_R1.fastq
        gzip ${sample_id}_host_contaminant_rRNA_removed_R2.fastq

    fi
    """
}

process Post_removal_rRNA_reads_fastqc {
    tag { sample_id }
    errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.project_id}/fastqc/rRNA_removal", mode: 'copy'
    label 'qc'

    input:
    // qc_files will be a list (grouped per sample)
    tuple val(sample_id), path(qc_files)

    output:
    tuple val(sample_id), file("*"), emit: post_rRNA_removal_fastqc_ch

    when:
    params.run_qc_stats

    script:
    """
    fastqc --memory 2000 --outdir . ${qc_files}
    """
}

process Multiqc_QC_rRNA_removal {
   
    errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.project_id}/", mode: 'copy'
    label 'qc'
    conda "$baseDir/env/aio_qc.yml"

    //cpus {cpus} // setting slurm allocation dynamically
    //memory {mem} // setting slurm allocation dynamically

    input:
    file(post_map2refseq_files)
    when:
    params.run_qc_stats
    script:
    """
    mkdir -p ${params.outdir}/${params.project_id}/multiqc/rRNA_removal/
    mkdir -p ${params.outdir}/${params.project_id}/qc_stats/targeted_read_mapping/
    multiqc ${params.outdir}/${params.project_id}/fastqc/rRNA_removal/ --data-format csv --outdir ${params.outdir}/${params.project_id}/multiqc/rRNA_removal/
    Rscript ${params.scripts}/create_remove_rRNA_qc_stats.R -r ${params.outdir}/${params.project_id}/multiqc/rRNA_removal/multiqc_data/multiqc_general_stats.csv \
                                                -o ${params.outdir}/${params.project_id}/qc_stats/rRNA_removal/
    """
}

process RiboDetector_Remove_rRNA_reads {
    tag { sample_id }

    publishDir { "${params.outdir}/${params.project_id}/${sample_id}/trim/quality_control/rRNA_removed/" }, mode: 'copy'

    label 'ribodetector'

    input:
    tuple val(sample_id), val(sr_file), val(lr_file), val(cpus), val(mem)

    //cpus { cpus }
    //memory { mem }

    output:
    tuple val(sample_id), file("${sample_id}_host_contaminant_rRNA_removed_pe.fastq.gz"), optional: true, emit: rRNA_host_remove_short
    tuple val(sample_id), file("${sample_id}_host_contaminant_rRNA_removed_LR.fastq.gz"), optional: true, emit: rRNA_host_remove_long

    tuple val(sample_id), file("${sample_id}_rRNA.fastq.gz"), optional: true, emit: rRNA_reads_short
    tuple val(sample_id), file("${sample_id}_rRNA_LR.fastq.gz"), optional: true, emit: rRNA_reads_long

    tuple val(sample_id), file("${sample_id}_host_contaminant_rRNA_removed_R1.fastq.gz"), optional: true, emit: rRNA_host_remove_short_r1
    tuple val(sample_id), file("${sample_id}_host_contaminant_rRNA_removed_R2.fastq.gz"), optional: true, emit: rRNA_host_remove_short_r2

    tuple val(sample_id), file("*"), optional: true, emit: rRNA_all_files_ch

    when:
    params.remove_rRNA_reads

    script:
    def mode = (params.ribodetector_mode ?: 'cpu').toString().toLowerCase()

    if (!(mode in ['cpu', 'gpu'])) {
        error "Invalid --ribodetector_mode '${params.ribodetector_mode}'. Use: cpu or gpu"
    }

    def ribodetector_cmd = mode == 'gpu' ? 'ribodetector' : 'ribodetector_cpu'

    def ribodetector_short_len = params.ribodetector_short_len ?: 100
    def ribodetector_long_len  = params.ribodetector_long_len  ?: 100
    def ribodetector_ensure    = params.ribodetector_ensure    ?: 'rrna'
    def ribodetector_chunk = params.ribodetector_chunk_size ?: (mode == 'gpu' ? 32 : 256)
    def ribodetector_memory_opt = mode == 'gpu'
        ? "-m ${params.ribodetector_gpu_memory ?: 3}"
        : ""

    """
    set -euo pipefail

    echo "sample_id: ${sample_id}"
    echo "sr_file: ${sr_file}"
    echo "lr_file: ${lr_file}"
    echo "ribodetector_mode: ${mode}"
    echo "ribodetector_cmd: ${ribodetector_cmd}"
    echo "ribodetector_chunk_size: ${ribodetector_chunk}"
    echo "ribodetector_memory_opt: ${ribodetector_memory_opt}"

    if [[ "${sr_file}" != "null" && -s "${sr_file}" ]]; then

        echo "Running Ribodetector on short/interleaved reads"

        ${ribodetector_cmd} \\
            -t ${task.cpus} \\
            -l ${ribodetector_short_len} \\
            ${ribodetector_memory_opt} \\
            -i "${sr_file}" \\
            -e ${ribodetector_ensure} \\
            --chunk_size ${ribodetector_chunk} \\
            --log ${sample_id}.ribodetector.short.log \\
            -o ${sample_id}_host_contaminant_rRNA_removed_pe.raw.fastq \\
            -r ${sample_id}_rRNA.raw.fastq

        repair.sh \\
            in=${sample_id}_host_contaminant_rRNA_removed_pe.raw.fastq \\
            out=${sample_id}_host_contaminant_rRNA_removed_pe.fastq \\
            overwrite=true

        reformat.sh \\
            in=${sample_id}_host_contaminant_rRNA_removed_pe.fastq \\
            out1=${sample_id}_host_contaminant_rRNA_removed_R1.fastq \\
            out2=${sample_id}_host_contaminant_rRNA_removed_R2.fastq \\
            overwrite=true

        pigz -p ${task.cpus} ${sample_id}_host_contaminant_rRNA_removed_pe.fastq
        pigz -p ${task.cpus} ${sample_id}_host_contaminant_rRNA_removed_R1.fastq
        pigz -p ${task.cpus} ${sample_id}_host_contaminant_rRNA_removed_R2.fastq

        if [[ -s ${sample_id}_rRNA.raw.fastq ]]; then
            mv ${sample_id}_rRNA.raw.fastq ${sample_id}_rRNA.fastq
            pigz -p ${task.cpus} ${sample_id}_rRNA.fastq
        fi
    else
        echo "No short/interleaved reads supplied for ${sample_id}; skipping short-read Ribodetector"
    fi


    if [[ "${lr_file}" != "null" && -s "${lr_file}" ]]; then

        echo "Running Ribodetector on long reads"

        ${ribodetector_cmd} \\
            -t ${task.cpus} \\
            -l ${ribodetector_long_len} \\
            ${ribodetector_memory_opt} \\
            -i "${lr_file}" \\
            -e ${ribodetector_ensure} \\
            --chunk_size ${ribodetector_chunk} \\
            --log ${sample_id}.ribodetector.long.log \\
            -o ${sample_id}_host_contaminant_rRNA_removed_LR.fastq \\
            -r ${sample_id}_rRNA_LR.fastq

        pigz -p ${task.cpus} ${sample_id}_host_contaminant_rRNA_removed_LR.fastq

        if [[ -s ${sample_id}_rRNA_LR.fastq ]]; then
            pigz -p ${task.cpus} ${sample_id}_rRNA_LR.fastq
        fi
    else
        echo "No long reads supplied for ${sample_id}; skipping long-read Ribodetector"
    fi
    """
}
