#!/usr/bin/ nextflow

nextflow.enable.dsl=2

/*
========================================================================================
   Map 2 Contigs
========================================================================================

*/

/* map reads to reference sequence and use for downstream analysis (Optional)*/

process Map_Reads_2_Contigs {
    tag { "${sample_id}_${assembler}" }
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/VS_supplemental_outputs", 
                mode: 'copy',
                overwrite: true,
                saveAs: { fn -> "${assembler}/${fn}" }
    //label 'optimized_map2reads' // can erase if new tuple works
    errorStrategy 'ignore'
    conda "$baseDir/env/md.yaml"

    input:
    tuple val(sample_id), file(fastq_1), file(fastq_2), file(long_read), file(contigs), val(mode), val(cpus), val(mem), val(assembler)

    cpus {cpus}  // setting slurm allocation dynamically
    memory {mem} // setting slurm allocation dynamically

    output:
    tuple val(sample_id), file("${sample_id}_assembly_mapping_LR_covstats.txt"), file("lr_count.txt"), optional: true, emit: map2assembly_ch_lr
    tuple val(sample_id), file("${sample_id}_assembly_mapping_SR_covstats.txt"), file("sr_count.txt"), optional: true, emit: map2assembly_ch_sr
    tuple val(sample_id), file("${sample_id}_assembly_mapped_LR.fastq"), file("${sample_id}_assembly_unmapped_LR.fastq"), optional: true, emit: mapping_reads2contigs_LR 
    tuple val(sample_id), file("${sample_id}_assembly_mapped_SR.fastq"), file("${sample_id}_assembly_unmapped_SR.fastq"), optional: true, emit: mapping_reads2contigs_SR

    when:
    params.vs

    script:
    def LR_QC_reads_com = "zgrep -c '^@' ${long_read}"
    //def LR_QC_count = LR_QC_reads_com.execute().toInteger()
    def LR_QC_count = LR_QC_reads_com
    def SR_QC_reads_com = "zgrep -c '^@' ${fastq_1}"
    //def SR_QC_count = SR_QC_reads_com.execute().toInteger() * 2
    def SR_QC_count = SR_QC_reads_com * 2
    
    """
    # TMP FIX
    #source /opt/conda/etc/profile.d/conda.sh
    #conda activate all_in_one_pipeline
    mkdir -p ${params.outdir}/${params.project_id}/${sample_id}/status_log/

    if [[ "${mode}" == "long" || "${mode}" == "hybrid" ]]; then
        LR_QC_reads_com=\$(zgrep -c '^@' ${long_read})
        echo \$LR_QC_reads_com > lr_count.txt
        minimap2 -ax map-ont -t $task.cpus ${contigs} \\
            ${long_read} > ${sample_id}_assembly_mapping_LR.sam
        
        samtools fastq -f 4 ${sample_id}_assembly_mapping_LR.sam > ${sample_id}_assembly_mapped_LR.fastq
        samtools fastq -F 4 ${sample_id}_assembly_mapping_LR.sam > ${sample_id}_assembly_unmapped_LR.fastq
        samtools sort ${sample_id}_assembly_mapping_LR.sam > ${sample_id}_assembly_mapping_LR_sorted.sam
        samtools index ${sample_id}_assembly_mapping_LR_sorted.sam
        samtools idxstats ${sample_id}_assembly_mapping_LR_sorted.sam > ${sample_id}_assembly_mapping_LR_covstats.txt

        pigz *.fastq
        echo "Mapping to Assembly: LR complete" > ${params.outdir}/${params.project_id}/${sample_id}/status_log/LR_map2assembly.finished
    fi

    if [[ "${mode}" == "short" || "${mode}" == "hybrid" ]]; then
        SR_QC_reads_com=\$(zgrep -c '^@' ${fastq_1})
        SR_QC_count=\$(( SR_QC_reads_com * 2 ))
        echo \$SR_QC_count > sr_count.txt
        reformat.sh in1=${fastq_1} in2=${fastq_2} out=${sample_id}_IR.fastq.gz ow=t
        minimap2 -ax sr -t $task.cpus ${contigs} \\
            ${sample_id}_IR.fastq.gz > ${sample_id}_assembly_mapping_SR.sam 
        
        samtools fastq -f 4 ${sample_id}_assembly_mapping_SR.sam > ${sample_id}_assembly_mapped_SR.fastq
        samtools fastq -F 4 ${sample_id}_assembly_mapping_SR.sam > ${sample_id}_assembly_unmapped_SR.fastq
        samtools sort ${sample_id}_assembly_mapping_SR.sam > ${sample_id}_assembly_mapping_SR_sorted.sam
        samtools index ${sample_id}_assembly_mapping_SR_sorted.sam
        samtools idxstats ${sample_id}_assembly_mapping_SR_sorted.sam > ${sample_id}_assembly_mapping_SR_covstats.txt

        pigz ${sample_id}_assembly_*_SR.fastq
        echo "Mapping to Assembly: SR complete" > ${params.outdir}/${params.project_id}/${sample_id}/status_log/SR_map2assembly.finished
    fi
    """
}

/*
========================================================================================
   Bash Filter Viruses
========================================================================================

*/

/* filter mmseqs for "high-confidence" viral calls and output filtered blast tables with associated fasta files and a final semi-quantitative read counts file*/

process Filter_Putative_Viruses {
    tag { "${sample_id}_${assembler}" }
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/VS_supplemental_outputs", 
                mode: 'copy',
                overwrite: true,
                saveAs: { fn -> "${assembler}/${fn}" }
    label 'optimized_map2reads'
    errorStrategy 'ignore'
    conda "$baseDir/env/vs.yml"
    input:
    tuple val(sample_id), file(mmseqs_parsed), val(assembler)

    //cpus {cpus}  // setting slurm allocation dynamically
    //memory {mem} // setting slurm allocation dynamically

    output:
    tuple val(sample_id), file("viral_blast.out"), val(assembler),  optional: true, emit: viral_blast_tab_ch

    when:
    params.vs

    script:
    """
    # TMP FIX
    #source /opt/conda/etc/profile.d/conda.sh
    #conda activate all_in_one_pipeline
    bash ${params.scripts}/commands_to_filter_blast.sh ${mmseqs_parsed}
    """
}





/*
========================================================================================
   Python generate VS output
========================================================================================

*/

/* filter mmseqs for "high-confidence" viral calls and output filtered blast tables with associated fasta files and a final semi-quantitative read counts file*/

process Generate_VS_Outputs {
    tag { "${sample_id}_${assembler}" }
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/VS_supplemental_outputs", 
                mode: 'copy',
                overwrite: true,
                saveAs: { fn -> "${assembler}/${fn}" }
    label 'optimized_map2reads'
    errorStrategy 'ignore'
    conda "$baseDir/env/vs.yml"

    input:
    tuple val(sample_id), file(viral_blast_tab), file(contigs_fasta), val(PE_covstats), val(PE_QC_reads), val(LR_covstats), val(LR_QC_reads), val(mode), val(assembler)
    //cpus {cpus}
    //memory {mem}

    output:
    tuple val(sample_id), file("*"),val(assembler), optional: true, emit: vs_reports_ch
    tuple val(sample_id), file("${sample_id}_accurate_read_counts.tsv"), val(assembler), optional: true, emit: accurate_read_counts_ch
    //tuple path("${sample_id}_accurate_read_counts.tsv"), val(assembler), optional: true, emit: accurate_read_counts_dir_ch
    tuple val(assembler), path("${sample_id}_accurate_read_counts.tsv"), optional: true, emit: accurate_read_counts_file_ch
    when:
    params.vs

    script:
    """
    # TMP FIX
    #source /opt/conda/etc/profile.d/conda.sh
    #conda activate all_in_one_pipeline
    if [[ "${mode}" == "short" ]]
    then
        python3 ${params.scripts}/filter_blast_viruses_working.py --sample-id ${sample_id} -t ${viral_blast_tab}  -c ${contigs_fasta} -p ${PE_covstats} -S ${PE_QC_reads} -n ${params.ncbi_taxa} -v ${params.virus_genome_sizes} -b ${params.bad_accessions_list}
    fi
    if [[ "${mode}" == "long" ]]
    then
        python3 ${params.scripts}/filter_blast_viruses_working.py --sample-id ${sample_id} -t ${viral_blast_tab}  -c ${contigs_fasta} -l ${LR_covstats} -L ${LR_QC_reads} -n ${params.ncbi_taxa} -v ${params.virus_genome_sizes} -b ${params.bad_accessions_list}
    fi
    if [[ "${mode}" == "hybrid" ]]
    then
        python3 ${params.scripts}/filter_blast_viruses_working.py --sample-id ${sample_id} -t ${viral_blast_tab}  -c ${contigs_fasta} -p ${PE_covstats} -S ${PE_QC_reads} -l ${LR_covstats} -L ${LR_QC_reads} -n ${params.ncbi_taxa} -v ${params.virus_genome_sizes} -b ${params.bad_accessions_list}
    fi
    """
}

process Merge_Arc {
  tag { "merge_all_${assembler}" }
  publishDir "${params.outdir}/${params.project_id}/VS_merged_arc_files/",
              mode: 'copy',
              overwrite: true,
              saveAs: { fn -> "${assembler}/${fn}" }  
  conda "$baseDir/env/vs.yml"
  errorStrategy 'ignore'
  input:
  //path accurate_read_counts_dir  // a list of all matched TSVs
  tuple val(assembler), path(arc_files)
  output:
  //path "${params.project_id}_total_hits_merged.csv"
  //path "${params.project_id}_normalized_family_hits_merged.csv"
  //path "${params.project_id}_normalized_RPM_hits_merged.csv"
  path "${params.project_id}_${assembler}_total_hits_merged.csv",              emit: total_csv
  path "${params.project_id}_${assembler}_normalized_family_hits_merged.csv",  emit: nfam_csv
  path "${params.project_id}_${assembler}_normalized_RPM_hits_merged.csv",     emit: nrpm_csv

  script:
  """
  # TMP FIX
  #source /opt/conda/etc/profile.d/conda.sh
  #conda activate all_in_one_pipeline
  set -euo pipefail
  rm -rf arc_${assembler}
  mkdir -p arc_${assembler}
  cp -f ${arc_files} arc_${assembler}
  #printf '%s\0' ${arc_files} | xargs -0 -I '{}' cp -L -n '{}' "arc_${assembler}/"
  python ${params.scripts}/combine_ARC_hits.py \
    -r arc_${assembler} \
    -p ${params.project_id}_${assembler} \
    -o .

  echo "[done] wrote merged CSVs"
  """
}


process Heatmap_Virusseeker {
  tag { tagName }
  publishDir "${params.outdir}/${params.project_id}/VS_merged_arc_files/VS_Heatmaps/", mode: 'copy', overwrite: true,
             saveAs: { fn -> "${tagName}/${fn}" }   // dynamic subdir per input
  conda "$baseDir/env/aio_qc.yml"
  errorStrategy 'ignore'
  input:
  tuple val(tagName), path(csv)

  output:
  tuple val(tagName), path("*_vs_heatmap.png"),               emit: heatmap
  tuple val(tagName), path("*_vs_scaleByFamily_heatmap.png"), emit: heatmap_scaled_by_family
  tuple val(tagName), path("*_vs_scaleBySample_heatmap.png"), emit: heatmap_scaled_by_sample


  script:
  """
  #eval "\$(command conda 'shell.bash' 'hook' 2> /dev/null)"
  #conda activate qc_report
  # TMP FIX
  #source /opt/conda/etc/profile.d/conda.sh
  #conda activate all_in_one_pipeline
  Rscript ${params.scripts}/VS_heatmap.R \
    -r ${csv} \
    -g Family \
    -o ./
  """
}
