#!/usr/bin/ nextflow


nextflow.enable.dsl=2



/*
========================================================================================
   Pull Unassigned MMseqs Contigs
========================================================================================

*/

/* map reads to reference sequence and use for downstream analysis (Optional)*/

process Pull_Unassigned_MMseqs_Blastx {
    tag { "${sample_id}_${assembler}" }
    publishDir { "${params.outdir}/${params.run_id}/${sample_id}/VS_supplemental_outputs" },
                mode: 'copy',
                overwrite: true,
                saveAs: { fn -> "${assembler}/${fn}" }
    errorStrategy 'ignore'
    conda "$baseDir/env/md.yaml"
    cpus { cpus }
    memory { mem }

    input:
    tuple val(sample_id), path(contigs), path(mmseqs_parsed), path(diamondview_contigs), val(mode), val(cpus), val(mem), val(assembler)

    output:
    tuple val(sample_id), file("${sample_id}_${assembler}_contigs_blastx_unassignedMMseqs.parsed"), optional: true, emit: unassigned_mmseqs_parsed_ch

    when:
    params.vs

    script:
    """
    awk -F',' '{print \$2}' ${mmseqs_parsed} | tail -n +2 | uniq > mmseqs_assigned_headers.txt
    seqkit grep -f mmseqs_assigned_headers.txt -v ${contigs} | seqkit seq -ni > unassigned_headers.txt
    if [[ -s unassigned_headers.txt ]];then
        head -1 ${diamondview_contigs} > ${sample_id}_${assembler}_contigs_blastx_unassignedMMseqs.parsed
        grep -f unassigned_headers.txt ${diamondview_contigs} >> ${sample_id}_${assembler}_contigs_blastx_unassignedMMseqs.parsed || > ${sample_id}_${assembler}_contigs_blastx_unassignedMMseqs.parsed
    else
        touch ${sample_id}_${assembler}_contigs_blastx_unassignedMMseqs.parsed
    fi
    """
}

/*
========================================================================================
   Pull Unassigned MMseqs Reads
========================================================================================

*/

/* map reads to reference sequence and use for downstream analysis (Optional)*/

process Pull_Unassigned_MMseqs_Reads {
    tag { "${sample_id}" }
    publishDir { "${params.outdir}/${params.run_id}/${sample_id}/VS_supplemental_outputs" },
                mode: 'copy',
                overwrite: true,
                saveAs: { fn -> "${assembler}/${fn}" }
    //label 'optimized_map2reads' // can erase if new tuple works
    errorStrategy 'ignore'
    conda "$baseDir/env/md.yaml"
    cpus { cpus }
    memory { mem }

    input:
    tuple val(sample_id),
	  val(read_type),
	  path(mmseqs_fasta),
	  path(mmseqs_parsed),
	  path(diamondview),
	  val(mode),
	  val(cpus),
          val(mem),
          val(assembler)

    output:
    tuple val(sample_id),
	  val(read_type),
	  val(assembler),
	  val(mode),
	  path("${sample_id}_${read_type}_blastx_unassignedMMseqs.out"),
	  emit: unassigned_mmseqs_blastx_ch

    when:
    params.vs

    script:
    """
    if [[ "${mode}" == "short" || "${mode}" == "hybrid" ]]; then
        awk -F',' '{print \$2}' ${mmseqs_parsed} | tail -n +2 | uniq > mmseqs_sr_assigned_headers.txt
        seqkit grep -f mmseqs_sr_assigned_headers.txt -v ${mmseqs_fasta} | seqkit seq -ni > sr_unassigned_headers.txt
        if [[ -s sr_unassigned_headers.txt ]];then
            head -1 ${diamondview} > ${sample_id}_${read_type}_blastx_unassignedMMseqs.out
            grep -f sr_unassigned_headers.txt ${diamondview} >> ${sample_id}_${read_type}_blastx_unassignedMMseqs.out || > ${sample_id}_${read_type}_blastx_unassignedMMseqs.out
        else
            touch ${sample_id}_${read_type}_blastx_unassignedMMseqs.out
        fi
    fi
    if [[ "${mode}" == "long" || "${mode}" == "hybrid" ]]; then
        awk -F',' '{print \$2}' ${mmseqs_parsed} | tail -n +2 | uniq > mmseqs_lr_assigned_headers.txt
        seqkit grep -f mmseqs_lr_assigned_headers.txt -v ${mmseqs_fasta} | seqkit seq -ni > lr_unassigned_headers.txt
	if [[ -s lr_unassigned_headers.txt ]];then
            head -1 ${diamondview} > ${sample_id}_${read_type}_blastx_unassignedMMseqs.out
            grep -f lr_unassigned_headers.txt ${diamondview} >> ${sample_id}_${read_type}_blastx_unassignedMMseqs.out || > ${sample_id}_${read_type}_blastx_unassignedMMseqs.out
        else
            touch ${sample_id}_${read_type}_blastx_unassignedMMseqs.out
        fi
    fi
    """
}


/*
========================================================================================
   Create diamondview for Reads
========================================================================================

*/

/* map reads to reference sequence and use for downstream analysis (Optional)*/

process Diamondview {
    tag { "${sample_id}_reads" }
    publishDir { "${params.outdir}/${params.run_id}/${sample_id}/blast/" },
                mode: 'copy',
                overwrite: true,
                saveAs: { fn -> "${fn}" }
    //label 'optimized_map2reads' // can erase if new tuple works
    errorStrategy 'ignore'
    conda "$baseDir/env/md.yaml"
    cpus { cpus }
    memory { mem }

    input:
    tuple val(sample_id),
          val(read_type),
          path(blastx_daa),
          val(mode),
          val(cpus),
          val(mem),
          val(assembler)

    output:
    tuple val(sample_id),
          val(read_type),
          val(assembler),
          val(mode),
          path("*_${read_type}reads_blastx_diamondview.tsv"),
          emit: diamondview_reads_ch

    when:
    params.vs

    script:
    """
    ln -sf ${blastx_daa} input_${read_type}.daa

    diamond view \\
        --threads ${task.cpus} \\
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle qlen slen \\
        --daa input_${read_type}.daa \\
        > ${sample_id}_${read_type}reads_blastx_diamondview.tsv
    """
}


/*
========================================================================================
   Pull Unmapped Short Reads
========================================================================================

*/

/* map reads to reference sequence and use for downstream analysis (Optional)*/

process Pull_and_Parse_Unmapped_Reads {
    tag { "${sample_id}_${assembler}" }
    publishDir { "${params.outdir}/${params.run_id}/${sample_id}/VS_supplemental_outputs/" },
                mode: 'copy',
                overwrite: true,
                saveAs: { fn -> "${assembler}/${fn}" }
    errorStrategy 'ignore'
    conda "$baseDir/env/vs.yml"
    label 'lowmem'
    input:
    tuple val(sample_id),
	  val(read_type),
          path(unassigned_blastx),
          val(mode),
          val(assembler)

    output:
    tuple val(sample_id),
	  val(read_type),
	  val(assembler),
          val(mode),
          path("*_blastx_unassignedMMseqs.parsed"),
	  emit: parsed_diamond_all_ch
    when:
    params.vs

    script:
    """
    export OMP_NUM_THREADS=${task.cpus}
    export OPENBLAS_NUM_THREADS=${task.cpus}
    export MKL_NUM_THREADS=${task.cpus}
    export NUMEXPR_NUM_THREADS=${task.cpus}

    out="${sample_id}_${read_type}_blastx_unassignedMMseqs.parsed"

    if [[ -s ${unassigned_blastx} ]]; then
        python ${params.scripts}/VS_MD_diamond_parser_linFilt_Mar2026_fast.py \\
            -i ${unassigned_blastx} \\
	    -s ${sample_id} \\
            -o "\$out" \\
            -t blastx \\
            -v ${params.vhunter} \\
            -n ${params.ncbi_taxa} \\
	    --no-sort
    else
	touch "\$out"
    fi
    """
}

/*
========================================================================================
   Fastq 2 Fasta
========================================================================================

*/

/* map reads to reference sequence and use for downstream analysis (Optional)*/

process Fastq_2_Fasta {
    tag { "${sample_id}_reads" }
    publishDir { "${params.outdir}/${params.run_id}/${sample_id}/VS_supplemental_outputs" },
                mode: 'copy',
                overwrite: true,
                saveAs: { fn -> "${assembler}/${fn}" }
    //label 'optimized_map2reads' // can erase if new tuple works
    errorStrategy 'ignore'
    conda "$baseDir/env/md.yaml"
    cpus { cpus }
    memory { mem }

    input:
    tuple val(sample_id), file(fastq_unmapped_sr), file(fastq_unmapped_lr), val(cpus), val(mem), val(mode), val(assembler)

    output:
    tuple val(sample_id), file("${sample_id}_assembly_unmapped_SR.fasta"), optional: true, emit: fastq2fasta_ch_sr
    tuple val(sample_id), file("${sample_id}_assembly_unmapped_LR.fasta"), optional: true, emit: fastq2fasta_ch_lr
    tuple val(sample_id), file("${sample_id}_assembly_unmapped_*.fasta"), optional: true, emit: fastq2fasta_ch_all

    when:
    params.vs

    script:
    """
    if [[ ${mode} == "short" || ${mode} == "hybrid" ]]; then
        seqkit fq2fa -j ${cpus} ${fastq_unmapped_sr} > ${sample_id}_assembly_unmapped_SR.fasta
    fi
    if [[ ${mode} == "long" || ${mode} == "hybrid" ]]; then
        seqkit fq2fa -j ${cpus} ${fastq_unmapped_lr} > ${sample_id}_assembly_unmapped_LR.fasta
    fi
    """
}


/*
========================================================================================
   Map 2 Contigs
========================================================================================

*/

/* map reads to reference sequence and use for downstream analysis (Optional)*/

process Map_Reads_2_Contigs {
    tag { "${sample_id}_${assembler}" }
    publishDir { "${params.outdir}/${params.run_id}/${sample_id}/VS_supplemental_outputs" }, 
                mode: 'copy',
                overwrite: true,
                saveAs: { fn -> "${assembler}/${fn}" }
    //label 'optimized_map2reads' // can erase if new tuple works
    errorStrategy 'ignore'
    conda "$baseDir/env/md.yaml"
    cpus { cpus }
    memory { mem }

    input:
    tuple val(sample_id), file(fastq_1), file(fastq_2), file(long_read), file(contigs), val(mode), val(cpus), val(mem), val(assembler)

    output:
    tuple val(sample_id), file("${sample_id}_assembly_mapping_LR_covstats.txt"), file("lr_count.txt"), optional: true, emit: map2assembly_ch_lr
    tuple val(sample_id), file("${sample_id}_assembly_mapping_SR_covstats.txt"), file("sr_count.txt"), optional: true, emit: map2assembly_ch_sr
    tuple val(sample_id), file("${sample_id}_assembly_mapped_LR.fastq.gz"), file("${sample_id}_assembly_unmapped_LR.fastq.gz"), optional: true, emit: mapping_reads2contigs_LR
    tuple val(sample_id), file("${sample_id}_assembly_mapped_SR.fastq.gz"), file("${sample_id}_assembly_unmapped_SR.fastq.gz"),	optional: true,	emit: mapping_reads2contigs_SR
    tuple val(sample_id), file("*"), optional: true, emit: map2assembly_ch_all


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
    #mkdir -p ${params.outdir}/${params.run_id}/${sample_id}/status_log/

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
    fi
    """
}


/*
========================================================================================
   Parse_MMseqs_Reads
========================================================================================

*/


process Parse_MMseqs_Reads {
  tag { "${sample_id} | ${assembler} | ${read_type}" }
  publishDir { "${params.outdir}/${params.run_id}/${sample_id}/VS_supplemental_outputs/${assembler}" },
              mode: 'copy',
              overwrite: true,
              saveAs: { fn -> "${assembler}/${fn}" }

  label 'lowmem'
  errorStrategy 'ignore'
  conda "${baseDir}/env/vs.yml"

  input:
  tuple val(sample_id), 
        val(assembler),
        val(read_type),
        path(mmseqs_out)

  output:
  tuple val(sample_id),
	val(assembler),
        val(read_type),
        path("${sample_id}_${assembler}_${read_type}.mmseqs.parsed"),
        emit: mmseqs_parsed_ch
script:
  """
  # TMP FIX
  #export OMP_NUM_THREADS=${task.cpus}
  #export OPENBLAS_NUM_THREADS=${task.cpus}
  #export MKL_NUM_THREADS=${task.cpus}
  #export NUMEXPR_NUM_THREADS=${task.cpus}
  # remove header
  tail -n +2 "${mmseqs_out}" > tmp && mv tmp "${mmseqs_out}"
  python ${params.scripts}/VS_MD_diamond_parser_linFilt_Mar2026_fast.py \\
    -i "${mmseqs_out}" \\
    -s ${sample_id}_${assembler}_${read_type} \\
    -t mmseqs \\
    -v ${params.vhunter} \\
    -n ${params.ncbi_taxa} \\
    --no-sort
  """
}


/*
process Parse_MMseqs_Reads {
  tag { "${sample_id}" }
  publishDir { "${params.outdir}/${params.run_id}/${sample_id}/VS_supplemental_outputs/" },
              mode: 'copy',
              overwrite: true,
              saveAs: { fn -> "${assembler}/${fn}" }

  label 'lowmem'
  errorStrategy 'ignore'
  conda "${baseDir}/env/vs.yml"

  input:
  tuple val(sample_id), file(mmseq_out_sr), file(mmseq_out_lr), val(cpus), val(mem), val(mode),val(assembler)

  cpus {cpus}  // setting slurm allocation dynamically
  memory {mem} // setting slurm allocation dynamically

  output:
  tuple val(sample_id), file("${sample_id}_sr.mmseqs.parsed"), optional: true, emit: mmseqs_sr_parsed_ch
  tuple val(sample_id), file("${sample_id}_lr.mmseqs.parsed"), optional: true, emit: mmseqs_lr_parsed_ch
  tuple val(sample_id), file("${sample_id}_*.mmseqs.parsed"), optional: true, emit: mmseqs_all_parsed_ch


script:
  """
  set -euo pipefail
  # TMP FIX
  export OMP_NUM_THREADS=${task.cpus}
  export OPENBLAS_NUM_THREADS=${task.cpus}
  export MKL_NUM_THREADS=${task.cpus}
  export NUMEXPR_NUM_THREADS=${task.cpus}
  if [[ "${mode}" == "short" || "${mode}" == "hybrid" ]]; then
    # remove header
    tail -n +2 "${mmseq_out_sr}" > tmp && mv tmp "${mmseq_out_sr}"

    python ${params.scripts}/VS_MD_diamond_parser_linFilt_Mar2026_fast.py \\
      -i "${mmseq_out_sr}" \\
      -s ${sample_id}_sr \\
      -t mmseqs \\
      -v ${params.vhunter} \\
      -n ${params.ncbi_taxa} \\
      --no-sort
  fi

  if [[ "${mode}" == "long" || "${mode}" == "hybrid" ]]; then
    # remove header
    tail -n +2 "${mmseq_out_lr}" > tmp && mv tmp "${mmseq_out_lr}"

    python ${params.scripts}/VS_MD_diamond_parser_linFilt_Mar2026_fast.py \\
      -i "${mmseq_out_lr}" \\
      -s ${sample_id}_lr \\
      -t mmseqs \\
      -v ${params.vhunter} \\
      -n ${params.ncbi_taxa} \\
      --no-sort
  fi

  """
}
*/


/*
========================================================================================
   Bash Filter Viruses
========================================================================================
*/

/* filter mmseqs for "high-confidence" viral calls and output filtered blast tables with associated fasta files and a final semi-quantitative read counts file*/

process Filter_Putative_Viruses_ALL {
    tag { "${sample_id}_${assembler}_${mode}" }

    publishDir { "${params.outdir}/${params.run_id}/${sample_id}/VS_supplemental_outputs" },
        mode: 'copy',
        overwrite: true,
        saveAs: { fn -> "${assembler}/${fn}" }

    label 'optimized_map2reads'
    errorStrategy 'ignore'

    input:
    tuple val(sample_id),
          path(mmseqs_contigs_parsed),
          path(diamond_contigs_parsed),
          val(mmseqs_short_reads_parsed),
          val(mmseqs_long_reads_parsed),
          val(diamond_short_reads_parsed),
          val(diamond_long_reads_parsed),
          val(assembler),
          val(mode)

    output:
    tuple val(sample_id),
          path("putative_viral_blast.csv"),
          val(assembler),
          optional: true,
          emit: viral_blast_tab_ch

    when:
    params.vs

    script:
    """
    data_count=0    

    head -1 ${mmseqs_contigs_parsed} > ${sample_id}_classified_parsed_blast.csv 

    add_if_nonempty () {
        f="\$1" 

        # Skip null-like values coming from optional short/long mode fields
        if [[ -z "\$f" || "\$f" == "null" || "\$f" == "[]" ]]; then
            return 0
        fi  

        # Skip anything that is not a regular file.
        # This prevents: wc: 'standard input': Is a directory
        if [[ ! -f "\$f" ]]; then
            echo "Skipping non-file input: \$f"
            return 0
        fi  

        # Skip empty files
        if [[ ! -s "\$f" ]]; then
            echo "Skipping empty file: \$f"
            return 0
        fi  

        line_count=\$(wc -l < "\$f")    

        if (( line_count > 1 )); then
            tail -n +2 "\$f" >> ${sample_id}_classified_parsed_blast.csv
            data_count=\$((data_count + 1))
        else
            echo "Skipping header-only file: \$f"
        fi
    }   

    if [[ "${mode}" == "short" ]]; then
        add_if_nonempty "${mmseqs_contigs_parsed}"
        add_if_nonempty "${diamond_contigs_parsed}"
        add_if_nonempty "${mmseqs_short_reads_parsed}"
        add_if_nonempty "${diamond_short_reads_parsed}"
    fi  

    if [[ "${mode}" == "long" ]]; then
        add_if_nonempty "${mmseqs_contigs_parsed}"
        add_if_nonempty "${diamond_contigs_parsed}"
        add_if_nonempty "${mmseqs_long_reads_parsed}"
        add_if_nonempty "${diamond_long_reads_parsed}"
    fi  

    if [[ "${mode}" == "hybrid" ]]; then
        add_if_nonempty "${mmseqs_contigs_parsed}"
        add_if_nonempty "${diamond_contigs_parsed}"
        add_if_nonempty "${mmseqs_short_reads_parsed}"
        add_if_nonempty "${mmseqs_long_reads_parsed}"
        add_if_nonempty "${diamond_short_reads_parsed}"
        add_if_nonempty "${diamond_long_reads_parsed}"
    fi  

    if (( data_count > 0 )); then
        bash ${params.scripts}/commands_to_filter_blast.sh ${sample_id}_classified_parsed_blast.csv || touch putative_viral_blast.csv
    else
        touch putative_viral_blast.csv
    fi  

    [[ -f putative_viral_blast.csv ]] || touch putative_viral_blast.csv
    """
}



process Generate_VS_Outputs {
    tag { "${sample_id}_${assembler}" }

    publishDir { "${params.outdir}/${params.run_id}/${sample_id}/VS_supplemental_outputs" },
        mode: 'copy',
        overwrite: true,
        saveAs: { fn -> "${assembler}/${fn}" }

    label 'optimized_map2reads'
    errorStrategy 'ignore'
    conda "$baseDir/env/vs.yml"

    input:
    tuple val(sample_id),
          path(putative_viral_output),
          path(contigs_fasta),
          val(unmapped_short_reads_fasta),
          val(unmapped_long_reads_fasta),
          val(PE_covstats),
          val(PE_QC_reads),
          val(LR_covstats),
          val(LR_QC_reads),
          val(mode),
          val(assembler)
    output:
    tuple val(sample_id), path("*"), val(assembler), optional: true, emit: vs_reports_ch
    tuple val(sample_id), path("${sample_id}_accurate_read_counts.tsv"), val(assembler), optional: true, emit: accurate_read_counts_ch
    tuple val(assembler), path("${sample_id}_accurate_read_counts.tsv"), optional: true, emit: accurate_read_counts_file_ch

    when:
    params.vs

    script:
    def pe_cov_arg = PE_covstats ? "-p ${PE_covstats}" : ""
    def pe_qc_arg  = PE_QC_reads ? "-S ${PE_QC_reads}" : ""
    def lr_cov_arg = LR_covstats ? "-l ${LR_covstats}" : ""
    def lr_qc_arg  = LR_QC_reads ? "-L ${LR_QC_reads}" : ""
    def bad_acc_arg = params.bad_accessions_list ? "-b ${params.bad_accessions_list}" : ""

    """
    if [[ "${mode}" == "short" ]]; then
        [[ -n "${unmapped_short_reads_fasta}" ]] || { echo "ERROR: missing unmapped_short_reads_fasta"; exit 1; }
        [[ -n "${PE_covstats}" ]] || { echo "ERROR: missing PE_covstats"; exit 1; }
        [[ -n "${PE_QC_reads}" ]] || { echo "ERROR: missing PE_QC_reads / -S input"; exit 1; }
    fi  

    if [[ "${mode}" == "long" ]]; then
        [[ -n "${unmapped_long_reads_fasta}" ]] || { echo "ERROR: missing unmapped_long_reads_fasta"; exit 1; }
        [[ -n "${LR_covstats}" ]] || { echo "ERROR: missing LR_covstats"; exit 1; }
        [[ -n "${LR_QC_reads}" ]] || { echo "ERROR: missing LR_QC_reads / -L input"; exit 1; }
    fi  

    if [[ "${mode}" == "hybrid" ]]; then
        [[ -n "${unmapped_short_reads_fasta}" ]] || { echo "ERROR: missing unmapped_short_reads_fasta"; exit 1; }
        [[ -n "${unmapped_long_reads_fasta}" ]] || { echo "ERROR: missing unmapped_long_reads_fasta"; exit 1; }
        [[ -n "${PE_covstats}" ]] || { echo "ERROR: missing PE_covstats"; exit 1; }
        [[ -n "${PE_QC_reads}" ]] || { echo "ERROR: missing PE_QC_reads / -S input"; exit 1; }
        [[ -n "${LR_covstats}" ]] || { echo "ERROR: missing LR_covstats"; exit 1; }
        [[ -n "${LR_QC_reads}" ]] || { echo "ERROR: missing LR_QC_reads / -L input"; exit 1; }
    fi

    if [[ "${mode}" == "short" ]]; then
        cat ${contigs_fasta} ${unmapped_short_reads_fasta} > contigs_and_reads.fasta

        python ${params.scripts}/filter_blast_viruses_parserDefCols.py \\
            --sample-id ${sample_id} \\
            -t ${putative_viral_output} \\
            -c contigs_and_reads.fasta \\
            ${pe_cov_arg} \\
            ${pe_qc_arg} \\
            -n ${params.ncbi_taxa} \\
            -v ${params.virus_genome_sizes} \\
            ${bad_acc_arg}
    fi

    if [[ "${mode}" == "long" ]]; then
        cat ${contigs_fasta} ${unmapped_long_reads_fasta} > contigs_and_reads.fasta

        python ${params.scripts}/filter_blast_viruses_parserDefCols.py \\
            --sample-id ${sample_id} \\
            -t ${putative_viral_output} \\
            -c ${contigs_fasta} \\
            ${lr_cov_arg} \\
            ${lr_qc_arg} \\
            -n ${params.ncbi_taxa} \\
            -v ${params.virus_genome_sizes} \\
            ${bad_acc_arg}
    fi

    if [[ "${mode}" == "hybrid" ]]; then
        cat ${contigs_fasta} ${unmapped_short_reads_fasta} ${unmapped_long_reads_fasta} > contigs_and_reads.fasta

        python ${params.scripts}/filter_blast_viruses_parserDefCols.py \\
            --sample-id ${sample_id} \\
            -t ${putative_viral_output} \\
            -c ${contigs_fasta} \\
            ${pe_cov_arg} \\
            ${pe_qc_arg} \\
            ${lr_cov_arg} \\
            ${lr_qc_arg} \\
            -n ${params.ncbi_taxa} \\
            -v ${params.virus_genome_sizes} \\
            ${bad_acc_arg}
    fi
    """
}


process Merge_Arc {
  tag { "merge_all_${assembler}" }
  publishDir "${params.outdir}/${params.run_id}/VS_merged_arc_files/",
              mode: 'copy',
              overwrite: true,
              saveAs: { fn -> "${assembler}/${fn}" }  
  conda "$baseDir/env/vs.yml"
  input:
  //path accurate_read_counts_dir  // a list of all matched TSVs
  tuple val(assembler), path(arc_files)
  output:
  //path "${params.run_id}_total_hits_merged.csv"
  //path "${params.run_id}_normalized_family_hits_merged.csv"
  //path "${params.run_id}_normalized_RPM_hits_merged.csv"
  path "${params.run_id}_${assembler}_total_hits_merged.csv",              emit: total_csv
  path "${params.run_id}_${assembler}_normalized_family_hits_merged.csv",  emit: nfam_csv
  path "${params.run_id}_${assembler}_normalized_RPM_hits_merged.csv",     emit: nrpm_csv

  script:
  """
  set -euo pipefail
  rm -rf arc_${assembler}
  mkdir -p arc_${assembler}
  cp -f ${arc_files} arc_${assembler}
  #printf '%s\0' ${arc_files} | xargs -0 -I '{}' cp -L -n '{}' "arc_${assembler}/"
  python ${params.scripts}/combine_ARC_hits.py \
    -r arc_${assembler} \
    -p ${params.run_id}_${assembler} \
    -o .

  echo "[done] wrote merged CSVs"
  """
}


process Heatmap_Virusseeker {
  tag { tagName }
  publishDir "${params.outdir}/${params.run_id}/VS_merged_arc_files/VS_Heatmaps/", mode: 'copy', overwrite: true,
             saveAs: { fn -> "${tagName}/${fn}" }   // dynamic subdir per input
  conda "$baseDir/env/aio_qc.yml"

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

  Rscript ${params.scripts}/VS_heatmap.R \
    -r ${csv} \
    -g Family \
    -o ./
  """
}
