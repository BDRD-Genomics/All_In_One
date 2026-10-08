#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
  BlastX
========================================================================================
*/

process BlastX_contigs {
    tag { "${assembler} | ${sample_id}" }
    publishDir { "${params.outdir}/${params.run_id}/${sample_id}/blast/" }, mode: 'copy'
    label 'optimized_blastx_contigs'
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), val(assembler), path(contigs_fasta)

    output:
    tuple val(sample_id), val(assembler), file("${sample_id}_${assembler}_blastx.daa"), emit: blastx_contigs_ch

    when:
    (params.dragonflye || params.spades || params.metaspades || params.unicycler)

    script:
    """
    diamond blastx ${params.diamond_args} \\
        --threads ${task.cpus} \\
        --db ${params.diamond_dbdir} \\
        --query ${contigs_fasta} \\
        --range-culling -F 15 \\
        --evalue 1e-5 \\
        --outfmt 100 \\
        --out ${sample_id}_${assembler}_blastx.daa
    """
}

process BlastX_reads {
    tag { "${sample_id}_${mode}" }
    publishDir { "${params.outdir}/${params.run_id}/${sample_id}/blast/" }, mode: 'copy'
    label 'optimized_blastx_reads'
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), file(fastq_file), val(mode)

    output:
    tuple val(sample_id), file("*"), emit: blastx_reads_ch
    tuple val(sample_id), file("${sample_id}_short_reads_blastx.daa"), optional: true, emit: short_reads_blastx_channel
    tuple val(sample_id), file("${sample_id}_long_reads_blastx.daa"), optional: true, emit: long_reads_blastx_channel

    when:
    params.agnostic_read_analysis

    script:
    """
    # A gzip file containing zero reads is still non-zero bytes on disk, so -s alone
    # is not enough. Check that the decompressed/plain FASTQ has actual content.
    has_reads=false

    if [[ -s "${fastq_file}" ]]; then
        if [[ "${fastq_file}" == *.gz ]]; then
            if gzip -cd "${fastq_file}" 2>/dev/null | head -n 1 | grep -q .; then
                has_reads=true
            fi
        else
            if head -n 1 "${fastq_file}" | grep -q .; then
                has_reads=true
            fi
        fi
    fi

    if [[ "${mode}" == "long" && "\$has_reads" == "true" ]]; then
        diamond blastx ${params.diamond_args} \\
            --threads ${task.cpus} \\
            --db ${params.diamond_dbdir} \\
            --query "${fastq_file}" \\
            --evalue 1e-5 \\
            --outfmt 100 \\
            --out "${sample_id}_long_reads_blastx.daa"

    elif [[ "${mode}" == "short" && "\$has_reads" == "true" ]]; then
        diamond blastx ${params.diamond_args} \\
            --threads ${task.cpus} \\
            --db ${params.diamond_dbdir} \\
            --query "${fastq_file}" \\
            --evalue 1e-5 \\
            --outfmt 100 \\
            --out "${sample_id}_short_reads_blastx.daa"

    else
        echo "SKIP: No valid reads provided for ${mode} mode" > "${sample_id}_blastx.skipped.txt"
    fi

    echo "OK" > "${sample_id}_blastx_reads.finished"
    """
}

// DAA2INFO contigs .daa file
process DAA2INFO_contigs_daa_file {
    tag { sample_id }
    errorStrategy 'ignore'
    publishDir { "${params.outdir}/${params.run_id}/${sample_id}/blast/" }, mode: 'copy'
    label 'normal'

    input:
    tuple val(sample_id), file(blastx_contigs)

    output:
    tuple val(sample_id), file("${sample_id}_c2c.txt"), emit: c2c_txt_file

    script:
    """
    daa2info \\
        -i ${blastx_contigs} \\
        -o ${sample_id}_c2c.txt \\
        -c2c Taxonomy -n -r -u
    """
}

/* Meganize Blastx Reads */

process Meganize_ShortReads_BlastX {
    tag { sample_id }
    errorStrategy 'ignore'
    publishDir { "${params.outdir}/${params.run_id}/${sample_id}/blast/meganized_reads" }, mode: 'copy'
    label 'megan'

    input:
    tuple val(sample_id), file(blastx_short)

    output:
    tuple val(sample_id), file("${sample_id}_short_reads_blastx_daa_summary_count.tsv"), emit: meganize_short_reads_ch
    tuple val(sample_id), file("*"), emit: all_shortreads_meganized_files
    tuple val(sample_id), file(blastx_short), emit: meganized_short_daa_ch

    when:
    params.agnostic_read_analysis || params.shortreads || params.hybrid

    script:
    """
    daa-meganizer \\
        --in ${blastx_short} \\
        --only Taxonomy \\
        --mapDB ${params.megan_mdb} \\
        --threads ${task.cpus} \\
        --minSupportPercent 0 \\
        --topPercent 0.5 \\
        --lcaAlgorithm weighted \\
        --longReads false \\
        --verbose

    paste <(daa2info -i ${blastx_short} -c2c Taxonomy | awk '{print \$1}') \\
          <(daa2info -i ${blastx_short} -p -c2c Taxonomy | awk '{print \$1,\$2}' FS='\t' OFS='\t') \\
          > ${sample_id}_short_reads_blastx_daa_summary_count.tsv

    daa2info \\
        -i ${blastx_short} \\
        -c2c Taxonomy \\
        -o ${sample_id}_shortreads_blastx.tsv

    ktImportTaxonomy \\
        -tax ${params.krona_db} \\
        -t 1 \\
        -m 2 \\
        ${sample_id}_shortreads_blastx.tsv \\
        -o ${sample_id}_shortreads_krona.html
    """
}

process Meganize_LongReads_BlastX {
    tag { sample_id }
    errorStrategy 'ignore'
    publishDir { "${params.outdir}/${params.run_id}/${sample_id}/blast/meganized_reads" }, mode: 'copy'
    label 'megan'

    input:
    tuple val(sample_id), file(blastx_long)

    output:
    tuple val(sample_id), file("${sample_id}_long_reads_blastx_daa_summary_count.tsv"), emit: meganize_long_reads_ch
    tuple val(sample_id), file("*"), emit: all_longreads_meganized_files
    tuple val(sample_id), file(blastx_long), emit: long_reads_daa_meganized_ch

    when:
    params.agnostic_read_analysis || params.longreads || params.hybrid

    script:
    """
    daa-meganizer \\
        --in ${blastx_long} \\
        --only Taxonomy \\
        --mapDB ${params.megan_mdb} \\
        --threads ${task.cpus} \\
        --minSupportPercent 0 \\
        --topPercent 0.5 \\
        --lcaAlgorithm weighted \\
        --longReads true \\
        --verbose

    paste <(daa2info -i ${blastx_long} -c2c Taxonomy | awk '{print \$1}') \\
          <(daa2info -i ${blastx_long} -p -c2c Taxonomy | awk '{print \$1,\$2}' FS='\t' OFS='\t') \\
          > ${sample_id}_long_reads_blastx_daa_summary_count.tsv

    daa2info \\
        -i ${blastx_long} \\
        -c2c Taxonomy \\
        -o ${sample_id}_longreads_blastx.tsv

    ktImportTaxonomy \\
        -tax ${params.krona_db} \\
        -t 1 \\
        -m 2 \\
        ${sample_id}_longreads_blastx.tsv \\
        -o ${sample_id}_longreads_blastx_krona.html
    """
}

process Meganize_BlastX_Contigs {
    tag { "${assembler} | ${sample_id}" }
    publishDir { "${params.outdir}/${params.run_id}/${sample_id}/blast/meganized_contigs/" }, mode: 'copy'
    label 'optimized_Meganize_Contigs_BlastX'
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), val(assembler), path(daa_file)

    output:
    tuple val(sample_id), val(assembler), file("${sample_id}_${assembler}_contigs_blastx_diamondview.tsv"), emit: meganized_blastx_contigs_ch
    tuple val(sample_id), val(assembler), file("${sample_id}_${assembler}_contigs_blastx_daa_summary_count.tsv"), emit: daa_summary_ch
    tuple val(sample_id), val(assembler), file("${sample_id}_${assembler}_contigs_blastx.tsv"), emit: megan_tax_tsv_ch
    tuple val(sample_id), val(assembler), file("${sample_id}_${assembler}_contigs_blastx_krona.html"), optional: true, emit: krona_html_ch
    tuple val(sample_id), val(assembler), file("${sample_id}_${assembler}_blastx.daa"), optional: true, emit: krona_daa_ch

    script:
    """
    daa-meganizer \\
        --in ${daa_file} \\
        --only Taxonomy \\
        --mapDB ${params.megan_mdb} \\
        --threads ${task.cpus} \\
        --topPercent 0.5 \\
        --minSupportPercent 0 \\
        --lcaAlgorithm longReads \\
        --longReads true \\
        --verbose

    # DAA summary counts (ID + counts)
    paste <(daa2info -i ${daa_file} -c2c Taxonomy | awk '{print \$1}') \\
          <(daa2info -i ${daa_file} -p -c2c Taxonomy | awk -F'\t' '{print \$1, \$2}' OFS='\t') \\
          > ${sample_id}_${assembler}_contigs_blastx_daa_summary_count.tsv

    diamond view \\
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle qlen slen \\
        --daa ${daa_file} \\
        > ${sample_id}_${assembler}_contigs_blastx_diamondview.tsv

    # MEGAN taxonomy TSV (Krona input)
    daa2info \\
        -i ${daa_file} \\
        -c2c Taxonomy \\
        -o ${sample_id}_${assembler}_contigs_blastx.tsv

    # Krona
    ktImportTaxonomy \\
        -tax ${params.krona_db} \\
        -t 1 \\
        -m 2 \\
        ${sample_id}_${assembler}_contigs_blastx.tsv \\
        -o ${sample_id}_${assembler}_contigs_blastx_krona.html
    """
}

process Parse_BlastX_Contigs {
    tag { "${assembler} | ${sample_id}" }
    publishDir { "${params.outdir}/${params.run_id}/${sample_id}/blast/" }, mode: 'copy'
    label 'lowmem'
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), val(assembler), val(diamondview_file)

    output:
    tuple val(sample_id), val(assembler), file("*"), emit: parse_blastx_ch

    script:
    """
    # TMP FIX
    export OMP_NUM_THREADS=${task.cpus}
    export OPENBLAS_NUM_THREADS=${task.cpus}
    export MKL_NUM_THREADS=${task.cpus}
    export NUMEXPR_NUM_THREADS=${task.cpus}
    export VECLIB_MAXIMUM_THREADS=${task.cpus}

    python ${params.scripts}/VS_MD_diamond_parser_linFilt_Mar2026_fast.py \\
        -i ${diamondview_file} \\
        -t blastx \\
        -v ${params.megan_mdb} \\
        -n ${params.ncbi_taxa}
    """
}
