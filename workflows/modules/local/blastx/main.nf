#!/usr/bin/ nextflow

nextflow.enable.dsl=2

/*
========================================================================================
   BlastX
========================================================================================
   Github   : 
   Contact  :     
----------------------------------------------------------------------------------------

*/

process BlastX_contigs {
    tag { "${assembler} | ${sample_id}" }
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/blast/", mode: 'copy'
    label 'optimized_blastx_contigs'
    conda "${baseDir}/env/diamond.yml"
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
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/blast/", mode: 'copy'
    label 'optimized_blastx_reads'
    conda "${baseDir}/env/diamond.yml"
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), val(fastq_file), val(mode)

    output:
    tuple val(sample_id), file("*"), emit: blastx_reads_ch
    tuple val(sample_id), file("${sample_id}_short_reads_blastx.daa"), optional: true, emit: short_reads_blastx_channel
    tuple val(sample_id), file("${sample_id}_long_reads_blastx.daa"), optional: true, emit: long_reads_blastx_channel
    
    when:
    params.agnostic_read_analysis

    script:
    """
    mkdir -p logs

    if [[ "$mode" == "long" && -s "$fastq_file" ]]; then
        echo "Running diamond on long read: $fastq_file" | tee logs/blastx_long.log
        diamond blastx ${params.diamond_args} \\
            --threads ${task.cpus} \\
            --db ${params.diamond_dbdir} \\
            --query "$fastq_file" \\
            --evalue 1e-5 \\
            --outfmt 100 \\
            --out "${sample_id}_long_reads_blastx.daa" >> BlastX_longreads.log 2>&1
    elif [[ "$mode" == "short" && -s "$fastq_file" ]]; then 
        echo "Running diamond on short reads: $fastq_file" | tee logs/blastx_short.log
        diamond blastx ${params.diamond_args} \\
            --threads ${task.cpus} \\
            --db ${params.diamond_dbdir} \\
            --query "$fastq_file" \\
            --evalue 1e-5 \\
            --outfmt 100 \\
            --out "${sample_id}_short_reads_blastx.daa" >> BlastX_shortreads.log 2>&1
    else
        echo "SKIP: No valid reads provided for $mode mode" > "${sample_id}_blastx.skipped.txt"
    fi

    echo "OK" > "${sample_id}_blastx_reads.finished"
    """
}

// DAA2INFO contigs .daa file

process DAA2INFO_contigs_daa_file {
    tag {sample_id}
    errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/blast/", mode: 'copy'
    conda "${baseDir}/env/diamond.yml"
    label 'normal'

    input:

    tuple val(sample_id), file(blastx_contigs)

    output:

    tuple val(sample_id), file("${sample_id}_c2c.txt"), emit: c2c_txt_file


    script:

    """
    if [ "${params.metaspades}" == "true" ]; then
      bash ${params.meganpath}/daa2info \
      -i ${params.outdir}/${params.project_id}/${sample_id}/blast/${sample_id}_metaspades_blastx.daa \
      -o ${sample_id}_c2c.txt \
      -c2c Taxonomy -n -r -u

    elif [ "${params.hybrid}" == "true" ]; then 
      bash ${params.meganpath}/daa2info \
      -i ${params.outdir}/${params.project_id}/${sample_id}/blast/${sample_id}_dragonflye_blastx.daa \
      -o ${sample_id}_c2c.txt \
      -c2c Taxonomy -n -r -u
    elif [ "${params.dragonflye_isolate}" == "true" ]; then 
      bash ${params.meganpath}/daa2info \
      -i ${params.outdir}/${params.project_id}/${sample_id}/blast/${sample_id}_dragonflye_isolate_contigs_blastx.daa \
      -o ${sample_id}_c2c.txt \
      -c2c Taxonomy -n -r -u
    elif [ "${params.unicycler}" == "true" ]; then 
      bash ${params.meganpath}/daa2info \
      -i ${params.outdir}/${params.project_id}/${sample_id}/blast/${sample_id}_unicycler_blastx.daa \
      -o ${sample_id}_c2c.txt \
      -c2c Taxonomy -n -r -u
    else 
      "skip this process"
  
    fi
    """
}

/* Meganize Blastx Reads */

process Meganize_ShortReads_BlastX {
    tag { sample_id }
    errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/blast/meganized_reads", mode: 'copy'
    label 'optimized_blastx_reads'
    conda "${baseDir}/env/diamond.yml"
    
    cpus { 64 }
    memory { '260 GB'}
    time '36h'

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
    ln -sf ${params.megandb}/ncbi.map ./ncbi.map 
    ln -sf ${params.megandb}/ncbi.tre ./ncbi.tre 

    ${params.meganpath}/daa-meganizer \\
        --in ${blastx_short} \\
        --only Taxonomy \\
        --mapDB ${params.megandb}/megan-map.db \\
        --threads ${task.cpus} \\
        --minSupportPercent 0 \\
        --topPercent 0.5 \\
        --lcaAlgorithm weighted \\
        --longReads false \\
        --verbose >> meganize_short_reads.log 2>&1

    paste <(${params.meganpath}/daa2info -P ${params.meganpath}/.MEGAN.def -i ${blastx_short} -c2c Taxonomy | awk '{print \$1}') \\
          <(${params.meganpath}/daa2info -P ${params.meganpath}/.MEGAN.def -i ${blastx_short} -p -c2c Taxonomy | awk '{print \$1,\$2}' FS='\\t' OFS='\\t') \\
          > ${sample_id}_short_reads_blastx_daa_summary_count.tsv

    ${params.meganpath}/daa2info -i ${blastx_short} \
        -c2c Taxonomy -o ${sample_id}_shortreads_blastx.tsv

    ktImportTaxonomy \
        -tax /export/database/krona/taxonomy \
        -t 1 -m 2 ${sample_id}_shortreads_blastx.tsv -o ${sample_id}_shortreads_krona.html
    """
}

process Meganize_LongReads_BlastX {
    tag { sample_id }
    errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/blast/meganized_reads", mode: 'copy'
    label 'optimized_blastx_reads'
    conda "${baseDir}/env/diamond.yml"
    cpus { 64 }
    memory { '260 GB'}
    time '36h'

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
    ln -sf ${params.megandb}/ncbi.map ./ncbi.map 
    ln -sf ${params.megandb}/ncbi.tre ./ncbi.tre 

    ${params.meganpath}/daa-meganizer \\
        --in ${blastx_long} \\
        --only Taxonomy \\
        --mapDB ${params.megandb}/megan-map.db \\
        --threads ${task.cpus} \\
        --minSupportPercent 0 \\
        --topPercent 0.5 \\
        --lcaAlgorithm weighted \\
        --longReads false \\
        --verbose >> meganize_longreads.log 2>&1

    paste <(${params.meganpath}/daa2info -P ${params.meganpath}/.MEGAN.def -i ${blastx_long} -c2c Taxonomy | awk '{print \$1}') \\
          <(${params.meganpath}/daa2info -P ${params.meganpath}/.MEGAN.def -i ${blastx_long} -p -c2c Taxonomy | awk '{print \$1,\$2}' FS='\\t' OFS='\\t') \\
          > ${sample_id}_long_reads_blastx_daa_summary_count.tsv

    ${params.meganpath}/daa2info -i ${blastx_long} \
        -c2c Taxonomy \
        -o ${sample_id}_longreads_blastx.tsv


    ktImportTaxonomy \
        -tax /export/database/krona/taxonomy \
        -t 1 -m 2 ${sample_id}_longreads_blastx.tsv -o ${sample_id}_longreads_blastx_krona.html

    """
}

process Meganize_BlastX_Contigs {
    tag { "${assembler} | ${sample_id}" }
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/blast/meganized_contigs/", mode: 'copy'
    label 'optimized_Meganize_Contigs_BlastX'
    conda "${baseDir}/env/diamond.yml"
    errorStrategy 'ignore'
    input:
    tuple val(sample_id), val(assembler), path(daa_file)

    output:
    tuple val(sample_id), val(assembler), file("${sample_id}_${assembler}_contigs_blastx_diamondview.tsv"), emit: meganized_blastx_contigs_ch
    tuple val(sample_id), val(assembler), file("${sample_id}_${assembler}_contigs_blastx_daa_summary_count.tsv"), emit: daa_summary_ch
    tuple val(sample_id), val(assembler), file("${sample_id}_${assembler}_contigs_blastx.tsv"), emit: megan_tax_tsv_ch
    tuple val(sample_id), val(assembler), file("${sample_id}_${assembler}_contigs_blastx_krona.html"), optional: true, emit: krona_html_ch
    tuple val(sample_id), val(assembler), file("${sample_id}_${assembler}_blastx.daa"), optional: true, emit: krona_daa_ch
    tuple val(sample_id), val(assembler), file("Meganize_BlastX_${assembler}.log"), optional: true, emit: log_ch

    script:
    """
    ln -sf ${params.megandb}/ncbi.map ./ncbi.map
    ln -sf ${params.megandb}/ncbi.tre ./ncbi.tre

    ${params.meganpath}/daa-meganizer \\
        --in ${daa_file} \\
        --mapDB ${params.megandb}/megan-map.db \\
        --threads ${task.cpus} \\
        --topPercent 0.5 \\
        --minSupportPercent 0 \\
        --lcaAlgorithm longReads \\
        --longReads true \\
        --verbose >> Meganize_BlastX_${assembler}.log 2>&1

    # DAA summary counts (ID + counts)
    paste <(${params.meganpath}/daa2info -P ${params.meganpath}/.MEGAN.def -i ${daa_file} -c2c Taxonomy | awk '{print \$1}') \\
          <(${params.meganpath}/daa2info -P ${params.meganpath}/.MEGAN.def -i ${daa_file} -p -c2c Taxonomy | awk -F'\\t' '{print \$1, \$2}' OFS='\\t') \\
          > ${sample_id}_${assembler}_contigs_blastx_daa_summary_count.tsv

    diamond view \\
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle qlen slen \\
        --daa ${daa_file} \\
        > ${sample_id}_${assembler}_contigs_blastx_diamondview.tsv

    # MEGAN taxonomy TSV (Krona input)
    ${params.meganpath}/daa2info \\
        -i ${daa_file} -c2c Taxonomy \\
        -o ${sample_id}_${assembler}_contigs_blastx.tsv

    # Krona 
    ktImportTaxonomy \\
        -tax /export/database/krona/taxonomy \\
        -t 1 -m 2 \\
        ${sample_id}_${assembler}_contigs_blastx.tsv \\
        -o ${sample_id}_${assembler}_contigs_blastx_krona.html
"""
}

process Parse_BlastX_Contigs {
    tag { "${assembler} | ${sample_id}" }
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/blast/", mode: 'copy'
    label 'lowmem'
    conda "${baseDir}/env/vs.yml"

    input:
    tuple val(sample_id), val(assembler), val(diamondview_file)

    output:
    tuple val(sample_id), val(assembler), file("*"), emit: parse_blastx_ch


    when:
    script:
    """
    # TMP FIX
    #
    #source /opt/conda/etc/profile.d/conda.sh
    #conda activate all_in_one_pipeline
    /opt/conda/envs/all_in_one_pipeline/bin/python ${params.scripts}/VS_MD_diamond_parser_linFilt_pandas_v4.py \\
      -i ${diamondview_file} \\
      -t blastx \\
      -r allRanks \\
      -v ${params.vhunter} \\
      -n ${params.ncbi_taxa} \\
      >> "Parse_BLASTX_${assembler}.log" 2>&1
    """
}
