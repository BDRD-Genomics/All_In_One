nextflow.enable.dsl = 2

/*
Reads Taxonomic Classifier Workflow

Input:
  classifier_reads_ch:
    tuple(sample_id, reads)

  short reads:
    tuple("${sample_id}_short", [r1, r2])

  long reads:
    tuple("${sample_id}_long", [lr])

Classifiers:
  - Kraken2
  - Bracken (runs from the Kraken2 report)
  - Sourmash
  - MetaPhlAn4
  - GOTTCHA2

Taxpasta:
  - Kraken2: supported
  - MetaPhlAn4: supported
  - Sourmash: native sketch/gather/kreport only; not sent to Taxpasta
  - GOTTCHA2: native output only for now
  */


/*
================================================================================
Classifier processes
================================================================================
*/

process KRAKEN2_CLASSIFY {

    tag "${sample_id}"

    label 'kraken2'

    publishDir "${params.outdir}/${params.project_id}/reads_taxonomic_classifier/kraken2", mode: 'copy'
    errorStrategy 'ignore'
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id),
          path("${sample_id}.kraken2.report.txt"),
          path("${sample_id}.kraken2.output.txt"),
          emit: kraken2_results

    script:
    def read_args

    if (reads.size() == 2) {
        read_args = "--paired ${reads[0]} ${reads[1]}"
    } else if (reads.size() == 1) {
        read_args = "${reads[0]}"
    } else {
        read_args = reads.join(" ")
    }

    """

    kraken2 \\
        --db ${params.kraken2_db} \\
        --threads ${task.cpus} \\
        --report ${sample_id}.kraken2.report.txt \\
        --output ${sample_id}.kraken2.output.txt \\
        ${params.kraken2_extra_args} \\
        ${read_args}
    """
}


process BRACKEN_ABUNDANCE {

    tag "${sample_id}"

    label 'kraken2'

    publishDir "${params.outdir}/${params.project_id}/reads_taxonomic_classifier/bracken", mode: 'copy'
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), path(kraken_report)

    output:
    tuple val(sample_id),
          path("${sample_id}.bracken.${params.bracken_level}.tsv"),
          path("${sample_id}.bracken.${params.bracken_level}.report.txt"),
          emit: bracken_results

    script:

    def read_length = sample_id.endsWith('_long') ?
        params.bracken_read_length_long :
        params.bracken_read_length_short

    def extra_args = params.bracken_extra_args ?: ""

    """
    set -euo pipefail

    bracken \
        -d ${params.kraken2_db} \
        -i ${kraken_report} \
        -o ${sample_id}.bracken.${params.bracken_level}.tsv \
        -w ${sample_id}.bracken.${params.bracken_level}.report.txt \
        -r ${read_length} \
        -l ${params.bracken_level} \
        -t ${params.bracken_threshold} \
        ${extra_args}
    """
}

process KRAKEN2_KRONA {

    tag "${sample_id}"

    label 'krona'

    publishDir "${params.outdir}/${params.project_id}/reads_taxonomic_classifier/kraken2", mode: 'copy'
    errorStrategy 'ignore'
    input:
    tuple val(sample_id), path(krakenReport)

    output:
    tuple val(sample_id),
          path("${sample_id}_krona.html"),
          emit: krona_results

    script:
    """
    ktImportTaxonomy -tax ${params.krona_db} -m 3 -t 5 ${krakenReport} -o ${sample_id}_krona.html
    """
}

process KRAKEN2_KRONA_SNAPSHOT {

    tag "${sample_id}"

    label 'krona'

    publishDir "${params.outdir}/${params.project_id}/reads_taxonomic_classifier/kraken2/snapshots", mode: 'copy'

    container "${params.krona_screenshot_container}"
    errorStrategy 'ignore'
    input:
    tuple val(sample_id), path(krona_html)

    output:
    tuple val(sample_id),
          path("${sample_id}.krona.snapshot.svg"),
          path("${sample_id}.krona.snapshot.png"),
          emit: kraken2_krona_snapshot

    script:
    """
    set -euo pipefail

    python ${params.scripts}/screenshot_krona_snapshot_button.py \\
      -i ${krona_html} \\
      -o ${sample_id}.krona.snapshot.svg \\
      -p ${sample_id}.krona.snapshot.png \\
      --width ${params.krona_snapshot_width} \\
      --height ${params.krona_snapshot_height} \\
      --wait-ms ${params.krona_snapshot_wait_ms} \\
      --popup-timeout-ms ${params.krona_snapshot_popup_timeout_ms} \\
      --download-timeout-ms ${params.krona_snapshot_download_timeout_ms} \\
      --png-width ${params.krona_snapshot_png_width} \
      --font-size ${params.krona_font}

    ls -lh \\
      ${sample_id}.krona.snapshot.svg \\
      ${sample_id}.krona.snapshot.png
    """
}

process SOURMASH_SKETCH {

    tag "${sample_id}"

    label 'sourmash'

    publishDir "${params.outdir}/${params.project_id}/reads_taxonomic_classifier/sourmash/sketches", mode: 'copy'
    errorStrategy 'ignore'
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id),
          path("${sample_id}.sig"),
          emit: sourmash_sketches

    script:
    def read_args = reads.join(" ")

    """

    sourmash sketch dna \\
        -p ${params.sourmash_scaled_param} \\
        --merge ${sample_id} \\
        -o ${sample_id}.sig \\
        ${read_args}
    """
}


process SOURMASH_GATHER {

    tag "${sample_id}"

    label 'sourmash'

    publishDir "${params.outdir}/${params.project_id}/reads_taxonomic_classifier/sourmash/gather", mode: 'copy'
    errorStrategy 'ignore'
    input:
    tuple val(sample_id), path(sig)

    output:
    tuple val(sample_id),
          path("${sample_id}.sourmash.gather.csv"),
          emit: sourmash_gather_results

    script:
    def extra_args = params.sourmash_extra_args ?: ""

    """

    sourmash gather \\
        ${sig} \\
        ${params.sourmash_db} \\
        -o ${sample_id}.sourmash.gather.csv \\
        ${extra_args}
    """
}


process SOURMASH_TAX_KREPORT {

    tag "${sample_id}"

    label 'sourmash'

    publishDir "${params.outdir}/${params.project_id}/reads_taxonomic_classifier/sourmash/kreport", mode: 'copy'
    errorStrategy 'ignore'
    input:
    tuple val(sample_id), path(gather_csv)

    output:
    tuple val(sample_id),
          path("${sample_id}.sourmash.kreport.txt"),
          emit: sourmash_kreport_results

    script:
    def extra_args = (
        params.sourmash_tax_extra_args == null ||
        params.sourmash_tax_extra_args == false ||
        params.sourmash_tax_extra_args == true
    ) ? "" : params.sourmash_tax_extra_args.toString()

    """

    sourmash tax metagenome \\
        --gather-csv ${gather_csv} \\
        --taxonomy ${params.sourmash_taxonomy} \\
        --output-format kreport \\
        --output-base ${sample_id}.sourmash \\
        ${extra_args}

    if [ ! -f "${sample_id}.sourmash.kreport.txt" ]; then
        echo "ERROR: Expected Sourmash kreport was not created: ${sample_id}.sourmash.kreport.txt" >&2
        echo "Files in work directory:" >&2
        ls -lah >&2
        exit 1
    fi
    """
}

process MASH {
    tag "${sample_id}"
    label 'mash'
    publishDir "${params.outdir}/${params.project_id}/reads_taxonomic_classifier/mash", mode: 'copy'
    errorStrategy 'ignore'
    input:
    tuple val(sample_id), path(trimmed_reads)

    output:
    tuple val(sample_id),
          path("${sample_id}.mash.txt"),
          emit: mash_results

    script:
    """
    mash screen -w -p ${task.cpus} ${params.mash_db} ${trimmed_reads} | sort -gr > ${sample_id}.mash.txt
    """
}


process METAPHLAN4_CLASSIFY {

    tag "${sample_id}"

    label 'metaphlan4'

    publishDir "${params.outdir}/${params.project_id}/reads_taxonomic_classifier/metaphlan4", mode: 'copy'
    errorStrategy 'ignore'
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id),
          path("${sample_id}.metaphlan4.profile.txt"),
          path("${sample_id}.metaphlan4.bowtie2.bz2"),
          emit: metaphlan4_results

    script:
    def input_reads = reads.join(",")
    def db_arg = params.metaphlan4_db_arg ?: "--bowtie2db"

    """

    metaphlan \\
        --offline \\
        ${input_reads} \\
        --input_type fastq \\
        --nproc ${task.cpus} \\
        ${db_arg} ${params.metaphlan4_db} \\
        --mapout ${sample_id}.metaphlan4.bowtie2.bz2 \\
        ${params.metaphlan4_extra_args} \\
        -o ${sample_id}.metaphlan4.profile.txt
    """
}


process GOTTCHA_CLASSIFY {

    tag "${sample_id}"

    label 'gottcha'

    publishDir "${params.outdir}/${params.project_id}/reads_taxonomic_classifier/gottcha", mode: 'copy'
    errorStrategy 'ignore'
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id),
          path("${sample_id}.gottcha.tsv"),
          emit: gottcha_results

    script:
    def read_args = reads.collect { "-i ${it}" }.join(" ")

    """
    set -euo pipefail

    if command -v gottcha2.py >/dev/null 2>&1; then
        GOTTCHA_CMD="gottcha2.py"
    elif command -v gottcha2 >/dev/null 2>&1; then
        GOTTCHA_CMD="gottcha2"
    else
        echo "ERROR: Could not find gottcha2.py or gottcha2 in PATH" >&2
        exit 1
    fi

    \${GOTTCHA_CMD} \\
        ${read_args} \\
        -d ${params.gottcha_db} \\
        -t ${task.cpus} \\
        -o ${sample_id}.gottcha \\
        ${params.gottcha_extra_args}

    if ls ${sample_id}.gottcha*.tsv >/dev/null 2>&1; then
        cp \$(ls ${sample_id}.gottcha*.tsv | head -n 1) ${sample_id}.gottcha.tsv
    else
        echo "ERROR: Could not find GOTTCHA TSV output for ${sample_id}" >&2
        echo "Files in work directory:" >&2
        ls -lah >&2
        exit 1
    fi
    """
}


/*
================================================================================
Taxpasta processes

Only tools with native/officially supported Taxpasta-compatible profiler output
are sent here:
  - Kraken2
  - MetaPhlAn4

Sourmash is not sent to Taxpasta. It can still produce a Sourmash kreport
when params.run_sourmash_kreport = true.
================================================================================
*/

process TAXPASTA_KRAKEN2 {

    tag "taxpasta_kraken2"

    label 'taxpasta'

    publishDir "${params.outdir}/${params.project_id}/reads_taxonomic_classifier/taxpasta", mode: 'copy'
    errorStrategy 'ignore'
    input:
    path kraken2_reports

    output:
    path "kraken2_taxpasta.tsv", emit: kraken2_taxpasta

    script:
    """
    set -euo pipefail

    n_reports=\$(ls *.kraken2.report.txt 2>/dev/null | wc -l)

    echo "Found \${n_reports} Kraken2 report(s) for Taxpasta"

    if [ "\${n_reports}" -eq 1 ]; then
        echo "One Kraken2 report found; using taxpasta standardise"

        taxpasta standardise \\
            --profiler kraken2 \\
            --output kraken2_taxpasta.tsv \\
            \$(ls *.kraken2.report.txt)

    elif [ "\${n_reports}" -gt 1 ]; then
        echo "Multiple Kraken2 reports found; using taxpasta merge"

        taxpasta merge \\
            --profiler kraken2 \\
            --output kraken2_taxpasta.tsv \\
            \$(ls *.kraken2.report.txt)

    else
        echo "ERROR: No Kraken2 reports found for Taxpasta" >&2
        exit 1
    fi
    """
}


process TAXPASTA_METAPHLAN4 {

    tag "taxpasta_metaphlan4"

    label 'taxpasta'

    publishDir "${params.outdir}/${params.project_id}/reads_taxonomic_classifier/taxpasta", mode: 'copy'
    errorStrategy 'ignore'
    input:
    path metaphlan_profiles

    output:
    path "metaphlan4_taxpasta.tsv", emit: metaphlan4_taxpasta

    script:
    """
    set -euo pipefail

    n_profiles=\$(ls *.metaphlan4.profile.txt 2>/dev/null | wc -l)

    echo "Found \${n_profiles} MetaPhlAn4 profile(s) for Taxpasta"

    if [ "\${n_profiles}" -eq 1 ]; then
        echo "One MetaPhlAn4 profile found; using taxpasta standardise"

        taxpasta standardise \\
            --profiler metaphlan \\
            --output metaphlan4_taxpasta.tsv \\
            \$(ls *.metaphlan4.profile.txt)

    elif [ "\${n_profiles}" -gt 1 ]; then
        echo "Multiple MetaPhlAn4 profiles found; using taxpasta merge"

        taxpasta merge \\
            --profiler metaphlan \\
            --output metaphlan4_taxpasta.tsv \\
            \$(ls *.metaphlan4.profile.txt)

    else
        echo "ERROR: No MetaPhlAn4 profiles found for Taxpasta" >&2
        exit 1
    fi
    """
}


/*
Summary
*/

process READS_TAXONOMIC_CLASSIFIER_SUMMARY {

    tag "reads_taxonomic_classifier_summary"

    label 'summary'

    publishDir "${params.outdir}/${params.project_id}/reads_taxonomic_classifier/summary", mode: 'copy'
    errorStrategy 'ignore'
    input:
    path kraken_reports
    path kraken_outputs
    path bracken_abundance_tables
    path bracken_reports
    path metaphlan_profiles
    path metaphlan_bowtie2
    path sourmash_sigs
    path sourmash_gathers
    path sourmash_kreports
    path gottcha_reports
    path taxpasta_tables

    output:
    path "reads_taxonomic_classifier_summary.tsv", emit: classifier_summary

    script:
    """
    set -euo pipefail

    echo -e "tool\\tfile_type\\tsample_or_table\\tfile" > reads_taxonomic_classifier_summary.tsv

    for f in ${kraken_reports.join(' ')}; do
        [ -f "\$f" ] && echo -e "kraken2\\treport\\t\$(basename "\$f" .kraken2.report.txt)\\t\$f" >> reads_taxonomic_classifier_summary.tsv
    done

    for f in ${kraken_outputs.join(' ')}; do
        [ -f "\$f" ] && echo -e "kraken2\\tper_read_output\\t\$(basename "\$f" .kraken2.output.txt)\\t\$f" >> reads_taxonomic_classifier_summary.tsv
    done

    for f in ${bracken_abundance_tables.join(' ')}; do
        [ -f "\$f" ] && echo -e "bracken\\tabundance_table\\t\$(basename "\$f" .bracken.${params.bracken_level}.tsv)\\t\$f" >> reads_taxonomic_classifier_summary.tsv
    done

    for f in ${bracken_reports.join(' ')}; do
        [ -f "\$f" ] && echo -e "bracken\\treestimated_report\\t\$(basename "\$f" .bracken.${params.bracken_level}.report.txt)\\t\$f" >> reads_taxonomic_classifier_summary.tsv
    done

    for f in ${metaphlan_profiles.join(' ')}; do
        [ -f "\$f" ] && echo -e "metaphlan4\\tprofile\\t\$(basename "\$f" .metaphlan4.profile.txt)\\t\$f" >> reads_taxonomic_classifier_summary.tsv
    done

    for f in ${metaphlan_bowtie2.join(' ')}; do
        [ -f "\$f" ] && echo -e "metaphlan4\\tbowtie2\\t\$(basename "\$f" .metaphlan4.bowtie2.bz2)\\t\$f" >> reads_taxonomic_classifier_summary.tsv
    done

    for f in ${sourmash_sigs.join(' ')}; do
        [ -f "\$f" ] && echo -e "sourmash\\tsketch\\t\$(basename "\$f" .sig)\\t\$f" >> reads_taxonomic_classifier_summary.tsv
    done

    for f in ${sourmash_gathers.join(' ')}; do
        [ -f "\$f" ] && echo -e "sourmash\\tgather\\t\$(basename "\$f" .sourmash.gather.csv)\\t\$f" >> reads_taxonomic_classifier_summary.tsv
    done

    for f in ${sourmash_kreports.join(' ')}; do
        [ -f "\$f" ] && echo -e "sourmash\\tkreport\\t\$(basename "\$f" .sourmash.kreport.txt)\\t\$f" >> reads_taxonomic_classifier_summary.tsv
    done

    for f in ${gottcha_reports.join(' ')}; do
        [ -f "\$f" ] && echo -e "gottcha2\\treport\\t\$(basename "\$f" .gottcha.tsv)\\t\$f" >> reads_taxonomic_classifier_summary.tsv
    done

    for f in ${taxpasta_tables.join(' ')}; do
        [ -f "\$f" ] && echo -e "taxpasta\\ttable\\t\$(basename "\$f")\\t\$f" >> reads_taxonomic_classifier_summary.tsv
    done
    """
}

process Exercise_Report_Taxonomy {
    publishDir "${params.outdir}/${params.project_id}/exercise_report", mode: 'copy'
    label 'lowmem'
    conda "$baseDir/env/pandas_env.yml"
    errorStrategy 'ignore'

    input:
    val(sample_id)

    output:
    path "*.docx", emit: report_out_ch

    script:
    """
    python ${params.scripts}/exercise_report_Taxonomy_only_v2.py -s ${sample_id} -i ${params.outdir}/${params.project_id}/ -o . -p ${params.project_id} -t ${baseDir}/templates/Taxonomy_Isolate_Template_v2.docx
    """

}



/*
Workflow
*/

workflow Reads_Taxonomic_Classifier_Workflow {

    take:
    classifier_reads_ch

    main:

    /*
    Run classifiers independently
    */

    if (params.run_kraken2) {
        KRAKEN2_CLASSIFY(classifier_reads_ch)

        if (params.run_bracken) {
            def bracken_input_ch = KRAKEN2_CLASSIFY.out.kraken2_results
                .map { sample_id, report, output -> tuple(sample_id, report) }

            BRACKEN_ABUNDANCE(bracken_input_ch)
        }

        def krona_input_ch = KRAKEN2_CLASSIFY.out.kraken2_results.map { sample_id, report, output -> tuple(sample_id, report) }
        KRAKEN2_KRONA(krona_input_ch)
        KRAKEN2_KRONA_SNAPSHOT(KRAKEN2_KRONA.out.krona_results)
    }

    if (params.run_sourmash) {
        SOURMASH_SKETCH(classifier_reads_ch)
        SOURMASH_GATHER(SOURMASH_SKETCH.out.sourmash_sketches)

        if (params.run_sourmash_kreport) {
            SOURMASH_TAX_KREPORT(SOURMASH_GATHER.out.sourmash_gather_results)
        }
    }

    if (params.run_mash) {
        MASH(classifier_reads_ch)
    }

    if (params.run_metaphlan4) {
        METAPHLAN4_CLASSIFY(classifier_reads_ch)
    }

    if (params.run_gottcha) {
        GOTTCHA_CLASSIFY(classifier_reads_ch)
    }


    /*
    Taxpasta

    Only native Taxpasta-compatible profiler outputs:
      - Kraken2
      - MetaPhlAn4

    Sourmash and GOTTCHA2 are not sent to Taxpasta.
    */

    taxpasta_tables_ch = Channel.empty()

    if (params.run_taxpasta && params.run_kraken2) {

        kraken2_reports_for_taxpasta_ch = KRAKEN2_CLASSIFY.out.kraken2_results
            .map { sample_id, report, output -> report }
            .collect()

        TAXPASTA_KRAKEN2(kraken2_reports_for_taxpasta_ch)

        taxpasta_tables_ch = taxpasta_tables_ch.mix(
            TAXPASTA_KRAKEN2.out.kraken2_taxpasta
        )
    }

    if (params.run_taxpasta && params.run_metaphlan4) {

        metaphlan_profiles_for_taxpasta_ch = METAPHLAN4_CLASSIFY.out.metaphlan4_results
            .map { sample_id, profile, bowtie2out -> profile }
            .collect()

        TAXPASTA_METAPHLAN4(metaphlan_profiles_for_taxpasta_ch)

        taxpasta_tables_ch = taxpasta_tables_ch.mix(
            TAXPASTA_METAPHLAN4.out.metaphlan4_taxpasta
        )
    }


    /*
    Summary
    */
    
    if (params.make_taxonomic_classifier_summary) {

        kraken_reports_summary_ch = params.run_kraken2 ?
            KRAKEN2_CLASSIFY.out.kraken2_results.map { sample_id, report, output -> report }.collect() :
            Channel.value([])

        kraken_outputs_summary_ch = params.run_kraken2 ?
            KRAKEN2_CLASSIFY.out.kraken2_results.map { sample_id, report, output -> output }.collect() :
            Channel.value([])

        bracken_abundance_summary_ch =
            (params.run_kraken2 && params.run_bracken) ?
                BRACKEN_ABUNDANCE.out.bracken_results.map { sample_id, abundance, report -> abundance }.collect() :
                Channel.value([])

        bracken_reports_summary_ch =
            (params.run_kraken2 && params.run_bracken) ?
                BRACKEN_ABUNDANCE.out.bracken_results.map { sample_id, abundance, report -> report }.collect() :
                Channel.value([])

        metaphlan_profiles_summary_ch = params.run_metaphlan4 ?
            METAPHLAN4_CLASSIFY.out.metaphlan4_results.map { sample_id, profile, bowtie2out -> profile }.collect() :
            Channel.value([])

        metaphlan_bowtie2_summary_ch = params.run_metaphlan4 ?
            METAPHLAN4_CLASSIFY.out.metaphlan4_results.map { sample_id, profile, bowtie2out -> bowtie2out }.collect() :
            Channel.value([])

        sourmash_sigs_summary_ch = params.run_sourmash ?
            SOURMASH_SKETCH.out.sourmash_sketches.map { sample_id, sig -> sig }.collect() :
            Channel.value([])

        sourmash_gathers_summary_ch = params.run_sourmash ?
            SOURMASH_GATHER.out.sourmash_gather_results.map { sample_id, gather -> gather }.collect() :
            Channel.value([])

        sourmash_kreports_summary_ch =
            (params.run_sourmash && params.run_sourmash_kreport) ?
                SOURMASH_TAX_KREPORT.out.sourmash_kreport_results.map { sample_id, kreport -> kreport }.collect() :
                Channel.value([])

        gottcha_reports_summary_ch = params.run_gottcha ?
            GOTTCHA_CLASSIFY.out.gottcha_results.map { sample_id, gottcha -> gottcha }.collect() :
            Channel.value([])

        taxpasta_tables_summary_ch = params.run_taxpasta ?
            taxpasta_tables_ch.collect() :
            Channel.value([])

        READS_TAXONOMIC_CLASSIFIER_SUMMARY(
            kraken_reports_summary_ch,
            kraken_outputs_summary_ch,
            bracken_abundance_summary_ch,
            bracken_reports_summary_ch,
            metaphlan_profiles_summary_ch,
            metaphlan_bowtie2_summary_ch,
            sourmash_sigs_summary_ch,
            sourmash_gathers_summary_ch,
            sourmash_kreports_summary_ch,
            gottcha_reports_summary_ch,
            taxpasta_tables_summary_ch
        )	    
	

        def exercise_report_ch = MASH.out.mash_results
            .join(KRAKEN2_KRONA_SNAPSHOT.out.kraken2_krona_snapshot)
            .map{ sid, mashout, svg, png -> sid }
            .map{ sid -> sid.replaceAll('_long|_short',"")}
	Exercise_Report_Taxonomy( exercise_report_ch )

    }




    /*
    Emits
    */

    emit:
    kraken2_results = params.run_kraken2 ? KRAKEN2_CLASSIFY.out.kraken2_results : Channel.empty()
    krona_results = params.run_kraken2 ? KRAKEN2_KRONA.out.krona_results : Channel.empty()

    bracken_results =
        (params.run_kraken2 && params.run_bracken) ?
            BRACKEN_ABUNDANCE.out.bracken_results :
            Channel.empty()

    sourmash_sketches = params.run_sourmash ? SOURMASH_SKETCH.out.sourmash_sketches : Channel.empty()

    sourmash_gather_results = params.run_sourmash ? SOURMASH_GATHER.out.sourmash_gather_results : Channel.empty()

    sourmash_kreport_results =
        (params.run_sourmash && params.run_sourmash_kreport) ?
            SOURMASH_TAX_KREPORT.out.sourmash_kreport_results :
            Channel.empty()

    mash_results = params.run_mash ? MASH.out.mash_results : Channel.empty()

    metaphlan4_results = params.run_metaphlan4 ? METAPHLAN4_CLASSIFY.out.metaphlan4_results : Channel.empty()

    gottcha_results = params.run_gottcha ? GOTTCHA_CLASSIFY.out.gottcha_results : Channel.empty()

    taxpasta_tables = params.run_taxpasta ? taxpasta_tables_ch : Channel.empty()
    taxonomy_exercise_report = Exercise_Report_Taxonomy.out.report_out_ch
}
