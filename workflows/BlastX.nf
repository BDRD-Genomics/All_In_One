#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {
    BlastX_contigs
} from './modules/local/blastx/main.nf'

workflow BlastX_Workflow {
    take:
    blast_input_ch

    main:
    blast_input_ch
        .view { "Starting BlastX on: $it" } \
        | BlastX_contigs

    emit:
    blastx_results_ch = BlastX_contigs.out.blastx_contigs_ch
}

