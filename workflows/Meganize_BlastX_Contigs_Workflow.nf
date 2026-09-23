#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {
    Meganize_BlastX_Contigs
} from './modules/local/blastx/main.nf'

workflow Meganize_BlastX_Contigs_Workflow {
    take:
    blastx_input_ch

    main:
    blastx_input_ch.view { " Meganizing: $it" }
    Meganize_BlastX_Contigs(blastx_input_ch)

    emit:
    meganized_ch = Meganize_BlastX_Contigs.out.meganize_blastx_contigs_ch
}
