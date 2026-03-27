#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {
    BlastX_contigs
} from './modules/local/blastx/main.nf'

workflow BlastX_Contigs_Workflow {

    take:
    blast_input_ch

    main:
    blast_input_ch.view { " Starting BlastX on: $it" }
    BlastX_contigs(blast_input_ch)

    // Ensure the output tuple is: (sample_id, assembler, file)
    def blastx_daa_ch = BlastX_contigs.out.blastx_daa_ch

    def dragonflye_blastx_contigs_ch = blastx_daa_ch.filter { it[1] == 'dragonflye' }
    def metaspades_blastx_contigs_ch = blastx_daa_ch.filter { it[1] == 'metaspades' }
    def spades_blastx_contigs_ch     = blastx_daa_ch.filter { it[1] == 'spades' }
    def unicycler_blastx_contigs_ch  = blastx_daa_ch.filter { it[1] == 'unicycler' }

    emit:
        dragonflye_blastx_contigs_ch = dragonflye_blastx_contigs_ch
        metaspades_blastx_contigs_ch = metaspades_blastx_contigs_ch
        spades_blastx_contigs_ch     = spades_blastx_contigs_ch
        unicycler_blastx_contigs_ch  = unicycler_blastx_contigs_ch
}
