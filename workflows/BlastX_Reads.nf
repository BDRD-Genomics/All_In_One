#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {
    BlastX_reads
} from './modules/local/blastx/main.nf'

workflow BlastX_Reads_Workflow {
    take:
    short_ch  // (sample_id, fastq_file, 'short')
    long_ch   // (sample_id, fastq_file, 'long')

    main:

    short_ch.view { " Short read input: $it" }
    long_ch.view  { " Long read input: $it" }

    blastx_reads_input_ch = short_ch.mix(long_ch)

    blastx_reads_input_ch.view { " Final BlastX input (one per read type): $it" }

    BlastX_reads(blastx_reads_input_ch)

    emit:
    blastx_reads_ch = BlastX_reads.out.blastx_reads_ch
    short_reads_blastx_channel = BlastX_reads.out.short_reads_blastx_channel
    long_reads_blastx_channel = BlastX_reads.out.long_reads_blastx_channel
}
