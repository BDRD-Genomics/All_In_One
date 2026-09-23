#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {
  Meganize_ShortReads_BlastX;
  Meganize_LongReads_BlastX;
} from '../modules/local/blastx/main.nf'

workflow Megan_ShortReads_WF {
  take:
    in_ch                    // (sid, path(daa)) — just one channel

  main:
    Meganize_ShortReads_BlastX(in_ch)      // no extra inputs, no joins, no grouping

  emit:
    meganize_short_reads_ch        = Meganize_ShortReads_BlastX.out.meganize_short_reads_ch
    all_shortreads_meganized_files = Meganize_ShortReads_BlastX.out.all_shortreads_meganized_files
    meganized_short_daa_ch         = Meganize_ShortReads_BlastX.out.meganized_short_daa_ch
}

workflow Megan_LongReads_WF {
  take:
    in_ch                    // (sid, path(daa)) — just one channel

  main:
    Meganize_LongReads_BlastX(in_ch)      // no extra inputs, no joins, no grouping

  emit:
    meganize_long_reads_ch        = Meganize_LongReads_BlastX.out.meganize_long_reads_ch
    all_longreads_meganized_files = Meganize_LongReads_BlastX.out.all_longreads_meganized_files
    meganized_long_daa_ch         = Meganize_LongReads_BlastX.out.long_reads_daa_meganized_ch
}