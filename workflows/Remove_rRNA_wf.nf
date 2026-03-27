#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {
    Remove_Common_Flora_rRNA_reads;
    Post_removal_rRNA_reads_fastqc;
    Multiqc_QC_rRNA_removal
} from './modules/local/remove_rRNA_from_reads/main.nf'

workflow Remove_rRNA_Reads_Workflow {
    take:
    input_ch  // tuple(sample_id, sr_file, lr_file, cpus, mem)

    main:
    input_ch.view { " Input to Remove_Common_Flora_rRNA_reads: $it" }
    Remove_Common_Flora_rRNA_reads(input_ch)
    //  Gather all FASTQs by sample for FastQC 
    def sr_rRNA_removed = Remove_Common_Flora_rRNA_reads.out.rRNA_host_remove_short.map { sid, f -> tuple(sid, f) }
    def lr_rRNA_removed = Remove_Common_Flora_rRNA_reads.out.rRNA_host_remove_long.map  { sid, f -> tuple(sid, f) }
    def lr_rRNA         = Remove_Common_Flora_rRNA_reads.out.rRNA_reads_long.map        { sid, f -> tuple(sid, f) }
    def sr_rRNA         = Remove_Common_Flora_rRNA_reads.out.rRNA_reads_short.map       { sid, f -> tuple(sid, f) }

    // Union of all produced FASTQs (some may be absent per sample)
    def all_fastqs_by_sample = sr_rRNA_removed
                                .mix(lr_rRNA_removed)
                                .mix(lr_rRNA)
                                .mix(sr_rRNA)
                                .groupTuple()                      // -> (sample_id, [files...])
                                .map { sid, files -> tuple(sid, files) }


    Post_removal_rRNA_reads_fastqc( all_fastqs_by_sample )

    Multiqc_QC_rRNA_removal(Post_removal_rRNA_reads_fastqc.out.post_rRNA_removal_fastqc_ch.collect())


    emit:
    rRNA_host_remove_short = Remove_Common_Flora_rRNA_reads.out.rRNA_host_remove_short
    rRNA_host_remove_short_r1 = Remove_Common_Flora_rRNA_reads.out.rRNA_host_remove_short_r1
    rRNA_host_remove_short_r2 = Remove_Common_Flora_rRNA_reads.out.rRNA_host_remove_short_r2
    rRNA_host_remove_long  = Remove_Common_Flora_rRNA_reads.out.rRNA_host_remove_long
    rRNA_all_files_ch      = Remove_Common_Flora_rRNA_reads.out.rRNA_all_files_ch
}
