#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {
    Remove_Common_Flora_rRNA_reads;
    RiboDetector_Remove_rRNA_reads;
    Post_removal_rRNA_reads_fastqc;
    Multiqc_QC_rRNA_removal
} from './modules/local/remove_rRNA_from_reads/main.nf'


workflow Remove_rRNA_Reads_Workflow {

    take:
    input_ch  // tuple(sample_id, sr_file, lr_file, cpus, mem)

    main:

    def rrna_method = (params.rRNA_method ?: 'bbmap_minimap2').toString().toLowerCase()

    input_ch.view { "Input to Remove_rRNA_Reads_Workflow using method '${rrna_method}': $it" }

    rRNA_host_remove_short_ch    = Channel.empty()
    rRNA_host_remove_short_r1_ch = Channel.empty()
    rRNA_host_remove_short_r2_ch = Channel.empty()
    rRNA_host_remove_long_ch     = Channel.empty()
    rRNA_reads_short_ch          = Channel.empty()
    rRNA_reads_long_ch           = Channel.empty()
    rRNA_all_files_ch            = Channel.empty()


    if (rrna_method == 'bbmap_minimap2') {

        Remove_Common_Flora_rRNA_reads(input_ch)

        rRNA_host_remove_short_ch    = Remove_Common_Flora_rRNA_reads.out.rRNA_host_remove_short
        rRNA_host_remove_short_r1_ch = Remove_Common_Flora_rRNA_reads.out.rRNA_host_remove_short_r1
        rRNA_host_remove_short_r2_ch = Remove_Common_Flora_rRNA_reads.out.rRNA_host_remove_short_r2
        rRNA_host_remove_long_ch     = Remove_Common_Flora_rRNA_reads.out.rRNA_host_remove_long

        rRNA_reads_short_ch          = Remove_Common_Flora_rRNA_reads.out.rRNA_reads_short
        rRNA_reads_long_ch           = Remove_Common_Flora_rRNA_reads.out.rRNA_reads_long

        rRNA_all_files_ch            = Remove_Common_Flora_rRNA_reads.out.rRNA_all_files_ch
    }

    else if (rrna_method == 'ribodetector') {

        RiboDetector_Remove_rRNA_reads(input_ch)

        rRNA_host_remove_short_ch    = RiboDetector_Remove_rRNA_reads.out.rRNA_host_remove_short
        rRNA_host_remove_short_r1_ch = RiboDetector_Remove_rRNA_reads.out.rRNA_host_remove_short_r1
        rRNA_host_remove_short_r2_ch = RiboDetector_Remove_rRNA_reads.out.rRNA_host_remove_short_r2
        rRNA_host_remove_long_ch     = RiboDetector_Remove_rRNA_reads.out.rRNA_host_remove_long

        rRNA_reads_short_ch          = RiboDetector_Remove_rRNA_reads.out.rRNA_reads_short
        rRNA_reads_long_ch           = RiboDetector_Remove_rRNA_reads.out.rRNA_reads_long

        rRNA_all_files_ch            = RiboDetector_Remove_rRNA_reads.out.rRNA_all_files_ch
    }

    else {
        error "Invalid --rRNA_method '${params.rRNA_method}'. Valid options are: bbmap_minimap2, ribodetector"
    }


    def sr_rRNA_removed = rRNA_host_remove_short_ch.map { sid, f -> tuple(sid, f) }
    def lr_rRNA_removed = rRNA_host_remove_long_ch.map  { sid, f -> tuple(sid, f) }
    def sr_rRNA         = rRNA_reads_short_ch.map       { sid, f -> tuple(sid, f) }
    def lr_rRNA         = rRNA_reads_long_ch.map        { sid, f -> tuple(sid, f) }

    def all_fastqs_by_sample = sr_rRNA_removed
                                .mix(lr_rRNA_removed)
                                .mix(sr_rRNA)
                                .mix(lr_rRNA)
                                .groupTuple()
                                .map { sid, files -> tuple(sid, files) }

    all_fastqs_by_sample.view { "Input to Post_removal_rRNA_reads_fastqc: $it" }

    Post_removal_rRNA_reads_fastqc(all_fastqs_by_sample)

    Multiqc_QC_rRNA_removal(
        Post_removal_rRNA_reads_fastqc.out.post_rRNA_removal_fastqc_ch.collect()
    )


    emit:
    rRNA_host_remove_short    = rRNA_host_remove_short_ch
    rRNA_host_remove_short_r1 = rRNA_host_remove_short_r1_ch
    rRNA_host_remove_short_r2 = rRNA_host_remove_short_r2_ch
    rRNA_host_remove_long     = rRNA_host_remove_long_ch

    rRNA_reads_short          = rRNA_reads_short_ch
    rRNA_reads_long           = rRNA_reads_long_ch

    rRNA_all_files_ch         = rRNA_all_files_ch
}