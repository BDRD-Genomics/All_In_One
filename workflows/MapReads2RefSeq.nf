#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { Map_Reads_2_RefSeq } from './modules/local/map2refseq/main.nf'
//include { Post_trim_fastqc   } from './modules/local/map2refseq/main.nf' 
include { Post_host_remove_fastqc   } from './modules/local/map2refseq/main.nf' 
include { Multiqc_QC_host_removal  } from './modules/local/map2refseq/main.nf' 
workflow MapReads_2_RefSeq {
    take:
    trimmed_short_ch
    trimmed_long_ch

    main:
    // Short reads -> resource-scaled tuple
    def short_ch = trimmed_short_ch.map { sample_id, fq1, fq2 ->
        def total_size = [fq1, fq2].findAll().sum { it.size() }
        def size_gb = total_size / 1e9
        def cpus = size_gb < 1 ? 32 : size_gb < 5 ? 32 : size_gb < 10 ? 64 : 128
        def mem  = size_gb < 1 ? '32GB' : size_gb < 5 ? '64GB' : size_gb < 10 ? '64GB' : '128GB'
        tuple(sample_id, fq1, fq2, null, 'short', cpus, mem)
    }

    // Long reads -> resource-scaled tuple
    def long_ch = trimmed_long_ch.map { sample_id, lr ->
        def size_gb = lr.size() / 1e9
        def cpus = size_gb < 1 ? 2 : size_gb < 5 ? 4 : size_gb < 10 ? 8 : 16
        def mem  = size_gb < 1 ? '8GB' : size_gb < 5 ? '16GB' : size_gb < 10 ? '32GB' : '64GB'
        tuple(sample_id, null, null, lr, 'long', cpus, mem)
    }

    // Merge both types & map
    def merged_ch = short_ch.mix(long_ch)

    // Run mapping
    Map_Reads_2_RefSeq( merged_ch )

    //  Gather all FASTQs by sample for FastQC 
    def sr_removed = Map_Reads_2_RefSeq.out.host_removed_fq_ch_sr.map { sid, f -> tuple(sid, f) }
    def lr_removed = Map_Reads_2_RefSeq.out.host_removed_fq_ch_lr.map { sid, f -> tuple(sid, f) }
    def lr_host    = Map_Reads_2_RefSeq.out.host_fq_ch_lr.map         { sid, f -> tuple(sid, f) }
    def sr_host    = Map_Reads_2_RefSeq.out.host_fq_ch_sr.map         { sid, f -> tuple(sid, f) }

    // Union of all produced FASTQs (some may be absent per sample)
    def all_fastqs_by_sample = sr_removed
                                .mix(lr_removed)
                                .mix(lr_host)
                                .mix(sr_host)
                                .groupTuple()                      // -> (sample_id, [files...])
                                .map { sid, files -> tuple(sid, files) }

    // Run FastQC once per sample on whatever files exist
    Post_host_remove_fastqc( all_fastqs_by_sample )

    Multiqc_QC_host_removal(Post_host_remove_fastqc.out.post_host_remove_fastqc_ch.collect())

    emit:
    host_removed_fq_ch_sr = Map_Reads_2_RefSeq.out.host_removed_fq_ch_sr
    host_removed_fq_ch_lr = Map_Reads_2_RefSeq.out.host_removed_fq_ch_lr
    posttrim_fastqc_ch    = Post_host_remove_fastqc.out.post_host_remove_fastqc_ch
    host_removed_fq_ch_r1 = Map_Reads_2_RefSeq.out.host_removed_fq_ch_r1
    host_removed_fq_ch_r2 = Map_Reads_2_RefSeq.out.host_removed_fq_ch_r2
}
