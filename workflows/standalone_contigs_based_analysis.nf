#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { CONTIG_based_analysis } from './contigs_based_analysis.nf'

workflow {
    Channel.fromPath(params.samplesheet).splitCsv(header: true).map{ row ->
        //tuple(row.sample_id, row.contigs, row.assembler)
        def fq1 = row.fastq_1?.contains('No_Read') ? null : file(row.fastq_1, checkIfExists: false)
        def fq2 = row.fastq_2?.contains('No_Read') ? null : file(row.fastq_2, checkIfExists: false)
        def lr  = row.long_read?.contains('No_Read') ? null : file(row.long_read, checkIfExists: false)
        tuple(row.sample_id, fq1, fq2, lr, row.contigs, row.mode, row.assembler)
    }.set {contigs_based_analysis_ch}

    CONTIG_based_analysis(contigs_based_analysis_ch)
}
