#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
Read Mapping Auto Workflow (Reference Pull + CoverM)

Input:
  readmapping_auto_input_ch: tuple(sample_id, kraken_report, bracken_table, mash_screen, blast_tsvs, long_read)

Pulls a local reference for every kraken2/bracken/mash/BLAST(contigs) hit, optionally dereplicates near-identical
references by ANI (params.readmapping_auto_cluster_dedup), then runs CoverM against the
final reference set.
*/

include { GET_REFERENCES; DEDUP_ANI; COVERM } from './modules/local/readmapping_auto/main.nf'

workflow ReadMapping_Auto_Workflow {
    take:
    readmapping_auto_input_ch

    main:

    def get_refs_input_ch = readmapping_auto_input_ch.map { sid, kraken_report, bracken_table, mash_screen, blast_tsvs, long_read ->
        tuple(sid, kraken_report, bracken_table, mash_screen, blast_tsvs)
    }
    GET_REFERENCES(get_refs_input_ch)

    def reference_fasta_ch = params.readmapping_auto_cluster_dedup ?
        DEDUP_ANI(GET_REFERENCES.out.references_ch).dedup_fasta_ch :
        GET_REFERENCES.out.references_ch.map { sid, fasta, tsv, genomes -> tuple(sid, fasta) }

    def long_read_ch = readmapping_auto_input_ch.map { sid, kraken_report, bracken_table, mash_screen, blast_tsvs, long_read ->
        tuple(sid, long_read)
    }

    COVERM(reference_fasta_ch.join(long_read_ch))

    emit:
    references_ch  = GET_REFERENCES.out.references_ch
    not_found_ch   = GET_REFERENCES.out.not_found_ch
    bam_ch         = COVERM.out.bam_ch
    stats_ch       = COVERM.out.stats_ch
}
