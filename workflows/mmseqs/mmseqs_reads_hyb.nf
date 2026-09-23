#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
MMSEQS: run EasySearch for read queries in a separate workflow.

Inputs:
  reads_all_ch  : (read_type, sid, path(fasta))
  target_db_ch  : single-value channel with mmseqs DB prefix

Outputs:
  short : (sid, path(tsv))
  long  : (sid, path(tsv))
  all   : (read_type, sid, path(tsv))
*/

include {
  MMSEQS_CREATEDB_QUERY as ES_CREATEDB
  MMSEQS_SEARCH         as ES_SEARCH
  MMSEQS_CONVERTALIS    as ES_CONVERT
} from '../modules/local/mmseq/mmseq_main.nf'

workflow MMSEQ_AllReads {

  take:
    reads_all_ch
    target_db_ch

  main:
    def SEP = '__'

    // normalize for the consolidated mmseq module
    def queries_ch = reads_all_ch.map { read_type, sid, fa ->
      def comp = "${read_type}${SEP}${sid}".toString()
      def qdb = "queryDB_${read_type}".toString()
      def adb = "alnDB_${read_type}".toString()
      tuple(comp, fa, qdb, adb, 'reads.mmseq.out')
    }

    def createdb_in = queries_ch.map { comp, fa, query_prefix, aln_prefix, out_suffix ->
      tuple(comp, fa, query_prefix)
    }

    def search_meta = queries_ch.map { comp, fa, query_prefix, aln_prefix, out_suffix ->
      tuple(comp, aln_prefix, out_suffix)
    }

    def built_q = ES_CREATEDB(createdb_in)

    def search_in = built_q.join(search_meta).map { comp, query_dir, query_prefix, aln_prefix, out_suffix ->
      tuple(comp, query_dir, query_prefix, aln_prefix)
    }

    def hits = ES_SEARCH(search_in, target_db_ch)

    def paired = built_q
      .join(hits)
      .join(queries_ch.map { comp, fa, query_prefix, aln_prefix, out_suffix ->
        tuple(comp, out_suffix)
      })
      .map { comp, query_dir, query_prefix, aln_dir, _qp, aln_prefix, _log, out_suffix ->
        tuple(comp, query_dir, aln_dir, query_prefix, aln_prefix, out_suffix)
      }

    def tsvs = ES_CONVERT(paired, target_db_ch)

    def all_out = tsvs.map { comp, tsv ->
      def parts = comp.toString().split(SEP, 2)
      def read_type = parts[0]
      def sid = parts.size() > 1 ? parts[1] : ''
      tuple(read_type, sid, tsv)
    }

    def short_out = all_out
      .filter { read_type, sid, tsv -> read_type == 'short' }
      .map    { read_type, sid, tsv -> tuple(sid, tsv) }

    def long_out = all_out
      .filter { read_type, sid, tsv -> read_type == 'long' }
      .map    { read_type, sid, tsv -> tuple(sid, tsv) }

  emit:
    mmseqs_short_reads_ch = short_out
    mmseqs_long_reads_ch  = long_out
    mseqs_all_reads_ch   = all_out
}
