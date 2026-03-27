#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
========================================================================================
   MMSEQ - run EasySearch per assembler
========================================================================================
*/

/*
Inputs:
  contigs_in_ch : (asm, sid, path(fasta))
  target_db_ch  : single-value channel of mmseqs DB prefix
Output:
  out           : (asm, sid, path(results.tsv))
*/

// Alias the underlying MMseqs processes HERE so they are declared ONCE
include {
  MMSEQS_CREATEDB_QUERY as ES_CREATEDB
  MMSEQS_SEARCH         as ES_SEARCH
  MMSEQS_CONVERTALIS    as ES_CONVERT
} from '../modules/local/mmseq/mmseq_main.nf'

workflow ContigSearch_AllAssemblers {
  take:
    contigs_in_ch
    target_db_ch

  main:
    def SEP = '__'

    // Build composite key to keep (asm, sid) distinct
    def queries_ch = contigs_in_ch.map { asm, sid, fa ->
      tuple("${asm}${SEP}${sid}".toString(), fa)
    }

    // Duplicate single-value DB channel safely
    def db_for_search  = target_db_ch.map { it }
    def db_for_convert = target_db_ch.map { it }

    // EasySearch inline (no extra include)
    def built_q    = ES_CREATEDB(queries_ch)                                 // (comp, qdb)
    def built_qdir = built_q.map { comp, qdb -> tuple(comp, qdb.parent) }    // (comp, qdir)
    def hits       = ES_SEARCH(built_qdir, db_for_search)                    // (comp, alnDB.dbtype, log)
    def paired     = built_qdir.join(hits).map { comp, qdir, adbtype, _log ->
                       tuple(comp, qdir, adbtype.parent) }                   // (comp, qdir, aln_dir)
    def tsvs       = ES_CONVERT(paired, db_for_convert)                      // (comp, tsv)

    // Split composite back to (asm, sid)
    def out_ch = tsvs.map { comp, tsv ->
      def parts = comp.toString().split(SEP, 2)
      def asm = parts[0]
      def sid = parts.size() > 1 ? parts[1] : ''
      tuple(asm, sid, tsv)
    }

  emit:
    out = out_ch
}