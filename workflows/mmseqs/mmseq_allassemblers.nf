#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
MMSEQS: run EasySearch for ALL assemblers in one shot.
Inputs:
  contigs_all_ch : (asm, sid, path(fasta))    // asm ∈ {dragonflye, metaspades, unicycler, ...}
  target_db_ch   : single-value channel with mmseqs DB prefix
Outputs:
  emit dragon    : (sid, path(tsv))
  emit meta      : (sid, path(tsv))
  emit uni       : (sid, path(tsv))
  emit all       : (asm, sid, path(tsv))
*/

// Declare/alias the processes ONCE here:
include {
  MMSEQS_CREATEDB_QUERY as ES_CREATEDB
  MMSEQS_SEARCH         as ES_SEARCH
  MMSEQS_CONVERTALIS    as ES_CONVERT
} from '../modules/local/mmseq/mmseq_main.nf'

workflow MMSEQ_AllAssemblers {
  take:
    contigs_all_ch    // (asm, sid, fasta)
    target_db_ch

  main:
    def SEP = '__'

    // Build composite key to keep (asm, sid) unique across samples & assemblers
    def queries_ch = contigs_all_ch.map { asm, sid, fa ->
      tuple("${asm}${SEP}${sid}".toString(), fa)
    }

    // Duplicate the single-value DB channel for fan-out safety
    def db_for_search  = target_db_ch.map { it }
    def db_for_convert = target_db_ch.map { it }

    // Inline "EasySearch"
    def built_q     = ES_CREATEDB(queries_ch)                                  // (comp, qdb)
    def built_qdir  = built_q.map { comp, qdb -> tuple(comp, qdb.parent) }     // (comp, qdir)
    def hits        = ES_SEARCH(built_qdir, db_for_search)                     // (comp, alnDB.dbtype, log)
    def paired      = built_qdir.join(hits).map { comp, qdir, adbtype, _log ->
                       tuple(comp, qdir, adbtype.parent) }                      // (comp, qdir, aln_dir)
    def tsvs        = ES_CONVERT(paired, db_for_convert)                        // (comp, tsv)

    // Split composite back to (asm, sid) and branch to per-assembler streams
    def all_out = tsvs.map { comp, tsv ->
      def parts = comp.toString().split(SEP, 2)
      def asm = parts[0]
      def sid = parts.size() > 1 ? parts[1] : ''
      tuple(asm, sid, tsv)                                                     // (asm, sid, tsv)
    }

    def split = all_out.branch {
      dragon : it[0] == 'dragonflye'
      medaka : it[0] == 'dragonflye_medaka'
      meta   : it[0] == 'metaspades'
      uni    : it[0] == 'unicycler'
      raven    : it[0] == 'raven'
      clc : it[0] == 'clc'
    }

  emit:
    dragon = split.dragon.map { asm, sid, tsv -> tuple(sid, tsv) }             // (sid, tsv)
    medaka = split.medaka.map { asm, sid, tsv -> tuple(sid, tsv) }             // (sid, tsv)
    meta   = split.meta  .map { asm, sid, tsv -> tuple(sid, tsv) }             // (sid, tsv)
    uni    = split.uni   .map { asm, sid, tsv -> tuple(sid, tsv) }             // (sid, tsv)
    raven  = split.raven   .map { asm, sid, tsv -> tuple(sid, tsv) }             // (sid, tsv)
    clc    = split.clc   .map { asm, sid, tsv -> tuple(sid, tsv) }              // (sid, tsv)
    all    = all_out                                                           // (asm, sid, tsv)
}
