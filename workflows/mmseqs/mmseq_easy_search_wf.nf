#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
MMSEQS Easy Search workflow (standalone)
- Takes fasta queries, builds per-query DBs, searches a target DB, converts alignments to TSV
- Input:
    queries_ch    -> (sample_id, path(query.fa))
    target_db_ch  -> path(target_db_prefix)  [single value channel]
- Output:
    out           -> (sample_id, path(results.tsv))
*/

// IMPORTANT: alias the processes so they don't clash with other includes
include {
  MMSEQS_CREATEDB_QUERY as ES_CREATEDB
  MMSEQS_SEARCH         as ES_SEARCH
  MMSEQS_CONVERTALIS    as ES_CONVERT
} from '../modules/local/mmseq/mmseq_main.nf'

workflow MMSEQ_EasySearch_Workflow {
  take:
    queries_ch
    target_db_ch

  main:
    // keep single-value channel uses simple; avoid into/dup if you like
    def db_for_search  = target_db_ch.map { it }
    def db_for_convert = target_db_ch.map { it }

    // 1) Build per-sample query DBs -> (sample_id, path(queryDB))
    def built_q = ES_CREATEDB(queries_ch)

    // 2) Search -> (sample_id, path(alnDB.dbtype), log)
    def built_q_dir = built_q.map { sid, qdb -> tuple(sid, qdb.parent) }
    def hits        = ES_SEARCH(built_q_dir, db_for_search)

    // 3) Pair queryDB_dir with alnDB_dir -> (sample_id, queryDB_dir, alnDB_dir)
    def q_and_hits = built_q_dir.join(hits).map { sid, qdir, adb_dbtype, _log ->
      tuple(sid, qdir, adb_dbtype.parent)
    }

    // 4) Convert to TSV -> (sample_id, path(tsv))
    def tsvs = ES_CONVERT(q_and_hits, db_for_convert)

  emit:
    out = tsvs
}