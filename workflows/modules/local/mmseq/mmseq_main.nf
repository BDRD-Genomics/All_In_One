#!/usr/bin/env nextflow



nextflow.enable.dsl=2

/*
MMSEQS process set used by MMSEQ_EasySearch_Workflow
Defines:
  - MMSEQS_CREATEDB_QUERY  : build per-query DB
  - MMSEQS_SEARCH          : search queryDB vs target DB
  - MMSEQS_CONVERTALIS     : convert alignment DB to TSV
*/

// Helpers for publishDir paths 
def cleanSid(id) {
  def s = id.toString()
  s.contains('__') ? s.split('__', 2)[1] : s
}
def compId(id) {
  id.toString()
}


process MMSEQS_CREATEDB_QUERY {
  tag {sample_id}
  errorStrategy 'ignore'
  publishDir { "${params.outdir}/${params.run_id}/${ cleanSid(sample_id) }/mmseqs/${ compId(sample_id) }" }, mode: 'copy'
  label 'lowmem'

  input:
  tuple val(sample_id), path(query_fa)

  output:
  tuple val(sample_id), path('queryDB_dir/')

  script:
  """
  set -euo pipefail
  mkdir -p queryDB_dir
  mmseqs createdb ${query_fa} queryDB_dir/queryDB ${params.mmseqs_createdb_opts ?: ''} 2>&1 | tee createdb_${sample_id}.log
  """
}

process MMSEQS_SEARCH {
  tag { sample_id }
  errorStrategy 'ignore'

  publishDir {
    "${params.outdir}/${params.run_id}/${cleanSid(sample_id)}/mmseqs/${compId(sample_id)}"
  }, mode: 'copy'

  label 'optimized_MMSEQ_SEARCH'

  input:
  tuple val(sample_id), path(qdb, stageAs: 'input_qdb')
  path target_db

  output:
  tuple val(sample_id), path('alnDB_dir/'), path("${sample_id}.log")

  script:
  """
  set -euo pipefail

  mkdir -p alnDB_dir_pre
  mkdir -p alnDB_dir

  if [[ -d input_qdb/queryDB_dir ]]; then
      cp -rL input_qdb/queryDB_dir real_qdb
  else
      cp -rL input_qdb real_qdb
  fi

  test -f real_qdb/queryDB.dbtype || {
      echo "ERROR: queryDB.dbtype not found"
      find real_qdb -maxdepth 3 -type f -print
      exit 1
  }

  mmseqs search \
      real_qdb/queryDB \
      ${params.mmseqs_nt_db} \
      alnDB_dir_pre/alnDB_pre \
      tmp_${sample_id} \
      --threads ${task.cpus} \
      --local-tmp ${params.tmp_dir} \
      ${params.mmseqs_search_opts ?: ''} \
      2>&1 | tee ${sample_id}.log

  mmseqs filterdb \
      alnDB_dir_pre/alnDB_pre \
      alnDB_dir/alnDB \
      --extract-lines ${params.mmseqs_search_max_seqs} \
      --filter-column 12 \
      --sort-entries 2
  """
}

process MMSEQS_CONVERTALIS {
  tag { sample_id }
  errorStrategy 'ignore'

  publishDir {
    "${params.outdir}/${params.run_id}/${cleanSid(sample_id)}/mmseqs/${compId(sample_id)}"
  }, mode: 'copy'

  label 'optimized_MMSEQ_CONTIG_SEARCH'

  input:
  tuple val(sample_id), path(queryDB_dir), path(alnDB_dir)
  path target_db

  output:
  tuple val(sample_id), path("${sample_id}.mmseq.out")

  script:
  def fmt = params.mmseqs_format ?: 'query,target,evalue,bits,alnlen,pident,qstart,qend,tstart,tend'

  """
  set -euo pipefail

  if [[ -d ${queryDB_dir}/queryDB_dir ]]; then
      cp -rL ${queryDB_dir}/queryDB_dir real_qdb
  else
      cp -rL ${queryDB_dir} real_qdb
  fi

  if [[ -d ${alnDB_dir}/alnDB_dir ]]; then
      cp -rL ${alnDB_dir}/alnDB_dir real_adb
  else
      cp -rL ${alnDB_dir} real_adb
  fi

  test -f real_qdb/queryDB.dbtype || {
      echo "ERROR: queryDB.dbtype not found"
      find real_qdb -maxdepth 3 -type f -print
      exit 1
  }

  test -f real_adb/alnDB.dbtype || {
      echo "ERROR: alnDB.dbtype not found"
      find real_adb -maxdepth 3 -type f -print
      exit 1
  }

  mmseqs convertalis \
      real_qdb/queryDB \
      ${params.mmseqs_nt_db} \
      real_adb/alnDB \
      ${sample_id}.mmseq.out \
      --threads ${task.cpus} \
      --format-mode 4 \
      --format-output "${fmt}" \
      ${params.mmseqs_convertalis_opts ?: ''}
  """
}

process Parse_MMSEQ_Contigs {
  tag { "${assembler} | ${sample_id}" }
  publishDir path: { "${params.outdir}/${params.run_id}/${sample_id}/blast/" }, mode: 'copy'
  label 'lowmem'
  errorStrategy 'ignore'

  input:
  tuple val(sample_id), val(assembler), path(mmseq_out)

  output:
  tuple val(sample_id), val(assembler), file("*.parsed"), emit: mmseqs_parsed_ch

  script:
  """
  set -euo pipefail
  # TMP FIX
  # remove header
  tail -n +2 "${mmseq_out}" > tmp && mv tmp "${mmseq_out}"

  python ${params.scripts}/VS_MD_diamond_parser_linFilt_Mar2026_fast.py \\
    -i "${mmseq_out}" \\
    -t mmseqs \\
    -v ${params.megan_mdb} \\
    -n ${params.ncbi_taxa} 
  """
}
