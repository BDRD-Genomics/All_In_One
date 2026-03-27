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
def cleanSid = { id -> 
  def s = id.toString()
  s.contains('__') ? s.split('__', 2)[1] : s
}
def compId = { id -> id.toString() }


process MMSEQS_CREATEDB_QUERY {
  tag {sample_id}
  errorStrategy 'ignore'
  //publishDir "${params.outdir}/${params.project_id}/${sample_id}/mmseq/", mode: 'copy'
  publishDir "${params.outdir}/${params.project_id}/${ cleanSid(sample_id) }/mmseqs/${ compId(sample_id) }", mode: 'copy'
  label 'lowmem'
  conda "$baseDir/env/mmseqs2.yml"

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
  tag {sample_id}
  errorStrategy 'ignore'
  //maxRetries 1
  //publishDir "${params.outdir}/${params.project_id}/${sample_id}/mmseq/", mode: 'copy'
  publishDir "${params.outdir}/${params.project_id}/${ cleanSid(sample_id) }/mmseqs/${ compId(sample_id) }", mode: 'copy'
  label 'optimized_MMSEQ_SEARCH'
  conda "$baseDir/env/mmseqs2.yml"

  input:
  tuple val(sample_id), path(qdb, stageAs: 'input_qdb')
  path target_db

  output:
  //tuple val(sample_id), val('alnDB'), path('alnDB*'), path("search_${sample_id}.log")
  //tuple val(sample_id), path('alnDB')      // emit a REAL path object
  //tuple val(sample_id), path('alnDB*'), path("search_${sample_id}.log")
  //tuple val(sample_id), path('alnDB.dbtype'), path("${sample_id}.log")
  //tuple val(sample_id), path('alnDB_dir/alnDB.dbtype'), path("${sample_id}.log")
  tuple val(sample_id), path('alnDB_dir/'), path("${sample_id}.log")

  script:
  """
  set -euo pipefail
  mkdir -p alnDB_dir
  #[[ -f input_qdb/queryDB.dbtype ]] \
  #  || { echo "ERROR: input_qdb/queryDB.dbtype not staged"; ls -la input_qdb/ || true; exit 1; }
  cp -rL input_qdb real_qdb
  # diagnostic
  ls -la real_qdb/
  mmseqs search \
        real_qdb/queryDB \
	${params.mmseqs_nt_db} \
	alnDB_dir/alnDB \
	tmp_${sample_id} \
	--threads ${task.cpus} \
        --local-tmp /database/tmp \
	${params.mmseqs_search_opts ?: ''} 2>&1 | tee ${sample_id}.log

  #[[ -f alnDB_dir/alnDB.dbtype ]] \
  #  || { echo "ERROR: alnDB.dbtype missing after search"; exit 1; }
  """
}

process MMSEQS_CONVERTALIS {
  tag {sample_id}
  errorStrategy 'ignore'
  //maxRetries 1
  //publishDir "${params.outdir}/${params.project_id}/${sample_id}/mmseq/", mode: 'copy'
  publishDir "${params.outdir}/${params.project_id}/${ cleanSid(sample_id) }/mmseqs/${ compId(sample_id) }", mode: 'copy'
  label 'optimized_MMSEQ_contigs_against_NT_metaspades'
  conda "$baseDir/env/mmseqs2.yml"

  input:
  tuple val(sample_id), path(queryDB_dir), path(alnDB_dir)
  path target_db

  output:
  tuple val(sample_id), path("${sample_id}.mmseq.out")

  script:
  def fmt = params.mmseqs_format ?: 'query,target,evalue,bits,alnlen,pident,qstart,qend,tstart,tend'
  """
  set -euo pipefail
  cp -rL ${queryDB_dir} real_qdb
  cp -rL ${alnDB_dir}/alnDB_dir   real_adb
  # diagnostic
  ls -la real_adb/
  mmseqs convertalis \
    real_qdb/queryDB \
    ${params.mmseqs_nt_db} \
    real_adb/alnDB \
    ${sample_id}.mmseq.out \
    --threads ${task.cpus} \
    --format-mode 4 \
    --format-output \
    "${fmt}" \
    ${params.mmseqs_convertalis_opts ?: ''} 
  """
}

process Parse_MMSEQ_Contigs {
  tag { "${assembler} | ${sample_id}" }
  publishDir path: { "${params.outdir}/${params.project_id}/${sample_id}/blast/" }, mode: 'copy'
  label 'lowmem'
  errorStrategy 'ignore'
  conda "${baseDir}/env/vs.yml"

  input:
  tuple val(sample_id), val(assembler), path(mmseq_out)

  output:
  tuple val(sample_id), val(assembler), file("*.parsed"), emit: mmseqs_parsed_ch
  tuple val(sample_id), val(assembler), file("*.log"), emit: mmseqs_log_ch

  script:
  """
  set -euo pipefail
  # TMP FIX
  #source /opt/conda/etc/profile.d/conda.sh
  #conda activate all_in_one_pipeline
  # remove header
  tail -n +2 "${mmseq_out}" > tmp && mv tmp "${mmseq_out}"

  /opt/conda/envs/all_in_one_pipeline/bin/python ${params.scripts}/VS_MD_diamond_parser_linFilt_pandas_v4.py \\
    -i "${mmseq_out}" \\
    -t mmseqs \\
    -r allRanks \\
    -v ${params.vhunter} \\
    -n ${params.ncbi_taxa} \\
    >> "Parse_MMSEQ_${assembler}.log" 2>&1
  """
}

