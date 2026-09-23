#!/usr/bin/env bash
# =============================================================================
# CLC Server Client Library
# =============================================================================
# Shell functions wrapping CLC Genomics Server Command Line Tools.
# Source this file in any script that needs to interact with CLC Server.
#
#   https://resources.qiagenbioinformatics.com/manuals/clcservercommandlinetools/current/User_Manual.pdf
#
# Architecture:
#   - Insert size is specified at import time (ngs_import_illumina); assembly
#     and mapping auto-detect paired distances from the read metadata.
#
# Required environment variables:
#   CLC_SERVER_HOST  - Server hostname or IP address
#   CLC_SERVER_PORT  - Server port (default: 7777)
#   CLC_SERVER_USER  - Username
#   CLC_SERVER_PASS  - Password or keystore token
#   CLC_DATA_ROOT    - Root server path for this pipeline's data (clc:// URL)
#
# Optional:
#   CLC_GRID         - Grid preset name for -G flag (e.g. "slurm").
#                      Must be set for all analysis commands to work.
#                      Leave empty only if the server has direct execution enabled.
#   CLC_JOB_TIMEOUT  - Seconds to wait for a grid job (default: 3600)
# =============================================================================

set -euo pipefail

# Defaults
CLC_SERVER_CLI="${CLC_SERVER_CLI:-clcserver}"
CLC_LOG_FILE="${CLC_LOG_FILE:-clc_server.log}"
CLC_JOB_TIMEOUT="${CLC_JOB_TIMEOUT:-3600}"
CLC_GRID="${CLC_GRID:-}"

# =============================================================================
# Internal helper functions
# =============================================================================

_clc_log() {
    local level="$1"; shift
    local msg="$*"
    local ts; ts=$(date '+%Y-%m-%d %H:%M:%S')
    echo "[${ts}] [${level}] ${msg}" | tee -a "${CLC_LOG_FILE}" >&2
}

_clc_check_env() {
    local missing=()
    [[ -z "${CLC_SERVER_HOST:-}" ]] && missing+=("CLC_SERVER_HOST")
    [[ -z "${CLC_SERVER_USER:-}" ]] && missing+=("CLC_SERVER_USER")
    [[ -z "${CLC_SERVER_PASS:-}" ]] && missing+=("CLC_SERVER_PASS")
    [[ -z "${CLC_DATA_ROOT:-}"   ]] && missing+=("CLC_DATA_ROOT")
    if [[ ${#missing[@]} -gt 0 ]]; then
        _clc_log "ERROR" "Missing required environment variables: ${missing[*]}"
        return 1
    fi
}

_clc_conn_args() {
    # Build base connection args array; add -G when CLC_GRID is set.
    # Used by sync commands (export, mapping stats). Async commands use _clc_submit_async.
    local -n _arr="$1"
    _arr=(
        "${CLC_SERVER_CLI}"
        -S "${CLC_SERVER_HOST}"
        -P "${CLC_SERVER_PORT:-7777}"
        -U "${CLC_SERVER_USER}"
        -W "${CLC_SERVER_PASS}"
    )
    [[ -n "${CLC_GRID:-}" ]] && _arr+=( -G "${CLC_GRID}" )
}

_clc_submit_async() {
    local action="$1"; shift
    _clc_check_env || return 1
    local cmd_arr=(
        "${CLC_SERVER_CLI}"
        -S "${CLC_SERVER_HOST}"
        -P "${CLC_SERVER_PORT:-7777}"
        -U "${CLC_SERVER_USER}"
        -W "${CLC_SERVER_PASS}"
    )
    [[ -n "${CLC_GRID:-}" ]] && cmd_arr+=( -G "${CLC_GRID}" )
    cmd_arr+=( -A "${action}" -Y "$@" )

    _clc_log "DEBUG" "Running ${cmd_arr[*]}"
    local raw
    raw=$("${cmd_arr[@]}" </dev/null 2>&1)
    _clc_log "DEBUG" "async raw output: ${raw}"

    # CLC Server returns the job ID on a line like: "Message: Process ID: sy0lwxtk1dqt3n"
    local job_id
    job_id=$(echo "${raw}" | grep -i 'Process ID:' | grep -oP '(?<=Process ID:\s)\S+' | tail -1)
    if [[ -z "${job_id}" ]]; then
        _clc_log "ERROR" "No job ID (Process ID) found in output — action '${action}' may have failed or is unavailable."
        _clc_log "ERROR" "Raw output: ${raw}"
        return 1
    fi
    echo "${job_id}"
}

_clc_retrieve() {
    # Get the output CLC URL(s) of a completed job using -R <job_id>.
    # Always exits 0 — callers guard with || retrieve_out="" so set -e never triggers.
    local job_id="$1"
    local cmd_arr=(
        "${CLC_SERVER_CLI}"
        -S "${CLC_SERVER_HOST}"
        -P "${CLC_SERVER_PORT:-7777}"
        -U "${CLC_SERVER_USER}"
        -W "${CLC_SERVER_PASS}"
        -R "${job_id}"
    )
    _clc_log "DEBUG" "Running ${cmd_arr[*]}"
    "${cmd_arr[@]}" </dev/null 2>&1 || true
}

_clc_result_url() {
    # Extract the Nth CLC URL from _clc_retrieve output (1-indexed, default 1).
    # Safe to call even when retrieve output is empty — returns empty string.
    # Usage: url=$(_clc_result_url "${retrieve_output}" [n])
    local raw="$1"
    local n="${2:-1}"
    echo "${raw}" | grep -oP 'clc://\S+' | sed -n "${n}p" || true
}

_clc_ls_data_url() {
    # List a CLC server folder and return the clcUrl of the first data object
    # matching type_filter. Skips Log objects when no filter is given.
    #
    # -A ls output per object (field name case varies; // separates objects):
    #   Type: CLC file (Sequence List)
    #   ClcUrl: clc://host:8443/<opaque-id>          ← no spaces, use this
    #   ClcUrl Simple: clc://host:8443/CLC_Data/...  ← may have spaces, skip
    #
    # Usage: url=$(_clc_ls_data_url <folder_clc_url> [type_substring])
    local folder_url="$1"
    local type_filter="${2:-}"

    local ls_out
    ls_out=$(
        "${CLC_SERVER_CLI}" \
            -S "${CLC_SERVER_HOST}" \
            -P "${CLC_SERVER_PORT:-7777}" \
            -U "${CLC_SERVER_USER}" \
            -W "${CLC_SERVER_PASS}" \
            -A ls -t "${folder_url}" 2>&1
    ) || true
    _clc_log "DEBUG" "_clc_ls_data_url[${type_filter:-any}] output: ${ls_out}"

    printf '%s\n' "${ls_out}" | python3 -c "
import sys
tf = sys.argv[1].lower() if len(sys.argv) > 1 and sys.argv[1] else ''
current_type = ''
for raw in sys.stdin:
    ll = raw.lower().rstrip()
    if ll.startswith('type:'):
        current_type = ll
    elif ll.startswith('clcurl:') and 'simple' not in ll:
        match = (tf in current_type) if tf else ('(log)' not in current_type)
        if match:
            # split on first colon only to preserve clc:// in the URL
            url = raw.split(':', 1)[1].strip()
            print(url)
            sys.exit(0)
" "${type_filter}" || true
}

clc_wait_job() {
    # Poll a CLC job until completed, failed, or timed out.
    # Usage: clc_wait_job <job_id> [timeout_seconds]
    local job_id="$1"
    local timeout="${2:-${CLC_JOB_TIMEOUT}}"

    if [[ -z "${job_id}" ]]; then
        _clc_log "ERROR" "clc_wait_job: empty job_id"
        return 1
    fi

    local start_time; start_time=$(date +%s)
    _clc_log "INFO" "Waiting for job ${job_id} (timeout=${timeout}s)"

    while true; do
        local elapsed=$(( $(date +%s) - start_time ))
        if [[ ${elapsed} -ge ${timeout} ]]; then
            _clc_log "ERROR" "Job ${job_id} timed out after ${timeout}s"
            return 1
        fi

        local cmd_arr=(
            "${CLC_SERVER_CLI}"
            -S "${CLC_SERVER_HOST}"
            -P "${CLC_SERVER_PORT:-7777}"
            -U "${CLC_SERVER_USER}"
            -W "${CLC_SERVER_PASS}"
            -I "${job_id}"
        )
        local status
        status=$("${cmd_arr[@]}" </dev/null 2>&1) || true

        if echo "${status}" | grep -qiE "is done: Yes"; then
            _clc_log "INFO" "Job ${job_id} completed (${elapsed}s)"
            return 0
        elif echo "${status}" | grep -qiE "has errors: Yes"; then
            _clc_log "ERROR" "Job ${job_id} failed: ${status}"
            _clc_log "ERROR" "Full status output: ${status}"
            return 1
        fi

        _clc_log "INFO" "Job ${job_id} running (${elapsed}s elapsed) — waiting 60s"
        sleep 60
    done
}

# =============================================================================
# Server filesystem operations (no -G — data ops run directly on the server)
# =============================================================================

clc_mkdir() {
    # Create a folder on the CLC Server.
    # Always creates a new folder; CLC auto-increments the name if one already exists.
    # Use clc_ensure_dir when you want "create only if missing" behaviour.
    # Usage: clc_mkdir <parent_clc_url> <folder_name>
    local parent_url="$1"
    local folder_name="$2"
    _clc_check_env || return 1
    local cmd_arr=(
        "${CLC_SERVER_CLI}"
        -S "${CLC_SERVER_HOST}"
        -P "${CLC_SERVER_PORT:-7777}"
        -U "${CLC_SERVER_USER}"
        -W "${CLC_SERVER_PASS}"
        -A mkdir
        -t "${parent_url}"
        -n "${folder_name}"
    )
    _clc_log "DEBUG" "Running ${cmd_arr[*]}"
    "${cmd_arr[@]}" </dev/null || true
    _clc_log "INFO" "Created folder: ${parent_url}/${folder_name}"
}

clc_ensure_dir() {
    # Create a folder only if it does not already exist.
    # Uses ls to check first; if ls succeeds the folder exists and we skip mkdir.
    # Usage: clc_ensure_dir <parent_clc_url> <folder_name>
    local parent_url="$1"
    local folder_name="$2"
    local full_url="${parent_url}/${folder_name}"
    _clc_check_env || return 1
    if "${CLC_SERVER_CLI}" \
            -S "${CLC_SERVER_HOST}" \
            -P "${CLC_SERVER_PORT:-7777}" \
            -U "${CLC_SERVER_USER}" \
            -W "${CLC_SERVER_PASS}" \
            -A ls -t "${full_url}" >/dev/null 2>&1; then
        _clc_log "DEBUG" "Folder already exists, skipping mkdir: ${full_url}"
    else
        clc_mkdir "${parent_url}" "${folder_name}"
    fi
}

clc_mkdir_p() {
    # Create a nested path of folders under CLC_DATA_ROOT (mkdir -p equivalent).
    # Uses clc_ensure_dir at each level so existing folders are never duplicated.
    # Usage: clc_mkdir_p "SampleA/qc/trimmed"
    local rel_path="$1"
    local current="${CLC_DATA_ROOT}"
    IFS='/' read -ra parts <<< "${rel_path}"
    for part in "${parts[@]}"; do
        [[ -z "${part}" ]] && continue
        clc_ensure_dir "${current}" "${part}"
        current="${current}/${part}"
    done
}

clc_rm() {
    # Delete a server object or folder (recursive).
    # Usage: clc_rm <clc_url>
    local server_url="$1"
    _clc_check_env || return 1
    _clc_log "INFO" "Deleting: ${server_url}"
    local cmd_arr=(
        "${CLC_SERVER_CLI}"
        -S "${CLC_SERVER_HOST}"
        -P "${CLC_SERVER_PORT:-7777}"
        -U "${CLC_SERVER_USER}"
        -W "${CLC_SERVER_PASS}"
        -A rm
        -t "${server_url}"
        -r
    )
    _clc_log "DEBUG" "Running ${cmd_arr[*]}"
    "${cmd_arr[@]}" </dev/null || true
}

clc_ls() {
    # List contents of a server folder.
    # Usage: clc_ls <clc_url>
    local server_url="${1:-${CLC_DATA_ROOT}}"
    _clc_check_env || return 1
    local cmd_arr=(
        "${CLC_SERVER_CLI}"
        -S "${CLC_SERVER_HOST}"
        -P "${CLC_SERVER_PORT:-7777}"
        -U "${CLC_SERVER_USER}"
        -W "${CLC_SERVER_PASS}"
        -A ls
        -t "${server_url}"
    )
    _clc_log "DEBUG" "Running ${cmd_arr[*]}"
    "${cmd_arr[@]}" </dev/null
}

# =============================================================================
# Import (upload local files to CLC Server)
# All import actions require CLC_GRID to be set.
# =============================================================================


clc_import_paired_reads() {
    # Import FASTQ reads to CLC Server using ngs_import_illumina.
    # Insert sizes are stored in the read metadata and used by downstream analysis.
    # Usage: clc_import_paired_reads <r1_local> <r2_local> <dest_clc_url>
    #                                [insert_lb=250] [insert_ub=500]
    # Returns (stdout): CLC URL of imported read set
    local r1_local="$1"
    local r2_local="$2"
    local dest_url="$3"
    local insert_lb="${4:-250}"
    local insert_ub="${5:-500}"

    _clc_check_env || return 1

    local job_id
    if [[ -n "${r2_local}" ]]; then
    	_clc_log "INFO" "Importing paired reads: ${r1_local} + ${r2_local} → ${dest_url}"
    	job_id=$(_clc_submit_async ngs_import_illumina \
                -f "${r1_local}" \
        	-f "${r2_local}" \
        	-d "${dest_url}" \
        	--paired-reads true \
        	--read-orientation FORWARD_REVERSE \
        	--min-distance "${insert_lb}" \
        	--max-distance "${insert_ub}")
    else
    	_clc_log "INFO" "IMporting single end reads: ${r1_local} -> ${dest_url}"
	job_id=$(_clc_submit_async ngs_import_illumina \
	       -f "${r1_local}" \
	       -d "${dest_url}" \
	       --paired-reads false)
    fi

    _clc_log "INFO" "Import job submitted: ${job_id}"
    clc_wait_job "${job_id}" || return 1
    _clc_retrieve "${job_id}" >/dev/null || true   # commits results to output folder

    local imported_url
    imported_url=$(_clc_ls_data_url "${dest_url}" "Sequence List") || true
    if [[ -z "${imported_url}" ]]; then
        _clc_log "ERROR" "Could not find Sequence List in ${dest_url} after import"
        return 1
    fi
    _clc_log "INFO" "Imported reads URL: ${imported_url}"
    echo "${imported_url}"
}

clc_import_long_reads() {
    # Import FASTQ reads to CLC Server using ngs_import_illumina.
    # Insert sizes are stored in the read metadata and used by downstream analysis.
    # Usage: clc_import_long_reads <r_local> <dest_clc_url>
    # Returns (stdout): CLC URL of imported read set
    local r_local="$1"
    local dest_url="$2"

    _clc_check_env || return 1

    local job_id
    _clc_log "INFO" "IMporting Long end reads: ${r_local} -> ${dest_url}"
	job_id=$(_clc_submit_async ngs_import_nanopore \
	       -f "${r_local}" \
	       -d "${dest_url}")

    _clc_log "INFO" "Import job submitted: ${job_id}"
    clc_wait_job "${job_id}" || return 1
    _clc_retrieve "${job_id}" >/dev/null || true   # commits results to output folder

    local imported_url
    imported_url=$(_clc_ls_data_url "${dest_url}" "Sequence List") || true
    if [[ -z "${imported_url}" ]]; then
        _clc_log "ERROR" "Could not find Sequence List in ${dest_url} after import"
        return 1
    fi
    _clc_log "INFO" "Imported reads URL: ${imported_url}"
    echo "${imported_url}"
}

clc_import_fasta() {
    # Import a FASTA file (reference/contigs) to CLC Server using ngs_import_fasta.
    # Usage: clc_import_fasta <local_fasta> <dest_clc_url>
    # Returns (stdout): CLC URL of imported sequence object
    local fasta_local="$1"
    local dest_url="$2"

    _clc_check_env || return 1
    _clc_log "INFO" "Importing FASTA: ${fasta_local} → ${dest_url}"

    local job_id
    job_id=$(_clc_submit_async ngs_import_fasta \
        -f "${fasta_local}" \
        -d "${dest_url}")

    _clc_log "INFO" "FASTA import job submitted: ${job_id}"
    clc_wait_job "${job_id}" || return 1
    _clc_retrieve "${job_id}" >/dev/null || true   # commits results to output folder

    local imported_url
    imported_url=$(_clc_ls_data_url "${dest_url}") || true
    if [[ -z "${imported_url}" ]]; then
        _clc_log "ERROR" "Could not find imported sequence in ${dest_url} after FASTA import"
        return 1
    fi
    _clc_log "INFO" "Imported FASTA URL: ${imported_url}"
    echo "${imported_url}"
}

# =============================================================================
# Export (download server objects to local files)
# Export runs synchronously but still requires CLC_GRID.
# =============================================================================

clc_export_fastq() {
    # Export a CLC read object as FASTQ to a local directory.
    # Usage: clc_export_fastq <server_url> <local_dir>
    local server_url="$1"
    local local_dir="$2"

    _clc_check_env || return 1
    mkdir -p "${local_dir}"
    _clc_log "INFO" "Exporting FASTQ: ${server_url} → ${local_dir}/"

    local cmd_arr
    _clc_conn_args cmd_arr
    cmd_arr+=(
        -A export
        -e fastq
        -i "${server_url}"
        -d "${local_dir}"
    )
    _clc_log "DEBUG" "Running ${cmd_arr[*]}"
    "${cmd_arr[@]}" </dev/null
}

clc_export_fasta() {
    # Export a CLC sequence object as FASTA to a local directory.
    # Usage: clc_export_fasta <server_url> <local_dir>
    local server_url="$1"
    local local_dir="$2"

    _clc_check_env || return 1
    mkdir -p "${local_dir}"
    _clc_log "INFO" "Exporting FASTA: ${server_url} → ${local_dir}/"

    local cmd_arr
    _clc_conn_args cmd_arr
    cmd_arr+=(
        -A export
        -e fasta
        -i "${server_url}"
        -d "${local_dir}"
    )
    _clc_log "DEBUG" "Running ${cmd_arr[*]}"
    "${cmd_arr[@]}" </dev/null
}

split_interleaved_fastq() {
    # Split an interleaved FASTQ file (from CLC FASTQ export) into R1 and R2.
    # Usage: split_interleaved_fastq <local_dir> <out_r1> <out_r2>
    local local_dir="$1"
    local out_r1="$2"
    local out_r2="$3"

    local interleaved_file
    interleaved_file=$(find "${local_dir}" -maxdepth 1 \( -name "*.fastq" -o -name "*.fq" \) 2>/dev/null | head -1)

    if [[ -z "${interleaved_file}" ]]; then
        _clc_log "ERROR" "No FASTQ file found in ${local_dir} to split"
        return 1
    fi

    _clc_log "INFO" "Splitting interleaved FASTQ: ${interleaved_file} → ${out_r1} + ${out_r2}"

    python3 - <<PYEOF
interleaved = "${interleaved_file}"
out_r1 = "${out_r1}"
out_r2 = "${out_r2}"

def iter_fastq(fh):
    while True:
        h = fh.readline()
        if not h:
            break
        s = fh.readline()
        p = fh.readline()
        q = fh.readline()
        if not (h and s and p and q):
            break
        yield h, s, p, q

with open(interleaved) as fh, \
     open(out_r1, 'w') as r1_fh, \
     open(out_r2, 'w') as r2_fh:
    i = 0
    for record in iter_fastq(fh):
        if i % 2 == 0:
            r1_fh.writelines(record)
        else:
            r2_fh.writelines(record)
        i += 1

print(f"Split {i} reads: {i//2} pairs into R1 and R2")
PYEOF

    if [[ ! -s "${out_r1}" ]] || [[ ! -s "${out_r2}" ]]; then
        _clc_log "ERROR" "Splitting failed: R1 or R2 is missing/empty"
        return 1
    fi
    _clc_log "INFO" "Split complete: $(wc -l < "${out_r1}") lines in R1, $(wc -l < "${out_r2}") lines in R2"
}

# =============================================================================
# CLC Tool wrappers
# All wrappers use _clc_submit_async (adds -G when CLC_GRID is set).
# =============================================================================

clc_trim_reads() {
    # Run CLC quality trimming on paired reads on the server.
    # Usage: clc_trim_reads <reads_url> <output_folder_url> [quality_cutoff=30] [min_length=50]
    # Returns (stdout): job ID
    #
    # quality_cutoff is Phred score; converted to error probability for --quality-limit.
    # Q30 → 0.001, Q20 → 0.01
    local reads_url="$1"
    local output_url="$2"
    local quality_cutoff="${3:-30}"
    local min_length="${4:-50}"

    local quality_limit
    quality_limit=$(python3 -c "print(10**(-${quality_cutoff}/10))")

    _clc_log "INFO" "Trimming: Q${quality_cutoff} (limit=${quality_limit}), min_L=${min_length}"

    _clc_submit_async trim \
        -i "${reads_url}" \
        -d "${output_url}" \
        --quality-limit "${quality_limit}" \
        --quality-trim true \
        --discard-short true \
        --discard-short-limit "${min_length}" \
        --ambiguous-limit 2 \
        --ambiguous-trim true
}

clc_assemble() {
    # Run CLC de novo assembly on the server.
    # Insert sizes are stored in the read metadata (set during ngs_import_illumina);
    # auto-detect is used here.
    # Usage: clc_assemble <reads_url> <output_folder_url> [word_size=64] [bubble_size=200]
    # Returns (stdout): job ID
    local reads_url="$1"
    local output_url="$2"
    local word_size="${3:-64}"
    local bubble_size="${4:-200}"

    _clc_log "INFO" "Assembling: wordsize=${word_size}, bubblesize=${bubble_size}"

    _clc_submit_async denovo_assembly \
        -i "${reads_url}" \
        -d "${output_url}" \
	--map-reads-to-contigs SIMPLE \
        --wordsize "${word_size}" \
        --wordsize-auto false \
        --bubblesize "${bubble_size}" \
        --bubblesize-auto false \
        --minimum-contig-length 200 \
        --auto-detect-paired-distances true \
        --perform-scaffolding false
}

clc_assemble_long() {
    # Run CLC LONG de novo assembly on the server.
    # Usage: clc_assemble_long <reads_url> <output_folder_url>
    # Returns (stdout): job ID
    local reads_url="$1"
    local output_url="$2"

    _clc_log "INFO" "Assembling: long reads"

    _clc_submit_async long_de_novo_assembly \
        -i "${reads_url}" \
        -d "${output_url}"
}

clc_map_reads() {
    # Map reads to a reference on the server.
    # Reference is passed via --references (not -i).
    # --collect-unmapped true includes unmapped reads in the output.
    # --output-mode CLUSTER produces stand-alone mapping objects needed for consensus extraction.
    # Usage: clc_map_reads <reads_url> <reference_url> <output_folder_url>
    # Returns (stdout): job ID
    local reads_url="$1"
    local reference_url="$2"
    local output_url="$3"

    _clc_log "INFO" "Mapping reads to reference (global alignment)"

    _clc_submit_async read_mapping \
        -i "${reads_url}" \
        --references "${reference_url}" \
        -d "${output_url}" \
        --global-alignment true \
        --auto-detect-paired-distances true \
        --collect-unmapped true \
        --output-mode CLUSTER
}

clc_extract_consensus() {
    # Extract consensus sequence from a CLC mapping object.
    # Usage: clc_extract_consensus <mapping_url> <output_folder_url>
    # Returns (stdout): job ID
    local mapping_url="$1"
    local output_url="$2"

    _clc_submit_async consensus_sequence_extraction \
        -i "${mapping_url}" \
        -d "${output_url}"
}

clc_export_unmapped_reads() {
    # Export unmapped reads from a read_mapping result using the second URL returned by
    # _clc_retrieve (read_mapping with --collect-unmapped true emits two objects:
    #   [0] the read mapping itself, [1] the unmapped sequence list).
    # Usage: clc_export_unmapped_reads <unmapped_clc_url> <local_dir>
    local unmapped_url="$1"
    local local_dir="$2"

    _clc_check_env || return 1
    mkdir -p "${local_dir}"
    _clc_log "INFO" "Exporting unmapped reads: ${unmapped_url} → ${local_dir}/"

    local cmd_arr
    _clc_conn_args cmd_arr
    cmd_arr+=(
        -A export
        -e fasta
        -i "${unmapped_url}"
        -d "${local_dir}"
    )
    _clc_log "DEBUG" "Running ${cmd_arr[*]}"
    "${cmd_arr[@]}" </dev/null
}

clc_mapping_stats() {
    # Get read mapping statistics from a CLC mapping object.
    # Runs synchronously; output written to stdout.
    # Usage: clc_mapping_stats <mapping_url>
    local mapping_url="$1"

    local cmd_arr
    _clc_conn_args cmd_arr
    cmd_arr+=(
        -A detailed_mapping_report
        -i "${mapping_url}"
    )
    _clc_log "DEBUG" "Running ${cmd_arr[*]}"
    "${cmd_arr[@]}" </dev/null
}

# =============================================================================
# High-level composite pipeline functions
# =============================================================================

clc_qc_pipeline() {
    # Complete CLC QC pipeline: upload → trim → export R1+R2 → cleanup
    #
    # Usage:
    #   clc_qc_pipeline <r1_local> <r2_local> <phage_name> <output_dir>
    #                   [quality_cutoff=30] [min_length=50]
    #                   [insert_lb=250] [insert_ub=500]
    #
    # Outputs (in <output_dir>/):
    #   <phage_name>.1.clc-trimmed.fastq
    #   <phage_name>.2.clc-trimmed.fastq
    #   <phage_name>.clc-qc.stats.txt

    local r1_local="$1"
    local r2_local="$2"
    local phage_name="$3"
    local output_dir="$4"
    local quality_cutoff="${5:-30}"
    local min_length="${6:-50}"
    local insert_lb="${7:-250}"
    local insert_ub="${8:-500}"

    _clc_check_env || return 1
    mkdir -p "${output_dir}"

    local work_url="${CLC_DATA_ROOT}/${phage_name}/qc"

    _clc_log "INFO" "=== CLC QC Pipeline: ${phage_name} ==="

    # 1. Create server workspace
    # clc_ensure_dir for the phage folder (never re-create if it exists — CLC would make phage-1)
    # clc_rm + clc_mkdir for the qc subfolder to overwrite any previous run cleanly
    _clc_log "INFO" "Step 1: Creating server workspace"
    clc_ensure_dir "${CLC_DATA_ROOT}" "${phage_name}"
    clc_rm "${work_url}"
    clc_mkdir "${CLC_DATA_ROOT}/${phage_name}" "qc"

    # 2. Import reads
    _clc_log "INFO" "Step 2: Importing raw reads"
    local reads_url="${work_url}/raw_reads"
    clc_mkdir "${work_url}" "raw_reads"
    local imported_url
    imported_url=$(clc_import_paired_reads "${r1_local}" "${r2_local}" "${reads_url}" \
        "${insert_lb}" "${insert_ub}")

    # Count input reads for QC percentage
    local in_r1_reads=0
    if [[ "${r1_local}" == *.gz ]]; then
        in_r1_reads=$(zcat "${r1_local}" | awk 'END{print NR/4}' 2>/dev/null || echo 0)
    else
        in_r1_reads=$(awk 'END{print NR/4}' "${r1_local}" 2>/dev/null || echo 0)
    fi
    local total_input=$(( in_r1_reads * 2 ))

    # 3. Trim reads
    _clc_log "INFO" "Step 3: Trimming reads"
    local trim_url="${work_url}/trimmed"
    clc_mkdir "${work_url}" "trimmed"
    local trim_job
    trim_job=$(clc_trim_reads "${imported_url}" "${trim_url}" "${quality_cutoff}" "${min_length}")
    clc_wait_job "${trim_job}" || { _clc_log "ERROR" "Trim job failed"; clc_rm "${work_url}"; return 1; }
    _clc_retrieve "${trim_job}" >/dev/null || true   # commits results to output folder

    # Find trimmed reads object via ls
    local trimmed_url
    trimmed_url=$(_clc_ls_data_url "${trim_url}" "Sequence List") || true
    if [[ -z "${trimmed_url}" ]]; then
        _clc_log "WARN" "No trimmed Sequence List in ${trim_url} — reads may all be high quality; using raw reads"
        trimmed_url="${imported_url}"
    fi
    _clc_log "INFO" "Trimmed reads URL: ${trimmed_url}"

    # 4. Export trimmed reads
    _clc_log "INFO" "Step 4: Exporting trimmed reads"
    local export_dir="${output_dir}/_clc_export_tmp"
    clc_export_fastq "${trimmed_url}" "${export_dir}"

    # 5. Split interleaved export into R1 + R2
    _clc_log "INFO" "Step 5: Splitting interleaved FASTQ"
    local out_r1="${output_dir}/${phage_name}.1.clc-trimmed.fastq"
    local out_r2="${output_dir}/${phage_name}.2.clc-trimmed.fastq"
    split_interleaved_fastq "${export_dir}" "${out_r1}" "${out_r2}"
    rm -rf "${export_dir}"

    # 6. Write stats file
    _clc_log "INFO" "Step 6: Writing QC stats"
    local r1_reads; r1_reads=$(awk 'END{print NR/4}' "${out_r1}")
    local total_reads=$(( r1_reads * 2 ))
    local pct_str="100.00"
    if [[ "${total_input}" -gt 0 ]]; then
        pct_str=$(python3 -c "print(f'{100*${total_reads}/${total_input}:.2f}')")
    fi
    echo "Output reads: ${total_reads} (${pct_str}%)" > "${output_dir}/${phage_name}.clc-qc.stats.txt"

    # 7. Cleanup server workspace
    _clc_log "INFO" "Step 7: Cleaning up server workspace"
    clc_rm "${work_url}"

    _clc_log "INFO" "=== CLC QC Pipeline complete: ${total_reads} reads ==="
}


clc_assembly_pipeline() {
    # Complete CLC assembly pipeline: upload → assemble → export contigs → cleanup
    #
    # Usage:
    #   clc_assembly_pipeline <r1_local> <r2_local> <phage_name> <output_dir>
    #                         [word_size=64] [bubble_size=200] [insert_lb=250] [insert_ub=500]
    #
    # Output: <output_dir>/<phage_name>_CLC_contigs.fasta

    local r1_local="$1"
    local r2_local="$2"
    local phage_name="$3"
    local output_dir="$4"
    local word_size="${5:-64}"
    local bubble_size="${6:-200}"
    local insert_lb="${7:-250}"
    local insert_ub="${8:-500}"

    _clc_check_env || return 1
    mkdir -p "${output_dir}"

    local work_url="${CLC_DATA_ROOT}/${phage_name}/assembly"

    _clc_log "INFO" "=== CLC Assembly Pipeline: ${phage_name} ==="

    # 1. Create workspace
    # clc_ensure_dir for the phage folder; clc_rm + clc_mkdir for the assembly subfolder
    _clc_log "INFO" "Step 1: Creating server workspace"
    clc_ensure_dir "${CLC_DATA_ROOT}" "${phage_name}"
    clc_rm "${work_url}"
    clc_mkdir "${CLC_DATA_ROOT}/${phage_name}" "assembly"

    # 2. Import reads (insert sizes stored in read metadata for auto-detect)
    _clc_log "INFO" "Step 2: Importing reads"
    local reads_url="${work_url}/reads"
    clc_mkdir "${work_url}" "reads"
    local imported_url
    imported_url=$(clc_import_paired_reads "${r1_local}" "${r2_local}" "${reads_url}" \
        "${insert_lb}" "${insert_ub}")

    # 3. Assemble
    _clc_log "INFO" "Step 3: Running assembly"
    local asm_url="${work_url}/assembly"
    clc_mkdir "${work_url}" "assembly"
    local asm_job
    asm_job=$(clc_assemble "${imported_url}" "${asm_url}" "${word_size}" "${bubble_size}")
    clc_wait_job "${asm_job}" || { _clc_log "ERROR" "Assembly job failed"; clc_rm "${work_url}"; return 1; }
    _clc_retrieve "${asm_job}" >/dev/null || true   # commits results to output folder

    # Find contig sequence object via ls
    local asm_result_url
    asm_result_url=$(_clc_ls_data_url "${asm_url}" "Sequence List") || true
    if [[ -z "${asm_result_url}" ]]; then
        _clc_log "ERROR" "Could not find assembly result in ${asm_url}"
        clc_rm "${work_url}"; return 1
    fi
    _clc_log "INFO" "Assembly result URL: ${asm_result_url}"

    # 4. Export contigs as FASTA
    _clc_log "INFO" "Step 4: Exporting contigs"
    local export_dir="${output_dir}/_clc_asm_export_tmp"
    clc_export_fasta "${asm_result_url}" "${export_dir}"

    local exported_fasta
    exported_fasta=$(find "${export_dir}" -maxdepth 1 \( -name "*.fasta" -o -name "*.fa" \) 2>/dev/null | head -1)
    if [[ -z "${exported_fasta}" ]]; then
        _clc_log "ERROR" "No FASTA file found in export directory ${export_dir}"
        clc_rm "${work_url}"
        return 1
    fi
    mv "${exported_fasta}" "${output_dir}/${phage_name}_CLC_contigs.fasta"
    rm -rf "${export_dir}"

    # 5. Cleanup server
    _clc_log "INFO" "Step 5: Cleaning up server workspace"
    clc_rm "${work_url}"

    _clc_log "INFO" "=== CLC Assembly Pipeline complete ==="
}

clc_long_assembly_pipeline() {
    # Complete LOng Reads CLC assembly pipeline: upload → assemble → export contigs → cleanup
    #
    # Usage:
    #   clc_long_assembly_pipeline <r_local> <phage_name> <output_dir>
    #
    # Output: <output_dir>/<phage_name>_CLC_contigs.fasta

    local r_local="$1"
    local phage_name="$2"
    local output_dir="$3"


    _clc_check_env || return 1
    mkdir -p "${output_dir}"

    local work_url="${CLC_DATA_ROOT}/${phage_name}/assembly"

    _clc_log "INFO" "=== CLC Assembly Pipeline: ${phage_name} ==="

    # 1. Create workspace
    # clc_ensure_dir for the phage folder; clc_rm + clc_mkdir for the assembly subfolder
    _clc_log "INFO" "Step 1: Creating server workspace"
    clc_ensure_dir "${CLC_DATA_ROOT}" "${phage_name}"
    clc_rm "${work_url}"
    clc_mkdir "${CLC_DATA_ROOT}/${phage_name}" "assembly"

    # 2. Import reads (insert sizes stored in read metadata for auto-detect)
    _clc_log "INFO" "Step 2: Importing reads"
    local reads_url="${work_url}/reads"
    clc_mkdir "${work_url}" "reads"
    local imported_url
    imported_url=$(clc_import_long_reads "${r_local}" "${reads_url}")

    # 3. Assemble
    _clc_log "INFO" "Step 3: Running assembly"
    local asm_url="${work_url}/assembly"
    clc_mkdir "${work_url}" "assembly"
    local asm_job
    asm_job=$(clc_assemble_long "${imported_url}" "${asm_url}")
    clc_wait_job "${asm_job}" || { _clc_log "ERROR" "Assembly job failed"; clc_rm "${work_url}"; return 1; }
    _clc_retrieve "${asm_job}" >/dev/null || true   # commits results to output folder

    # Find contig sequence object via ls
    local asm_result_url
    asm_result_url=$(_clc_ls_data_url "${asm_url}" "Sequence List") || true
    if [[ -z "${asm_result_url}" ]]; then
        _clc_log "ERROR" "Could not find assembly result in ${asm_url}"
        clc_rm "${work_url}"; return 1
    fi
    _clc_log "INFO" "Assembly result URL: ${asm_result_url}"

    # 4. Export contigs as FASTA
    _clc_log "INFO" "Step 4: Exporting contigs"
    local export_dir="${output_dir}/_clc_asm_export_tmp"
    clc_export_fasta "${asm_result_url}" "${export_dir}"

    local exported_fasta
    exported_fasta=$(find "${export_dir}" -maxdepth 1 \( -name "*.fasta" -o -name "*.fa" \) 2>/dev/null | head -1)
    if [[ -z "${exported_fasta}" ]]; then
        _clc_log "ERROR" "No FASTA file found in export directory ${export_dir}"
        clc_rm "${work_url}"
        return 1
    fi
    mv "${exported_fasta}" "${output_dir}/${phage_name}_CLC_contigs.fasta"
    rm -rf "${export_dir}"

    # 5. Cleanup server
    _clc_log "INFO" "Step 5: Cleaning up server workspace"
    clc_rm "${work_url}"

    _clc_log "INFO" "=== CLC Long Reads Assembly Pipeline complete ==="
}

clc_mapping_pipeline() {
    # Complete CLC mapping pipeline:
    #   upload reads + reference → map (with unmapped collection)
    #   → stats → extract consensus → export → cleanup
    #
    # Usage:
    #   clc_mapping_pipeline <r1_local> <r2_local> <reference_fasta>
    #                        <phage_name> <output_dir>
    #                        [insert_lb=250] [insert_ub=500]
    #
    # Outputs (in <output_dir>/):
    #   <phage_name>_corrected_contigs.fasta
    #   <phage_name>_unmapped_reads.fasta  (empty stub — unmapped reads are in mapping output)
    #   <phage_name>_mapping_stats.txt

    local r1_local="$1"
    local r2_local="$2"
    local reference_fasta="$3"
    local phage_name="$4"
    local output_dir="$5"
    local insert_lb="${6:-250}"
    local insert_ub="${7:-500}"

    _clc_check_env || return 1
    mkdir -p "${output_dir}"

    local work_url="${CLC_DATA_ROOT}/${phage_name}/mapping"

    _clc_log "INFO" "=== CLC Mapping Pipeline: ${phage_name} ==="

    # 1. Create workspace
    # clc_ensure_dir for the phage folder; clc_rm + clc_mkdir for the mapping subfolder
    _clc_log "INFO" "Step 1: Creating server workspace"
    clc_ensure_dir "${CLC_DATA_ROOT}" "${phage_name}"
    clc_rm "${work_url}"
    clc_mkdir "${CLC_DATA_ROOT}/${phage_name}" "mapping"

    # 2. Import reads
    _clc_log "INFO" "Step 2: Importing reads"
    local reads_url="${work_url}/reads"
    clc_mkdir "${work_url}" "reads"
    local reads_imported
    reads_imported=$(clc_import_paired_reads "${r1_local}" "${r2_local}" "${reads_url}" \
        "${insert_lb}" "${insert_ub}")

    # 3. Import reference (FASTA contigs)
    _clc_log "INFO" "Step 3: Importing reference"
    local ref_url="${work_url}/reference"
    clc_mkdir "${work_url}" "reference"
    local ref_imported
    ref_imported=$(clc_import_fasta "${reference_fasta}" "${ref_url}")

    # 4. Map reads to reference
    _clc_log "INFO" "Step 4: Mapping reads"
    local map_url="${work_url}/mapping"
    clc_mkdir "${work_url}" "mapping"
    local map_job
    map_job=$(clc_map_reads "${reads_imported}" "${ref_imported}" "${map_url}")
    clc_wait_job "${map_job}" || { _clc_log "ERROR" "Mapping job failed"; clc_rm "${work_url}"; return 1; }
    _clc_retrieve "${map_job}" >/dev/null || true   # commits results to output folder

    # Find objects in mapping output folder via ls.
    # read_mapping --collect-unmapped true produces two objects:
    #   Read Mapping  → used for consensus extraction and stats
    #   Sequence List → unmapped reads
    local map_result_url
    map_result_url=$(_clc_ls_data_url "${map_url}" "Read Mapping") || true
    local unmapped_result_url
    unmapped_result_url=$(_clc_ls_data_url "${map_url}" "Sequence List") || true
    if [[ -z "${map_result_url}" ]]; then
        _clc_log "ERROR" "Could not find Read Mapping object in ${map_url}"
        clc_rm "${work_url}"; return 1
    fi
    _clc_log "INFO" "Mapping result URL: ${map_result_url}"
    _clc_log "INFO" "Unmapped reads URL: ${unmapped_result_url:-<none>}"

    # 5. Mapping statistics
    _clc_log "INFO" "Step 5: Collecting mapping statistics"
    clc_mapping_stats "${map_result_url}" \
        > "${output_dir}/${phage_name}_mapping_stats.txt" 2>/dev/null || \
        echo "Mapped reads 0 (0.00%)" > "${output_dir}/${phage_name}_mapping_stats.txt"

    # 6. Extract consensus (corrected contigs)
    _clc_log "INFO" "Step 6: Extracting consensus"
    local cons_url="${work_url}/consensus"
    clc_mkdir "${work_url}" "consensus"
    local cons_job
    cons_job=$(clc_extract_consensus "${map_result_url}" "${cons_url}")
    clc_wait_job "${cons_job}" || _clc_log "WARN" "Consensus extraction may have failed"
    _clc_retrieve "${cons_job}" >/dev/null || true   # commits results to output folder

    local cons_result_url
    cons_result_url=$(_clc_ls_data_url "${cons_url}") || true
    [[ -z "${cons_result_url}" ]] && { _clc_log "WARN" "Could not find consensus object in ${cons_url}"; cons_result_url="${cons_url}"; }

    local cons_export_dir="${output_dir}/_cons_tmp"
    clc_export_fasta "${cons_result_url}" "${cons_export_dir}" || true
    local cons_fasta
    cons_fasta=$(find "${cons_export_dir}" -maxdepth 1 \( -name "*.fasta" -o -name "*.fa" \) 2>/dev/null | head -1)
    if [[ -n "${cons_fasta}" ]]; then
        mv "${cons_fasta}" "${output_dir}/${phage_name}_corrected_contigs.fasta"
    else
        _clc_log "WARN" "Consensus export failed; copying original reference as fallback"
        cp "${reference_fasta}" "${output_dir}/${phage_name}_corrected_contigs.fasta"
    fi
    rm -rf "${cons_export_dir}"

    # 7. Unmapped reads — second output object from read_mapping (--collect-unmapped true)
    _clc_log "INFO" "Step 7: Exporting unmapped reads"
    if [[ -n "${unmapped_result_url}" ]]; then
        local unmap_export_dir="${output_dir}/_unmap_tmp"
        clc_export_unmapped_reads "${unmapped_result_url}" "${unmap_export_dir}" || true
        local unmap_fasta
        unmap_fasta=$(find "${unmap_export_dir}" -maxdepth 1 \( -name "*.fasta" -o -name "*.fa" \) 2>/dev/null | head -1)
        if [[ -n "${unmap_fasta}" ]]; then
            mv "${unmap_fasta}" "${output_dir}/${phage_name}_unmapped_reads.fasta"
        else
            _clc_log "WARN" "No unmapped reads exported (possibly all reads mapped)"
            touch "${output_dir}/${phage_name}_unmapped_reads.fasta"
        fi
        rm -rf "${unmap_export_dir}"
    else
        _clc_log "WARN" "No unmapped reads URL in mapping result (possibly all reads mapped)"
        touch "${output_dir}/${phage_name}_unmapped_reads.fasta"
    fi

    # 8. Cleanup server workspace
    _clc_log "INFO" "Step 8: Cleaning up server workspace"
    clc_rm "${work_url}"

    _clc_log "INFO" "=== CLC Mapping Pipeline complete ==="
}


# =============================================================================
# Utility
# =============================================================================

clc_rm_phage() {
    # Delete the entire phage folder from the CLC server (use after a pipeline run or to reset).
    # Usage: clc_rm_phage <phage_name>
    local phage_name="$1"
    _clc_check_env || return 1
    clc_rm "${CLC_DATA_ROOT}/${phage_name}"
}

clc_preflight() {
    # Test connectivity and verify CLC_DATA_ROOT is accessible.
    _clc_check_env || return 1
    _clc_log "INFO" "Testing CLC Server connection: ${CLC_SERVER_HOST}:${CLC_SERVER_PORT:-7777}"
    local ls_arr=(
        "${CLC_SERVER_CLI}"
        -S "${CLC_SERVER_HOST}"
        -P "${CLC_SERVER_PORT:-7777}"
        -U "${CLC_SERVER_USER}"
        -W "${CLC_SERVER_PASS}"
        -A ls
        -t "${CLC_DATA_ROOT}"
    )
    _clc_log "DEBUG" "Running ${ls_arr[*]}"
    "${ls_arr[@]}" </dev/null 2>&1 | head -20 || {
        _clc_log "ERROR" "Connection test failed. Check CLC_SERVER_HOST/PORT/USER/PASS/DATA_ROOT."
        return 1
    }
    _clc_log "INFO" "Connection OK. CLC_GRID=${CLC_GRID:-<not set>}"
}
