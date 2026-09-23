#!/usr/bin/env bash
set -euo pipefail

# Submit ge_screen jobs for multiple assemblies with matching read files.
#
# Behavior:
#   - multiple assemblies from the same sample are processed in separate directories:
#       <outdir>/<sample>/<assembly_stem>/
#   - captures processing job IDs from ge_screen.sh submission output
#   - optionally submits dashboard generation with a SLURM dependency so it waits
#     for all processing jobs to finish
#
# Dashboard dependency behavior:
#   --submit-dashboard
#      submits dashboard job after all processing jobs are submitted
#   --dashboard-dependency afterok      default; dashboard runs only if all jobs succeed
#   --dashboard-dependency afterany     dashboard runs after all jobs finish, even failed jobs
#   --dashboard-dependency singleton    dashboard uses singleton dependency instead of job IDs
#
# Example:
#   bash run_ge_screen_improved_v3.sh \
#     --assembly-dir /path/to/all_assemblies \
#     --reads-dir /path/to/raw_reads/all_barcodes \
#     --outdir /path/to/ge_screen \
#     --partition ss \
#     --threads 16 \
#     --mem 64G \
#     --submit-dashboard \
#     --dashboard-dependency afterok

usage() {
  cat <<'EOF'
Submit ge_screen jobs for multiple assembly/read pairs.

Defaults:
  --scriptdir                  /mnt/genomics/common_projects/COMMEX_Exercise/scripts/ge_screen
  --ge-screen-script           <scriptdir>/ge_screen.sh
  --assembly-dir               current directory
  --reads-dir                  current directory
  --outdir                     ge_screen_results
  --patterns                   *_dragonflye.fasta,*_dragonflye_contigs.fasta,*_EDGE.fasta
  --read-ext                   .fastq.gz
  --partition                  ss
  --threads                    16
  --mem                        64G
  --time                       not passed to ge_screen.sh unless set
  --minimap2-preset            auto
  --amrfinderplus              true
  --amrfinder-db               not set
  --amrfinder-organism         not set
  --amrfinder-annotation-format prodigal
  --dry-run                    false
  --submit-dashboard           false
  --dashboard-dependency       afterok
  --dashboard-embed-plotly     false
  --dashboard-max-nt-rows      100
  --dashboard-max-univec-rows  100
  --dashboard-max-coverage-rows 250
  --dashboard-max-abricate-rows 500
  --dashboard-max-amrfinder-rows 500

Options:
  --scriptdir DIR
  --ge-screen-script FILE
  --assembly-dir DIR
  --reads-dir DIR
  --outdir DIR
  --patterns CSV
  --read-ext EXT
  --partition NAME
  --threads N
  --mem MEM
  --time HH:MM:SS
  --minimap2-preset PRESET
  --amrfinderplus
  --no-amrfinderplus
  --amrfinder-db DIR
  --amrfinder-organism NAME
  --amrfinder-annotation-format FORMAT
  --submit-dashboard
  --dashboard-dependency MODE  afterok|afterany|singleton|none
  --dashboard-script FILE
  --dashboard-wrapper FILE
  --dashboard-output FILE
  --dashboard-summary FILE
  --dashboard-embed-plotly
  --dashboard-max-nt-rows N
  --dashboard-max-univec-rows N
  --dashboard-max-coverage-rows N
  --dashboard-max-abricate-rows N
  --dashboard-max-amrfinder-rows N
  --dry-run
  -h, --help

Output layout:
  <outdir>/<sample>/<assembly_stem>/

Read matching:
  Assembly suffixes removed before appending --read-ext:
    _dragonflye_contigs.fasta
    _dragonflye.fasta
    _EDGE.fasta
Example:
  ST3-2_0001_DODSAFE_SS_fastcat_barcode0013_dragonflye.fasta
  -> read stem:
  ST3-2_0001_DODSAFE_SS_fastcat_barcode0013
  -> expected read:
  ST3-2_0001_DODSAFE_SS_fastcat_barcode0013.fastq.gz
EOF
}

SCRIPT_DIR="/mnt/genomics/common_projects/COMMEX_Exercise/scripts/ge_screen"
GE_SCREEN_SCRIPT=""
ASSEMBLY_DIR="."
READS_DIR="."
OUTDIR="ge_screen_results"
PATTERNS="*_dragonflye.fasta,*_dragonflye_contigs.fasta,*_EDGE.fasta"
READ_EXT=".fastq.gz"
PARTITION="ss"
THREADS="16"
MEM="64G"
TIME_LIMIT=""
MINIMAP2_PRESET="auto"
RUN_AMRFINDER="true"
AMRFINDER_DB=""
AMRFINDER_ORGANISM=""
AMRFINDER_ANNOTATION_FORMAT="prodigal"
DRY_RUN="false"
SUBMIT_DASHBOARD="false"
DASHBOARD_DEPENDENCY="afterok"
DASHBOARD_SCRIPT=""
DASHBOARD_WRAPPER=""
DASHBOARD_OUTPUT=""
DASHBOARD_SUMMARY=""
DASHBOARD_EMBED_PLOTLY="false"
DASHBOARD_MAX_NT_ROWS="100"
DASHBOARD_MAX_UNIVEC_ROWS="100"
DASHBOARD_MAX_COVERAGE_ROWS="250"
DASHBOARD_MAX_ABRICATE_ROWS="500"
DASHBOARD_MAX_AMRFINDER_ROWS="500"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --scriptdir) SCRIPT_DIR="${2:?}"; shift 2 ;;
    --ge-screen-script) GE_SCREEN_SCRIPT="${2:?}"; shift 2 ;;
    --assembly-dir) ASSEMBLY_DIR="${2:?}"; shift 2 ;;
    --reads-dir) READS_DIR="${2:?}"; shift 2 ;;
    --outdir) OUTDIR="${2:?}"; shift 2 ;;
    --patterns) PATTERNS="${2:?}"; shift 2 ;;
    --read-ext) READ_EXT="${2:?}"; shift 2 ;;
    --partition) PARTITION="${2:?}"; shift 2 ;;
    --threads) THREADS="${2:?}"; shift 2 ;;
    --mem) MEM="${2:?}"; shift 2 ;;
    --time) TIME_LIMIT="${2:?}"; shift 2 ;;
    --minimap2-preset) MINIMAP2_PRESET="${2:?}"; shift 2 ;;
    --amrfinderplus) RUN_AMRFINDER="true"; shift ;;
    --no-amrfinderplus) RUN_AMRFINDER="false"; shift ;;
    --amrfinder-db) AMRFINDER_DB="${2:?}"; shift 2 ;;
    --amrfinder-organism) AMRFINDER_ORGANISM="${2:?}"; shift 2 ;;
    --amrfinder-annotation-format) AMRFINDER_ANNOTATION_FORMAT="${2:?}"; shift 2 ;;
    --submit-dashboard) SUBMIT_DASHBOARD="true"; shift ;;
    --dashboard-dependency) DASHBOARD_DEPENDENCY="${2:?}"; shift 2 ;;
    --dashboard-script) DASHBOARD_SCRIPT="${2:?}"; shift 2 ;;
    --dashboard-wrapper) DASHBOARD_WRAPPER="${2:?}"; shift 2 ;;
    --dashboard-output) DASHBOARD_OUTPUT="${2:?}"; shift 2 ;;
    --dashboard-summary) DASHBOARD_SUMMARY="${2:?}"; shift 2 ;;
    --dashboard-embed-plotly) DASHBOARD_EMBED_PLOTLY="true"; shift ;;
    --dashboard-max-nt-rows) DASHBOARD_MAX_NT_ROWS="${2:?}"; shift 2 ;;
    --dashboard-max-univec-rows) DASHBOARD_MAX_UNIVEC_ROWS="${2:?}"; shift 2 ;;
    --dashboard-max-coverage-rows) DASHBOARD_MAX_COVERAGE_ROWS="${2:?}"; shift 2 ;;
    --dashboard-max-abricate-rows) DASHBOARD_MAX_ABRICATE_ROWS="${2:?}"; shift 2 ;;
    --dashboard-max-amrfinder-rows) DASHBOARD_MAX_AMRFINDER_ROWS="${2:?}"; shift 2 ;;
    --dry-run) DRY_RUN="true"; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

case "$DASHBOARD_DEPENDENCY" in
  afterok|afterany|singleton|none) ;;
  *)
    echo "ERROR: --dashboard-dependency must be one of: afterok, afterany, singleton, none" >&2
    exit 2
    ;;
esac

is_int() {
  [[ "$1" =~ ^-?[0-9]+$ ]]
}

for pair in \
  "dashboard-max-nt-rows:$DASHBOARD_MAX_NT_ROWS" \
  "dashboard-max-univec-rows:$DASHBOARD_MAX_UNIVEC_ROWS" \
  "dashboard-max-coverage-rows:$DASHBOARD_MAX_COVERAGE_ROWS" \
  "dashboard-max-abricate-rows:$DASHBOARD_MAX_ABRICATE_ROWS" \
  "dashboard-max-amrfinder-rows:$DASHBOARD_MAX_AMRFINDER_ROWS"
do
  name="${pair%%:*}"
  value="${pair#*:}"
  if ! is_int "$value"; then
    echo "ERROR: --$name must be an integer. Use -1 for all rows." >&2
    exit 2
  fi
done

ASSEMBLY_DIR="$(realpath "$ASSEMBLY_DIR")"
READS_DIR="$(realpath "$READS_DIR")"
OUTDIR="$(realpath -m "$OUTDIR")"
SCRIPT_DIR="$(realpath "$SCRIPT_DIR")"

if [[ -z "$GE_SCREEN_SCRIPT" ]]; then
  GE_SCREEN_SCRIPT="$SCRIPT_DIR/ge_screen.sh"
fi
GE_SCREEN_SCRIPT="$(realpath "$GE_SCREEN_SCRIPT")"

if [[ ! -d "$ASSEMBLY_DIR" ]]; then
  echo "ERROR: assembly directory not found: $ASSEMBLY_DIR" >&2
  exit 2
fi

if [[ ! -d "$READS_DIR" ]]; then
  echo "ERROR: reads directory not found: $READS_DIR" >&2
  exit 2
fi

if [[ ! -f "$GE_SCREEN_SCRIPT" ]]; then
  echo "ERROR: ge_screen script not found: $GE_SCREEN_SCRIPT" >&2
  exit 2
fi

mkdir -p "$OUTDIR"

assembly_stem_from_name() {
  local name="$1"
  case "$name" in
    *_dragonflye_contigs.fasta) echo "${name%_dragonflye_contigs.fasta}" ;;
    *_dragonflye.fasta)        echo "${name%_dragonflye.fasta}" ;;
    *_EDGE.fasta)              echo "${name%_EDGE.fasta}" ;;
    *.fasta)                   echo "${name%.fasta}" ;;
    *.fa)                      echo "${name%.fa}" ;;
    *.fna)                     echo "${name%.fna}" ;;
    *)                         echo "$name" ;;
  esac
}

assembler_label_from_name() {
  local name="$1"
  case "$name" in
    *_dragonflye_contigs.fasta) echo "dragonflye_contigs" ;;
    *_dragonflye.fasta)        echo "dragonflye" ;;
    *_EDGE.fasta)              echo "EDGE" ;;
    *.fasta)                   echo "fasta" ;;
    *.fa)                      echo "fa" ;;
    *.fna)                     echo "fna" ;;
    *)                         echo "assembly" ;;
  esac
}

sample_from_stem() {
  local stem="$1"
  echo "$stem" | cut -d'_' -f1-5
}

extract_job_id() {
  # Accept common sbatch output forms:
  #   Submitted batch job 51729
  #   51729
  #   51729;cluster
  local text="$1"
  local id=""
  id="$(printf '%s\n' "$text" | awk '
    /Submitted batch job/ {print $NF}
    /^[0-9]+(;.*)?$/ {sub(/;.*/, "", $1); print $1}
  ' | tail -n 1)"
  if [[ "$id" =~ ^[0-9]+$ ]]; then
    printf '%s\n' "$id"
    return 0
  fi
  return 1
}

find_reads_for_stem() {
  local stem="$1"
  local candidate="$READS_DIR/${stem}${READ_EXT}"

  if [[ -f "$candidate" ]]; then
    echo "$candidate"
    return 0
  fi

  local fallbacks=(
    "$READS_DIR/${stem}.fastq.gz"
    "$READS_DIR/${stem}.fq.gz"
    "$READS_DIR/${stem}.fastq"
    "$READS_DIR/${stem}.fq"
  )

  local f
  for f in "${fallbacks[@]}"; do
    if [[ -f "$f" ]]; then
      echo "$f"
      return 0
    fi
  done

  shopt -s nullglob
  local hits=( "$READS_DIR/${stem}"*.fastq.gz "$READS_DIR/${stem}"*.fq.gz "$READS_DIR/${stem}"*.fastq "$READS_DIR/${stem}"*.fq )
  shopt -u nullglob

  if [[ ${#hits[@]} -eq 1 ]]; then
    echo "${hits[0]}"
    return 0
  elif [[ ${#hits[@]} -gt 1 ]]; then
    printf 'ERROR_MULTIPLE'
    printf '\t%s' "${hits[@]}"
    printf '\n'
    return 1
  fi

  return 1
}

# Collect assemblies safely.
IFS=',' read -r -a pattern_array <<< "$PATTERNS"
assemblies=()
shopt -s nullglob
for pattern in "${pattern_array[@]}"; do
  pattern="${pattern#"${pattern%%[![:space:]]*}"}"
  pattern="${pattern%"${pattern##*[![:space:]]}"}"
  for f in "$ASSEMBLY_DIR"/$pattern; do
    [[ -f "$f" ]] && assemblies+=( "$f" )
  done
done
shopt -u nullglob

# De-duplicate assemblies while preserving order.
deduped=()
declare -A seen=()
for f in "${assemblies[@]}"; do
  key="$(realpath "$f")"
  if [[ -z "${seen[$key]:-}" ]]; then
    deduped+=( "$key" )
    seen[$key]=1
  fi
done
assemblies=( "${deduped[@]}" )

if [[ ${#assemblies[@]} -eq 0 ]]; then
  echo "ERROR: no assemblies matched in $ASSEMBLY_DIR using patterns: $PATTERNS" >&2
  exit 1
fi

echo "ge_screen batch submission"
echo "  assembly-dir:          $ASSEMBLY_DIR"
echo "  reads-dir:             $READS_DIR"
echo "  outdir:                $OUTDIR"
echo "  ge-screen-script:      $GE_SCREEN_SCRIPT"
echo "  patterns:              $PATTERNS"
echo "  assemblies found:      ${#assemblies[@]}"
echo "  partition:             $PARTITION"
echo "  threads:               $THREADS"
echo "  mem:                   $MEM"
echo "  time:                  ${TIME_LIMIT:-not set}"
echo "  minimap2 preset:       $MINIMAP2_PRESET"
echo "  amrfinderplus:         $RUN_AMRFINDER"
echo "  amrfinder-db:          ${AMRFINDER_DB:-not set}"
echo "  amrfinder-organism:    ${AMRFINDER_ORGANISM:-not set}"
echo "  amrfinder-annotation-format: $AMRFINDER_ANNOTATION_FORMAT"
echo "  dry-run:               $DRY_RUN"
echo "  submit-dashboard:      $SUBMIT_DASHBOARD"
echo "  dashboard-dependency:  $DASHBOARD_DEPENDENCY"
echo "  dashboard-embed-plotly: $DASHBOARD_EMBED_PLOTLY"
echo "  dashboard-max-nt-rows: $DASHBOARD_MAX_NT_ROWS"
echo "  dashboard-max-univec-rows: $DASHBOARD_MAX_UNIVEC_ROWS"
echo "  dashboard-max-coverage-rows: $DASHBOARD_MAX_COVERAGE_ROWS"
echo "  dashboard-max-abricate-rows: $DASHBOARD_MAX_ABRICATE_ROWS"
echo "  dashboard-max-amrfinder-rows: $DASHBOARD_MAX_AMRFINDER_ROWS"
echo

submitted=0
skipped=0
failed=0
job_ids=()

summary_tsv="$OUTDIR/ge_screen_submission_summary.tsv"
printf "assembly\tstem\tsample\tassembler\treads\toutdir\tstatus\tjob_id\n" > "$summary_tsv"

for assembly in "${assemblies[@]}"; do
  name="$(basename "$assembly")"
  stem="$(assembly_stem_from_name "$name")"
  assembler="$(assembler_label_from_name "$name")"
  sample="$(sample_from_stem "$stem")"

  # Separate directory per assembly, even when multiple assemblies share a sample.
  sample_outdir="$OUTDIR/$sample/$stem"

  reads=""
  reads_result=""
  if reads_result="$(find_reads_for_stem "$stem" 2>/dev/null)"; then
    reads="$reads_result"
  else
    if [[ "$reads_result" == ERROR_MULTIPLE* ]]; then
      echo "SKIP: multiple possible read files for $name:"
      IFS=$'\t' read -r _tag -a multi_hits <<< "$reads_result"
      for hit in "${multi_hits[@]}"; do
        [[ -n "$hit" ]] && echo "  $hit"
      done
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$assembly" "$stem" "$sample" "$assembler" "" "$sample_outdir" "MULTIPLE_READ_MATCHES" "" >> "$summary_tsv"
    else
      echo "SKIP: no reads found for $name using stem '$stem' in $READS_DIR"
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$assembly" "$stem" "$sample" "$assembler" "" "$sample_outdir" "MISSING_READS" "" >> "$summary_tsv"
    fi
    skipped=$((skipped + 1))
    continue
  fi

  mkdir -p "$sample_outdir"

  cmd=(
    bash "$GE_SCREEN_SCRIPT"
    --assembly "$assembly"
    --reads "$reads"
    --outdir "$sample_outdir"
    --partition "$PARTITION"
    --threads "$THREADS"
    --mem "$MEM"
    --minimap2-preset "$MINIMAP2_PRESET"
  )

  if [[ "$RUN_AMRFINDER" == "true" ]]; then
    cmd+=( --amrfinderplus )
  else
    cmd+=( --no-amrfinderplus )
  fi

  if [[ -n "$AMRFINDER_DB" ]]; then
    cmd+=( --amrfinder-db "$AMRFINDER_DB" )
  fi

  if [[ -n "$AMRFINDER_ORGANISM" ]]; then
    cmd+=( --amrfinder-organism "$AMRFINDER_ORGANISM" )
  fi

  if [[ -n "$AMRFINDER_ANNOTATION_FORMAT" ]]; then
    cmd+=( --amrfinder-annotation-format "$AMRFINDER_ANNOTATION_FORMAT" )
  fi

  if [[ -n "$TIME_LIMIT" ]]; then
    cmd+=( --time "$TIME_LIMIT" )
  fi

  echo "Submitting:"
  echo "  assembly:  $assembly"
  echo "  assembler: $assembler"
  echo "  reads:     $reads"
  echo "  sample:    $sample"
  echo "  outdir:    $sample_outdir"

  if [[ "$DRY_RUN" == "true" ]]; then
    printf '  DRY-RUN command:'
    printf ' %q' "${cmd[@]}"
    printf '\n'
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$assembly" "$stem" "$sample" "$assembler" "$reads" "$sample_outdir" "DRY_RUN" "" >> "$summary_tsv"
  else
    set +e
    submit_output="$("${cmd[@]}" 2>&1)"
    submit_status=$?
    set -e

    printf '%s\n' "$submit_output"

    if [[ $submit_status -eq 0 ]]; then
      submitted=$((submitted + 1))
      job_id=""
      if job_id="$(extract_job_id "$submit_output")"; then
        job_ids+=( "$job_id" )
        echo "  captured job id: $job_id"
      else
        echo "  WARNING: submission succeeded but no SLURM job ID could be parsed from output." >&2
      fi
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$assembly" "$stem" "$sample" "$assembler" "$reads" "$sample_outdir" "SUBMITTED" "$job_id" >> "$summary_tsv"
    else
      failed=$((failed + 1))
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$assembly" "$stem" "$sample" "$assembler" "$reads" "$sample_outdir" "FAILED_TO_SUBMIT" "" >> "$summary_tsv"
    fi
  fi
  echo
done

echo "Submission summary:"
echo "  submitted: $submitted"
echo "  skipped:   $skipped"
echo "  failed:    $failed"
echo "  job IDs:   ${job_ids[*]:-none}"
echo "  summary:   $summary_tsv"

if [[ "$SUBMIT_DASHBOARD" == "true" ]]; then
  if [[ "$DRY_RUN" == "true" ]]; then
    echo
    echo "DRY-RUN: dashboard dependency would be based on job IDs: ${job_ids[*]:-none}"
  fi

  if [[ "$DASHBOARD_DEPENDENCY" != "none" && "$DASHBOARD_DEPENDENCY" != "singleton" && ${#job_ids[@]} -eq 0 ]]; then
    echo "WARNING: --submit-dashboard requested, but no processing job IDs were captured." >&2
    echo "         Dashboard will not be submitted with job dependency." >&2
    echo "         Use --dashboard-dependency none to submit it immediately, or check ge_screen.sh submission output." >&2
    exit 0
  fi

  if [[ -z "$DASHBOARD_WRAPPER" ]]; then
    if [[ -f "$SCRIPT_DIR/submit_ge_screen_dashboard.sh" ]]; then
      DASHBOARD_WRAPPER="$SCRIPT_DIR/submit_ge_screen_dashboard.sh"
    elif [[ -f "$SCRIPT_DIR/submit_ge_screen_dashboard_compact.sh" ]]; then
      DASHBOARD_WRAPPER="$SCRIPT_DIR/submit_ge_screen_dashboard_compact.sh"
    elif [[ -f "$PWD/submit_ge_screen_dashboard.sh" ]]; then
      DASHBOARD_WRAPPER="$PWD/submit_ge_screen_dashboard.sh"
    elif [[ -f "$PWD/submit_ge_screen_dashboard_compact.sh" ]]; then
      DASHBOARD_WRAPPER="$PWD/submit_ge_screen_dashboard_compact.sh"
    else
      echo "WARNING: --submit-dashboard requested but no dashboard submission wrapper was found." >&2
      echo "         Use --dashboard-wrapper /path/to/submit_ge_screen_dashboard.sh" >&2
      exit 0
    fi
  fi

  DASHBOARD_WRAPPER="$(realpath "$DASHBOARD_WRAPPER")"

  if [[ -z "$DASHBOARD_OUTPUT" ]]; then
    DASHBOARD_OUTPUT="$OUTDIR/ge_screen_dashboard.html"
  fi

  if [[ -z "$DASHBOARD_SUMMARY" ]]; then
    DASHBOARD_SUMMARY="$OUTDIR/ge_screen_dashboard_summary.tsv"
  fi

  dash_cmd=(
    bash "$DASHBOARD_WRAPPER"
    --input-dir "$OUTDIR"
    --output-html "$DASHBOARD_OUTPUT"
    --summary-tsv "$DASHBOARD_SUMMARY"
    --threads 2
    --mem 8G
    --max-nt-rows "$DASHBOARD_MAX_NT_ROWS"
    --max-univec-rows "$DASHBOARD_MAX_UNIVEC_ROWS"
    --max-coverage-rows "$DASHBOARD_MAX_COVERAGE_ROWS"
    --max-abricate-rows "$DASHBOARD_MAX_ABRICATE_ROWS"
    --max-amrfinder-rows "$DASHBOARD_MAX_AMRFINDER_ROWS"
  )

  if [[ -n "$PARTITION" ]]; then
    dash_cmd+=( --partition "$PARTITION" )
  fi

  if [[ -n "$DASHBOARD_SCRIPT" ]]; then
    dash_cmd+=( --dashboard-script "$DASHBOARD_SCRIPT" )
  fi

  if [[ "$DASHBOARD_EMBED_PLOTLY" == "true" ]]; then
    dash_cmd+=( --embed-plotly )
  fi

  # Pass dependency through environment to the dashboard wrapper if it supports EXTRA_SBATCH_ARGS.
  dependency_arg=""
  case "$DASHBOARD_DEPENDENCY" in
    afterok|afterany)
      IFS=: eval 'joined_job_ids="${job_ids[*]}"'
      dependency_arg="--dependency=${DASHBOARD_DEPENDENCY}:${joined_job_ids}"
      ;;
    singleton)
      dependency_arg="--dependency=singleton"
      ;;
    none)
      dependency_arg=""
      ;;
  esac

  echo
  echo "Submitting dashboard generation:"
  printf '  '
  if [[ -n "$dependency_arg" ]]; then
    printf 'EXTRA_SBATCH_ARGS=%q ' "$dependency_arg"
  fi
  printf ' %q' "${dash_cmd[@]}"
  printf '\n'

  if [[ "$DRY_RUN" == "true" ]]; then
    echo "  DRY-RUN: dashboard not submitted."
  else
    if [[ -n "$dependency_arg" ]]; then
      EXTRA_SBATCH_ARGS="$dependency_arg" "${dash_cmd[@]}"
    else
      "${dash_cmd[@]}"
    fi
  fi
fi
