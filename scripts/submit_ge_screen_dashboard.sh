#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Submit ge_screen dashboard generation to SLURM.

This wrapper is designed for ge_screen_dashboard_compact.py and supports
the compact-dashboard row-limit arguments.

Defaults:
  --input-dir          current directory
  --output-html        <input-dir>/ge_screen_dashboard.html
  --summary-tsv        <input-dir>/ge_screen_dashboard_summary.tsv
  --dashboard-script   ge_screen_dashboard.py in current directory,
                       otherwise ge_screen_dashboard.py next to this wrapper
  --title              Genetic Engineering Screen Dashboard
  --partition          not set; uses cluster default partition
  --threads            2
  --mem                8G
  --time               01:00:00
  --job-name           ge_screen_dashboard
  --dependency         none
  --max-nt-rows        100
  --max-univec-rows    100
  --max-coverage-rows  250
  --max-abricate-rows  500
  --max-amrfinder-rows 500
  --embed-plotly       false

Options:
  --input-dir DIR
  --output-html FILE
  --summary-tsv FILE
  --no-summary-tsv
  --dashboard-script FILE
  --title TEXT

  SLURM:
    --partition NAME
    --threads N
    --mem MEM
    --time HH:MM:SS
    --job-name NAME
    --dependency DEP       e.g. afterok:12345:12346 or afterany:12345:12346

  Compact dashboard:
    --max-nt-rows N        default 100; use -1 for all rows, not recommended
    --max-univec-rows N    default 100; use -1 for all rows
    --max-coverage-rows N  default 250; use -1 for all rows
    --max-abricate-rows N  default 500; use -1 for all rows
    --max-amrfinder-rows N  default 500; use -1 for all rows
    --embed-plotly         embed Plotly JS for offline viewing; adds several MB

Environment:
  EXTRA_SBATCH_ARGS        Extra sbatch argument string, e.g.
                           "--dependency=afterok:123:124"

Examples:
  # Submit using current directory and compact defaults
  bash submit_ge_screen_dashboard.sh

  # Submit for a specific result root
  bash submit_ge_screen_dashboard.sh \
    --input-dir /mnt/genomics/common_projects/COMMEX_Exercise/25Apr2026/ge_screen/ge_screen_results \
    --partition ss

  # Submit after processing jobs finish successfully
  bash submit_ge_screen_dashboard.sh \
    --input-dir /path/to/ge_screen_results \
    --partition ss \
    --dependency afterok:12345:12346

  # Offline dashboard with embedded Plotly
  bash submit_ge_screen_dashboard.sh \
    --input-dir /path/to/ge_screen_results \
    --embed-plotly

  # More compact dashboard
  bash submit_ge_screen_dashboard.sh \
    --input-dir /path/to/ge_screen_results \
    --max-nt-rows 50 \
    --max-univec-rows 50 \
    --max-coverage-rows 100
EOF
}

INPUT_DIR="."
OUTPUT_HTML=""
SUMMARY_TSV=""
NO_SUMMARY_TSV="false"
DASHBOARD_SCRIPT=""
TITLE="Genetic Engineering Screen Dashboard"

PARTITION=""
THREADS="2"
MEM="8G"
TIME_LIMIT="01:00:00"
JOB_NAME="ge_screen_dashboard"
DEPENDENCY=""

MAX_NT_ROWS="100"
MAX_UNIVEC_ROWS="100"
MAX_COVERAGE_ROWS="250"
MAX_ABRICATE_ROWS="500"
MAX_AMRFINDER_ROWS="500"
EMBED_PLOTLY="false"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input-dir) INPUT_DIR="${2:?}"; shift 2 ;;
    --output-html) OUTPUT_HTML="${2:?}"; shift 2 ;;
    --summary-tsv) SUMMARY_TSV="${2:?}"; shift 2 ;;
    --no-summary-tsv) NO_SUMMARY_TSV="true"; shift ;;
    --dashboard-script) DASHBOARD_SCRIPT="${2:?}"; shift 2 ;;
    --title) TITLE="${2:?}"; shift 2 ;;

    --partition) PARTITION="${2:?}"; shift 2 ;;
    --threads) THREADS="${2:?}"; shift 2 ;;
    --mem) MEM="${2:?}"; shift 2 ;;
    --time) TIME_LIMIT="${2:?}"; shift 2 ;;
    --job-name) JOB_NAME="${2:?}"; shift 2 ;;
    --dependency) DEPENDENCY="${2:?}"; shift 2 ;;

    --max-nt-rows) MAX_NT_ROWS="${2:?}"; shift 2 ;;
    --max-univec-rows) MAX_UNIVEC_ROWS="${2:?}"; shift 2 ;;
    --max-coverage-rows) MAX_COVERAGE_ROWS="${2:?}"; shift 2 ;;
    --max-abricate-rows) MAX_ABRICATE_ROWS="${2:?}"; shift 2 ;;
    --max-amrfinder-rows) MAX_AMRFINDER_ROWS="${2:?}"; shift 2 ;;
    --embed-plotly) EMBED_PLOTLY="true"; shift ;;

    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

is_int() {
  [[ "$1" =~ ^-?[0-9]+$ ]]
}

for pair in \
  "max-nt-rows:$MAX_NT_ROWS" \
  "max-univec-rows:$MAX_UNIVEC_ROWS" \
  "max-coverage-rows:$MAX_COVERAGE_ROWS" \
  "max-abricate-rows:$MAX_ABRICATE_ROWS" \
  "max-amrfinder-rows:$MAX_AMRFINDER_ROWS"
do
  name="${pair%%:*}"
  value="${pair#*:}"
  if ! is_int "$value"; then
    echo "ERROR: --$name must be an integer. Use -1 for all rows." >&2
    exit 2
  fi
done

INPUT_DIR="$(realpath "$INPUT_DIR")"

if [[ -z "$OUTPUT_HTML" ]]; then
  OUTPUT_HTML="$INPUT_DIR/ge_screen_dashboard.html"
fi

if [[ "$NO_SUMMARY_TSV" == "false" && -z "$SUMMARY_TSV" ]]; then
  SUMMARY_TSV="$INPUT_DIR/ge_screen_dashboard_summary.tsv"
fi

if [[ -z "$DASHBOARD_SCRIPT" ]]; then
  if [[ -f "ge_screen_dashboard.py" ]]; then
    DASHBOARD_SCRIPT="$PWD/ge_screen_dashboard.py"
  else
    WRAPPER_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    DASHBOARD_SCRIPT="$WRAPPER_DIR/ge_screen_dashboard.py"
  fi
fi

DASHBOARD_SCRIPT="$(realpath "$DASHBOARD_SCRIPT")"

if [[ ! -d "$INPUT_DIR" ]]; then
  echo "ERROR: input directory not found: $INPUT_DIR" >&2
  exit 2
fi

if [[ ! -f "$DASHBOARD_SCRIPT" ]]; then
  echo "ERROR: dashboard script not found: $DASHBOARD_SCRIPT" >&2
  echo "       Expected ge_screen_dashboard.py, or pass --dashboard-script /path/to/ge_screen_dashboard.py" >&2
  exit 2
fi

mkdir -p "$(dirname "$OUTPUT_HTML")"
LOG_DIR="$(dirname "$OUTPUT_HTML")/dashboard_slurm_logs"
mkdir -p "$LOG_DIR"

SBATCH_ARGS=(
  --parsable
  --job-name "$JOB_NAME"
  --cpus-per-task "$THREADS"
  --mem "$MEM"
  --time "$TIME_LIMIT"
  --output "$LOG_DIR/%x.%j.out"
  --error "$LOG_DIR/%x.%j.err"
)

if [[ -n "$PARTITION" ]]; then
  SBATCH_ARGS+=(--partition "$PARTITION")
fi

if [[ -n "$DEPENDENCY" ]]; then
  SBATCH_ARGS+=(--dependency "$DEPENDENCY")
fi

if [[ -n "${EXTRA_SBATCH_ARGS:-}" ]]; then
  # Accept a simple extra sbatch argument string, e.g.
  # EXTRA_SBATCH_ARGS="--dependency=afterok:123:124"
  # shellcheck disable=SC2206
  EXTRA_ARGS=( ${EXTRA_SBATCH_ARGS} )
  SBATCH_ARGS+=( "${EXTRA_ARGS[@]}" )
fi

CMD=(
  python "$DASHBOARD_SCRIPT"
  --input-dir "$INPUT_DIR"
  --output-html "$OUTPUT_HTML"
  --title "$TITLE"
  --max-nt-rows "$MAX_NT_ROWS"
  --max-univec-rows "$MAX_UNIVEC_ROWS"
  --max-coverage-rows "$MAX_COVERAGE_ROWS"
  --max-abricate-rows "$MAX_ABRICATE_ROWS"
  --max-amrfinder-rows "$MAX_AMRFINDER_ROWS"
)

if [[ "$NO_SUMMARY_TSV" == "true" ]]; then
  CMD+=(--no-summary-tsv)
elif [[ -n "$SUMMARY_TSV" ]]; then
  CMD+=(--summary-tsv "$SUMMARY_TSV")
fi

if [[ "$EMBED_PLOTLY" == "true" ]]; then
  CMD+=(--embed-plotly)
fi

echo "Submitting ge_screen dashboard job with:"
echo "  input-dir:          $INPUT_DIR"
echo "  output-html:        $OUTPUT_HTML"
if [[ "$NO_SUMMARY_TSV" == "true" ]]; then
  echo "  summary-tsv:        disabled"
else
  echo "  summary-tsv:        $SUMMARY_TSV"
fi
echo "  dashboard-script:   $DASHBOARD_SCRIPT"
echo "  title:              $TITLE"
echo "  partition:          ${PARTITION:-cluster default}"
echo "  threads:            $THREADS"
echo "  mem:                $MEM"
echo "  time:               $TIME_LIMIT"
echo "  dependency:         ${DEPENDENCY:-${EXTRA_SBATCH_ARGS:-none}}"
echo "  max-nt-rows:        $MAX_NT_ROWS"
echo "  max-univec-rows:    $MAX_UNIVEC_ROWS"
echo "  max-coverage-rows:  $MAX_COVERAGE_ROWS"
echo "  max-abricate-rows:  $MAX_ABRICATE_ROWS"
echo "  max-amrfinder-rows: $MAX_AMRFINDER_ROWS"
echo "  embed-plotly:       $EMBED_PLOTLY"
echo "  logs:               $LOG_DIR"

jobid="$(
  sbatch "${SBATCH_ARGS[@]}" --wrap "$(printf '%q ' "${CMD[@]}")"
)"

echo "Submitted dashboard job: $jobid"
echo "Dashboard target: $OUTPUT_HTML"
echo "SLURM logs: $LOG_DIR"
