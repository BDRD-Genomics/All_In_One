#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

usage() {
  cat <<'USAGE'
Usage:
  ge_screen_full_fixed.sh \
    --assembly 'assemblies/*.fasta' \
    [--reads 'reads/*.fastq.gz' | (--reads1 'reads/*_R1.fastq.gz' --reads2 'reads/*_R2.fastq.gz')] \
    [--outdir engineering_screen] \
    [--query-fasta query.fasta] \
    [--partition batch] [--threads 16] [--mem 32G] [--time 12:00:00] \
    [--minimap2-preset auto|sr|map-ont|map-pb|asm5|asm10|asm20] \
    [--amrfinderplus|--no-amrfinderplus] [--amrfinder-db DIR] [--amrfinder-organism NAME]

Required:
  --assembly         Assembly path or quoted wildcard pattern

Read input modes (choose one):
  --reads            Single read file path/pattern (ONT or single-end)
  --reads1           Paired-end R1 path/pattern
  --reads2           Paired-end R2 path/pattern

Optional:
  --outdir           Parent output directory [engineering_screen]
  --query-fasta      FASTA to BLAST against nt/UniVec [defaults to assembly]
  --partition        SLURM partition [batch]
  --threads          Threads / cpus-per-task [16]
  --mem              SLURM memory [32G]
  --time             SLURM walltime [12:00:00]
  --minimap2-preset  auto chooses map-ont for --reads and sr for --reads1/--reads2 [auto]
  --amrfinderplus     Run AMRFinderPlus with --plus after Prodigal gene calling [default]
  --no-amrfinderplus  Disable AMRFinderPlus
  --amrfinder-db      Optional AMRFinderPlus database directory passed with -d
  --amrfinder-organism Optional AMRFinderPlus organism name passed with -O
  --amrfinder-annotation-format AMRFinderPlus --annotation_format value [prodigal]
  -h, --help         Show this help

Behavior:
  - Submit mode: expands globs and submits one SLURM job per matched assembly.
  - Run mode: internal only, executed under sbatch.
USAGE
}

SCRIPT_PATH="$(readlink -f "$0" 2>/dev/null || realpath "$0")"

ASSEMBLY_PATTERN=""
READS_PATTERN=""
READS1_PATTERN=""
READS2_PATTERN=""
OUTDIR="engineering_screen"
QUERY_FASTA=""
PARTITION="batch"
THREADS="16"
MEM="32G"
TIME_LIMIT="12:00:00"
MINIMAP2_PRESET="auto"
RUN_AMRFINDER="true"
AMRFINDER_DB=""
AMRFINDER_ORGANISM=""
AMRFINDER_ANNOTATION_FORMAT="prodigal"
RUN_MODE=0

# -------- parse args --------
if [[ "${1:-}" == "__run" ]]; then
  RUN_MODE=1
  shift
fi

while [[ $# -gt 0 ]]; do
  case "$1" in
    --assembly) ASSEMBLY_PATTERN="$2"; shift 2 ;;
    --reads) READS_PATTERN="$2"; shift 2 ;;
    --reads1) READS1_PATTERN="$2"; shift 2 ;;
    --reads2) READS2_PATTERN="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --query-fasta) QUERY_FASTA="$2"; shift 2 ;;
    --partition) PARTITION="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --mem) MEM="$2"; shift 2 ;;
    --time) TIME_LIMIT="$2"; shift 2 ;;
    --minimap2-preset) MINIMAP2_PRESET="$2"; shift 2 ;;
    --amrfinderplus) RUN_AMRFINDER="true"; shift ;;
    --no-amrfinderplus) RUN_AMRFINDER="false"; shift ;;
    --amrfinder-db) AMRFINDER_DB="$2"; shift 2 ;;
    --amrfinder-organism) AMRFINDER_ORGANISM="$2"; shift 2 ;;
    --amrfinder-annotation-format) AMRFINDER_ANNOTATION_FORMAT="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *)
      echo "ERROR: Unknown argument: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
done

join_by_tab() {
  local IFS=$'\t'
  printf '%s' "$*"
}

resolve_glob() {
  local pattern="$1"
  local -n out_ref=$2
  out_ref=()
  local matches=( $pattern )
  if [[ ${#matches[@]} -eq 0 ]]; then
    return 1
  fi
  out_ref=( "${matches[@]}" )
  return 0
}

sanitize_stem() {
  local p="$1"
  local b
  b="$(basename "$p")"
  b="${b%.fasta}"
  b="${b%.fa}"
  b="${b%.fna}"
  b="${b%.fas}"
  printf '%s' "$b"
}

# -------- run mode --------
if [[ "$RUN_MODE" -eq 1 ]]; then
  # Values arrive via --export from sbatch
  : "${JOB_ASSEMBLY:?JOB_ASSEMBLY not set}"
  : "${JOB_OUTDIR:?JOB_OUTDIR not set}"
  : "${JOB_READ_MODE:?JOB_READ_MODE not set}"
  : "${JOB_MINIMAP2_PRESET:?JOB_MINIMAP2_PRESET not set}"

  ASSEMBLY="$JOB_ASSEMBLY"
  SAMPLE_OUTDIR="$JOB_OUTDIR"
  QUERY_FASTA_RUN="${JOB_QUERY_FASTA:-$ASSEMBLY}"
  READ_MODE="$JOB_READ_MODE"
  PRESET="$JOB_MINIMAP2_PRESET"
  RUN_AMRFINDER_RUN="${JOB_RUN_AMRFINDER:-true}"
  AMRFINDER_DB_RUN="${JOB_AMRFINDER_DB:-}"
  AMRFINDER_ORGANISM_RUN="${JOB_AMRFINDER_ORGANISM:-}"
  AMRFINDER_ANNOTATION_FORMAT_RUN="${JOB_AMRFINDER_ANNOTATION_FORMAT:-prodigal}"

  mkdir -p "$SAMPLE_OUTDIR"

  unset PERL5LIB PERL_LOCAL_LIB_ROOT PERL_MB_OPT PERL_MM_OPT

  assembly_stem="$(sanitize_stem "$ASSEMBLY")"
  query_stem="$(sanitize_stem "$QUERY_FASTA_RUN")"

  echo "Host: $(hostname)"
  echo "Assembly: $ASSEMBLY"
  echo "Query FASTA: $QUERY_FASTA_RUN"
  echo "Output dir: $SAMPLE_OUTDIR"
  echo "Read mode: $READ_MODE"
  echo "minimap2 preset: $PRESET"
  echo "AMRFinderPlus: $RUN_AMRFINDER_RUN"
  if [[ -n "$AMRFINDER_DB_RUN" ]]; then echo "AMRFinderPlus DB: $AMRFINDER_DB_RUN"; fi
  if [[ -n "$AMRFINDER_ORGANISM_RUN" ]]; then echo "AMRFinderPlus organism: $AMRFINDER_ORGANISM_RUN"; fi
  echo "AMRFinderPlus annotation format: $AMRFINDER_ANNOTATION_FORMAT_RUN"
  echo "Threads: ${SLURM_CPUS_PER_TASK:-$THREADS}"

  echo "[1/6] BLAST against NCBI core_nt"
  blastn \
    -query "$QUERY_FASTA_RUN" \
    -db /export/database/blastdb/core_nt/core_nt \
    -task megablast \
    -evalue 1e-20 \
    -max_target_seqs 20 \
    -num_threads "${SLURM_CPUS_PER_TASK:-$THREADS}" \
    -outfmt '6 qseqid qlen qstart qend sacc stitle sstart send length slen pident qcovs sstrand gaps evalue bitscore score staxids' \
    > "$SAMPLE_OUTDIR/${query_stem}_vs_nt.tsv"

  echo "[2/6] Screen against UniVec"
  blastn \
    -query "$QUERY_FASTA_RUN" \
    -db /export/database/blastdb/UniVec_DB/UniVec_Core \
    -task blastn \
    -evalue 1e-10 \
    -max_target_seqs 50 \
    -num_threads "${SLURM_CPUS_PER_TASK:-$THREADS}" \
    -outfmt '6 qseqid qlen qstart qend sacc stitle sstart send length slen pident qcovs sstrand gaps evalue bitscore score' \
    > "$SAMPLE_OUTDIR/${query_stem}_vs_univec.tsv" || true

  echo "[3/6] Plasmid / AMR / VF screens with abricate"
  abricate --db plasmidfinder "$QUERY_FASTA_RUN" > "$SAMPLE_OUTDIR/abricate.plasmidfinder.tsv" || true
  abricate --db card          "$QUERY_FASTA_RUN" > "$SAMPLE_OUTDIR/abricate.card.tsv" || true
  abricate --db vfdb          "$QUERY_FASTA_RUN" > "$SAMPLE_OUTDIR/abricate.vfdb.tsv" || true
  abricate --summary "$SAMPLE_OUTDIR"/abricate.*.tsv > "$SAMPLE_OUTDIR/abricate.summary.tsv" || true

  echo "[4/6] Map reads to assembly for depth"
  bam="$SAMPLE_OUTDIR/${assembly_stem}.readmap.bam"
  case "$READ_MODE" in
    single)
      IFS=$'\t' read -r -a READS_ARR <<< "${JOB_READS_STR:-}"
      if [[ ${#READS_ARR[@]} -eq 0 || -z "${READS_ARR[0]:-}" ]]; then
        echo "ERROR: No reads supplied in single-read mode" >&2
        exit 1
      fi
      minimap2 -ax "$PRESET" -t "${SLURM_CPUS_PER_TASK:-$THREADS}" "$ASSEMBLY" "${READS_ARR[@]}" \
        | samtools sort -@ "${SLURM_CPUS_PER_TASK:-$THREADS}" -o "$bam"
      ;;
    paired)
      IFS=$'\t' read -r -a READS1_ARR <<< "${JOB_READS1_STR:-}"
      IFS=$'\t' read -r -a READS2_ARR <<< "${JOB_READS2_STR:-}"
      if [[ ${#READS1_ARR[@]} -eq 0 || ${#READS2_ARR[@]} -eq 0 ]]; then
        echo "ERROR: No reads supplied in paired mode" >&2
        exit 1
      fi
      if [[ ${#READS1_ARR[@]} -ne ${#READS2_ARR[@]} ]]; then
        echo "ERROR: Paired read counts differ in run mode (${#READS1_ARR[@]} vs ${#READS2_ARR[@]})" >&2
        exit 1
      fi
      map_inputs=()
      for idx in "${!READS1_ARR[@]}"; do
        map_inputs+=( "${READS1_ARR[$idx]}" "${READS2_ARR[$idx]}" )
      done
      minimap2 -ax "$PRESET" -t "${SLURM_CPUS_PER_TASK:-$THREADS}" "$ASSEMBLY" "${map_inputs[@]}" \
        | samtools sort -@ "${SLURM_CPUS_PER_TASK:-$THREADS}" -o "$bam"
      ;;
    *)
      echo "ERROR: Unsupported JOB_READ_MODE=$READ_MODE" >&2
      exit 1
      ;;
  esac
  samtools index "$bam"
  samtools coverage "$bam" > "$SAMPLE_OUTDIR/${assembly_stem}.coverage.tsv"
  samtools depth -aa "$bam" > "$SAMPLE_OUTDIR/${assembly_stem}.depth.tsv"

  echo "[5/7] ORF calling with Prodigal"
  proteins_faa="$SAMPLE_OUTDIR/${query_stem}.proteins.faa"
  genes_fna="$SAMPLE_OUTDIR/${query_stem}.genes.fna"
  prodigal_gff="$SAMPLE_OUTDIR/${query_stem}.prodigal.gff"

  prodigal -p meta -i "$QUERY_FASTA_RUN" \
    -a "$proteins_faa" \
    -d "$genes_fna" \
    -o "$prodigal_gff" \
    -f gff || true

  echo "[6/7] AMRFinderPlus AMR/stress/virulence screen"
  amrfinder_out="$SAMPLE_OUTDIR/${query_stem}.amrfinderplus.tsv"
  if [[ "$RUN_AMRFINDER_RUN" == "true" ]]; then
    if command -v amrfinder >/dev/null 2>&1; then
      if [[ -s "$proteins_faa" && -s "$prodigal_gff" ]]; then
        amrfinder_cmd=(
          amrfinder
          -p "$proteins_faa"
          -n "$QUERY_FASTA_RUN"
          -g "$prodigal_gff"
          --plus
          --annotation_format "$AMRFINDER_ANNOTATION_FORMAT_RUN"
          --threads "${SLURM_CPUS_PER_TASK:-$THREADS}"
          -o "$amrfinder_out"
        )

        if [[ -n "$AMRFINDER_DB_RUN" ]]; then
          amrfinder_cmd+=( -d "$AMRFINDER_DB_RUN" )
        fi

        if [[ -n "$AMRFINDER_ORGANISM_RUN" ]]; then
          amrfinder_cmd+=( -O "$AMRFINDER_ORGANISM_RUN" )
        fi

        "${amrfinder_cmd[@]}" || {
          echo "WARNING: AMRFinderPlus failed for $QUERY_FASTA_RUN" >&2
          : > "$amrfinder_out"
        }
      else
        echo "WARNING: Prodigal protein/GFF output missing or empty; skipping AMRFinderPlus" >&2
        : > "$amrfinder_out"
      fi
    else
      echo "WARNING: amrfinder executable not found; skipping AMRFinderPlus" >&2
      : > "$amrfinder_out"
    fi
  else
    echo "AMRFinderPlus disabled by --no-amrfinderplus"
    : > "$amrfinder_out"
  fi

  echo "[7/7] Done. Review:"
  echo "  $SAMPLE_OUTDIR/${query_stem}_vs_nt.tsv"
  echo "  $SAMPLE_OUTDIR/${query_stem}_vs_univec.tsv"
  echo "  $SAMPLE_OUTDIR/abricate.summary.tsv"
  echo "  $SAMPLE_OUTDIR/${assembly_stem}.coverage.tsv"
  echo "  $proteins_faa"
  echo "  $prodigal_gff"
  echo "  $amrfinder_out"
  exit 0
fi

# -------- submit mode --------
if [[ -z "$ASSEMBLY_PATTERN" ]]; then
  echo "ERROR: --assembly is required." >&2
  usage >&2
  exit 1
fi

if [[ -n "$READS_PATTERN" && ( -n "$READS1_PATTERN" || -n "$READS2_PATTERN" ) ]]; then
  echo "ERROR: Use either --reads or --reads1/--reads2, not both." >&2
  exit 1
fi

if [[ -z "$READS_PATTERN" && ( -z "$READS1_PATTERN" || -z "$READS2_PATTERN" ) ]]; then
  echo "ERROR: Provide --reads or both --reads1 and --reads2." >&2
  exit 1
fi

assemblies=()
if ! resolve_glob "$ASSEMBLY_PATTERN" assemblies; then
  echo "ERROR: No assemblies matched pattern: $ASSEMBLY_PATTERN" >&2
  exit 1
fi

READ_MODE="single"
reads=()
reads1=()
reads2=()
if [[ -n "$READS_PATTERN" ]]; then
  if ! resolve_glob "$READS_PATTERN" reads; then
    echo "ERROR: No reads matched pattern: $READS_PATTERN" >&2
    exit 1
  fi
  READ_MODE="single"
else
  if ! resolve_glob "$READS1_PATTERN" reads1; then
    echo "ERROR: No reads1 matched pattern: $READS1_PATTERN" >&2
    exit 1
  fi
  if ! resolve_glob "$READS2_PATTERN" reads2; then
    echo "ERROR: No reads2 matched pattern: $READS2_PATTERN" >&2
    exit 1
  fi
  if [[ ${#reads1[@]} -ne ${#reads2[@]} ]]; then
    echo "ERROR: --reads1 and --reads2 matched different counts (${#reads1[@]} vs ${#reads2[@]})." >&2
    exit 1
  fi
  READ_MODE="paired"
fi

case "$MINIMAP2_PRESET" in
  auto)
    if [[ "$READ_MODE" == "paired" ]]; then
      EFFECTIVE_PRESET="sr"
    else
      EFFECTIVE_PRESET="map-ont"
    fi
    ;;
  sr|map-ont|map-pb|asm5|asm10|asm20)
    EFFECTIVE_PRESET="$MINIMAP2_PRESET"
    ;;
  *)
    echo "ERROR: Unsupported --minimap2-preset value: $MINIMAP2_PRESET" >&2
    exit 1
    ;;
esac

mkdir -p "$OUTDIR"

echo "Assemblies matched (${#assemblies[@]}):"
printf '  %s\n' "${assemblies[@]}"
if [[ "$READ_MODE" == "single" ]]; then
  echo "Reads matched (${#reads[@]}):"
  printf '  %s\n' "${reads[@]}"
else
  echo "Reads1 matched (${#reads1[@]}):"
  printf '  %s\n' "${reads1[@]}"
  echo "Reads2 matched (${#reads2[@]}):"
  printf '  %s\n' "${reads2[@]}"
fi
echo "minimap2 preset: $MINIMAP2_PRESET -> $EFFECTIVE_PRESET"
echo "AMRFinderPlus: $RUN_AMRFINDER"
if [[ -n "$AMRFINDER_DB" ]]; then echo "AMRFinderPlus DB: $AMRFINDER_DB"; fi
if [[ -n "$AMRFINDER_ORGANISM" ]]; then echo "AMRFinderPlus organism: $AMRFINDER_ORGANISM"; fi
echo "AMRFinderPlus annotation format: $AMRFINDER_ANNOTATION_FORMAT"

READS_STR=""
READS1_STR=""
READS2_STR=""
if [[ "$READ_MODE" == "single" ]]; then
  READS_STR="$(join_by_tab "${reads[@]}")"
else
  READS1_STR="$(join_by_tab "${reads1[@]}")"
  READS2_STR="$(join_by_tab "${reads2[@]}")"
fi

for assembly in "${assemblies[@]}"; do
  if [[ ! -r "$assembly" ]]; then
    echo "WARNING: Skipping unreadable assembly: $assembly" >&2
    continue
  fi
  assembly_stem="$(sanitize_stem "$assembly")"
  sample_outdir="$OUTDIR/$assembly_stem"
  mkdir -p "$sample_outdir"
  query_fasta="$assembly"
  if [[ -n "$QUERY_FASTA" ]]; then
    query_fasta="$QUERY_FASTA"
  fi
  job_out="$sample_outdir/ge_screen.%j.out"
  job_err="$sample_outdir/ge_screen.%j.err"

  echo "Submitting SLURM job for assembly: $assembly"
  jobid=$(sbatch \
    --parsable \
    --job-name "ge_screen_${assembly_stem}" \
    --partition "$PARTITION" \
    --cpus-per-task "$THREADS" \
    --mem "$MEM" \
    --time "$TIME_LIMIT" \
    --output "$job_out" \
    --error "$job_err" \
    --export=ALL,JOB_ASSEMBLY="$assembly",JOB_OUTDIR="$sample_outdir",JOB_QUERY_FASTA="$query_fasta",JOB_READ_MODE="$READ_MODE",JOB_MINIMAP2_PRESET="$EFFECTIVE_PRESET",JOB_READS_STR="$READS_STR",JOB_READS1_STR="$READS1_STR",JOB_READS2_STR="$READS2_STR",JOB_RUN_AMRFINDER="$RUN_AMRFINDER",JOB_AMRFINDER_DB="$AMRFINDER_DB",JOB_AMRFINDER_ORGANISM="$AMRFINDER_ORGANISM",JOB_AMRFINDER_ANNOTATION_FORMAT="$AMRFINDER_ANNOTATION_FORMAT" \
    "$SCRIPT_PATH" __run)
  echo "Submitted batch job $jobid"
done
