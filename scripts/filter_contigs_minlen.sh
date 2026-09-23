#!/usr/bin/env bash
set -euo pipefail

# Filter contigs in a FASTA file by minimum length (nucleotides).
# Keeps sequences with length >= MIN_LEN (default: 250).
#
# Usage:
#   ./filter_contigs_minlen.sh <input.fasta> <output.fasta>
#
# Examples:
#   ./filter_contigs_minlen.sh contigs.fasta contigs.min250.fasta
#   MIN_LEN=500 ./filter_contigs_minlen.sh contigs.fasta contigs.min500.fasta

IN="${1:-}"
OUT="${2:-}"

if [[ -z "${IN}" || -z "${OUT}" ]]; then
  echo "Usage: $0 <input.fasta> <output.fasta>" >&2
  exit 1
fi

if [[ ! -f "${IN}" ]]; then
  echo "ERROR: Input file not found: ${IN}" >&2
  exit 1
fi

MIN_LEN="${MIN_LEN:-250}"

python3 - "${IN}" "${OUT}" "${MIN_LEN}" <<'PY'
import sys

inp = sys.argv[1]
outp = sys.argv[2]
min_len = int(sys.argv[3])

def write_record(out_handle, header, seq, width=60):
    out_handle.write(header + "\n")
    for i in range(0, len(seq), width):
        out_handle.write(seq[i:i+width] + "\n")

total = kept = 0
header = None
seq_parts = []

with open(inp, "r") as f, open(outp, "w") as out:
    for line in f:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                seq = "".join(seq_parts)
                total += 1
                if len(seq) >= min_len:
                    write_record(out, header, seq)
                    kept += 1
            header = line
            seq_parts = []
        else:
            seq_parts.append(line)

    if header is not None:
        seq = "".join(seq_parts)
        total += 1
        if len(seq) >= min_len:
            write_record(out, header, seq)
            kept += 1

removed = total - kept
print(f"Total contigs: {total}")
print(f"Kept (>= {min_len} nt): {kept}")
print(f"Removed (< {min_len} nt): {removed}")
PY
