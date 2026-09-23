#!/bin/bash
set -euo pipefail

bamfile="$1"

#source "$(conda info --base)/etc/profile.d/conda.sh"
#conda activate samtools

samtools split -f '%!.bam' "$bamfile"

# Strip everything before known ONT prefix
for f in *.bam; do
    newname=$(echo "$f" | sed 's/^.*_dna_r10\.4\.1_e8\.2_400bps_sup@v5\.0\.0_//')
    if [[ "$f" != "$newname" ]]; then
        mv "$f" "$newname"
    fi
done

