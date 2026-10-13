#!/bin/bash

if [[ $# -ne 3 ]]; then
	echo "Usage: $0 <APPTAINER_TMPDIR> <APPTAINER_CACHEDIR> <SIF_DIR>" >&2
	exit 1
fi

export APPTAINER_TMPDIR=$1
export APPTAINER_CACHEDIR=$2
sif_dir=$3

mkdir -p "$APPTAINER_TMPDIR" "$APPTAINER_CACHEDIR" "$sif_dir"

while IFS= read -r image; do
    name="${image##*/}"
    name="${name%:latest}"
    sif="${sif_dir}/${name}.sif"

    if [[ -s "$sif" ]]; then
        echo "Skipping existing: $sif"
        continue
    fi

    echo "Building: $image -> $sif"
    if apptainer build "${sif}.partial" "docker-daemon://${image}"; then
        mv "${sif}.partial" "$sif"
    else
        echo "FAILED: $image" >&2
        rm -f "${sif}.partial"
    fi
done < <(
    docker image ls --format '{{.Repository}}:{{.Tag}}' |
        awk '/(^|\/)allinone-[^:]+:latest$/' |
        sort -u
)
