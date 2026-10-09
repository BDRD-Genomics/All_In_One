#!/usr/bin/env bash

set -u

REGISTRY="ghcr.io/bdrd-genomics"
TAG="latest"

IMAGES=(
    allinone-aio_qc
    allinone-amrfinder_plus
    allinone-autocycler
    allinone-bbmap
    allinone-blast
    allinone-busco
    allinone-checkm-genome
    allinone-diamond
    allinone-dragonflye
    allinone-fastp_long
    allinone-fastqc
    allinone-interleave
    allinone-kraken2
    allinone-krona
    allinone-md
    allinone-megan
    allinone-mmseqs2
    allinone-multiqc
    allinone-myloasm
    allinone-nanoplot
    allinone-nextflow
    allinone-phispy
    allinone-prokka
    allinone-python_reports
    allinone-python_utils
    allinone-readmapping_auto
    allinone-rgi
    allinone-ribodetector
    allinone-spades
    allinone-unicycler
    allinone-vs
)

FAILED=()

for IMAGE in "${IMAGES[@]}"; do
    FULL_IMAGE="${REGISTRY}/${IMAGE}:${TAG}"

    echo
    echo "=================================================="
    echo "Pulling: ${FULL_IMAGE}"
    echo "=================================================="

    if docker pull "${FULL_IMAGE}"; then
        echo "SUCCESS: ${FULL_IMAGE}"
    else
        echo "FAILED: ${FULL_IMAGE}"
        FAILED+=("${FULL_IMAGE}")
    fi
done

echo
echo "=================================================="
echo "SUMMARY"
echo "=================================================="

echo "Total images: ${#IMAGES[@]}"
echo "Failed:       ${#FAILED[@]}"

if (( ${#FAILED[@]} > 0 )); then
    echo
    echo "Failed images:"
    printf '  %s\n' "${FAILED[@]}"
    exit 1
else
    echo
    echo "All images pulled successfully."
fi
