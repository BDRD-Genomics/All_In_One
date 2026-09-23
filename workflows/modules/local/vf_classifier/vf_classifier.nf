#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process VF_CLASSIFIER {
    tag { "${assembler} | ${sample_id}" }
    publishDir { "${params.outdir}/${params.project_id}/${sample_id}/vf_classifier/${assembler}" }, mode: 'copy'
    errorStrategy 'terminate'

    input:
    tuple val(assembler), val(sample_id), path(input_file)

    output:
    tuple val(sample_id), val(assembler), path("${sample_id}_annotated.csv"), emit: annotated
    tuple val(sample_id), val(assembler), path("${sample_id}_summary.txt"), emit: summary
    tuple val(sample_id), val(assembler), path("${sample_id}_vf_calls.tsv"), emit: model_calls
    tuple val(sample_id), val(assembler), path("${sample_id}_proteins.faa"), optional: true, emit: proteins
    tuple val(sample_id), val(assembler), path("${sample_id}_genes.gff"), optional: true, emit: genes
    tuple val(sample_id), val(assembler), path("${sample_id}_vfdb.tsv"), optional: true, emit: vfdb_calls
    tuple val(sample_id), val(assembler), path("${sample_id}_per_contig.tsv"), optional: true, emit: per_contig
    tuple val(sample_id), val(assembler), path("${sample_id}_package"), emit: package_dir
    tuple val(sample_id), val(assembler), path("${sample_id}_3-review-no-homology.csv"), optional: true, emit: review

    script:
    def vfroot = params.vf_classifier_root
    def pmode = params.vf_mode == 'metagenome' ? 'meta' : 'single'
    """
    set -euo pipefail

    PREFIX="${sample_id}"
    PROTEINS="${input_file}"

    if [[ "${params.vf_input_type}" == "nucleotide" ]]; then
        prodigal \
            -i "${input_file}" \
            -a "\${PREFIX}_proteins.faa" \
            -f gff \
            -o "\${PREFIX}_genes.gff" \
            -p "${pmode}" \
            -q
        PROTEINS="\${PREFIX}_proteins.faa"
    fi

    python ${vfroot}/scripts/inference.py \
        --config ${vfroot}/config.yml \
        --checkpoint ${vfroot}/checkpoint \
        --input "\${PROTEINS}" \
        --output "\${PREFIX}_vf_calls.tsv" \
        --threshold ${params.vf_threshold}

    VFDB_ARG=""
    if [[ "${params.vf_homology}" == "yes" ]]; then
        python ${vfroot}/scripts/build_truth.py \
            --query "\${PROTEINS}" \
            --vfdb "${params.vfdb}" \
            --out "\${PREFIX}_vfdb.tsv" \
            --threads ${task.cpus}
        VFDB_ARG="--vfdb \${PREFIX}_vfdb.tsv"
    fi

    BYCONTIG_ARG=""
    if [[ "${params.vf_mode}" == "metagenome" ]]; then
        BYCONTIG_ARG="--by-contig \${PREFIX}_per_contig.tsv"
    fi

    python ${vfroot}/scripts/interpret_calls.py \
        --model "\${PREFIX}_vf_calls.tsv" \
        \${VFDB_ARG} \
        --out "\${PREFIX}_annotated.csv" \
        --summary "\${PREFIX}_summary.txt" \
        \${BYCONTIG_ARG} \
        --threshold ${params.vf_threshold}

    python ${params.scripts}/extract_review.py \
        --annotated "\${PREFIX}_annotated.csv" \
        --out "\${PREFIX}_3-review-no-homology.csv" || true

    mkdir -p "\${PREFIX}_package"
    cp -f "\${PREFIX}_vf_calls.tsv" "\${PREFIX}_annotated.csv" "\${PREFIX}_summary.txt" "\${PREFIX}_package/"
    for extra in "\${PREFIX}_proteins.faa" "\${PREFIX}_genes.gff" "\${PREFIX}_vfdb.tsv" "\${PREFIX}_per_contig.tsv" "\${PREFIX}_3-review-no-homology.csv"; do
        [[ -f "\${extra}" ]] && cp -f "\${extra}" "\${PREFIX}_package/"
    done
    """
}
