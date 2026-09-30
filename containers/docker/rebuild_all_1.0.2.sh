#!/usr/bin/env bash
set -u
set -o pipefail

REGISTRY="ghcr.io/bdrd-genomics"
TAG="1.0.2"
SOURCE_URL="https://github.com/BDRD-Genomics/All_In_One"

passed=()
build_failed=()
ps_failed=()
push_failed=()
skipped=()

BASE="${REGISTRY}/allinone-conda-base:latest"

echo "============================================================"
echo "All_In_One Docker rebuild + validation + push"
echo "Release: ${TAG}"
echo "Base:    ${BASE}"
echo "Source:  ${SOURCE_URL}"
echo "============================================================"
echo


# ------------------------------------------------------------
# Validate the base image first
# ------------------------------------------------------------

echo "Checking base image..."

if ! docker image inspect "$BASE" >/dev/null 2>&1; then
    echo "ERROR: Missing local base image:"
    echo "  ${BASE}"
    exit 1
fi

if ! docker run --rm \
    --entrypoint /bin/sh \
    "$BASE" \
    -c '
        command -v ps >/dev/null 2>&1 || exit 10
        ps -eo pid,ppid,comm >/dev/null 2>&1 || exit 11
    '
then
    echo "ERROR: Base image does not have working ps:"
    echo "  ${BASE}"
    exit 1
fi

echo "PASS: base image"
echo


# ------------------------------------------------------------
# Discover application Dockerfiles
# ------------------------------------------------------------

mapfile -t dockerfiles < <(
    find . \
        -maxdepth 1 \
        -type f \
        -name 'Dockerfile.*' \
        ! -name '*.bak' \
        -printf '%f\n' \
    | sort
)


# ------------------------------------------------------------
# Build → test → push
# ------------------------------------------------------------

for dockerfile in "${dockerfiles[@]}"; do

    name="${dockerfile#Dockerfile.}"

    case "$name" in
        base|conda-base|almalinux9.autoactivate|tool-template)
            echo "SKIP: ${dockerfile}"
            skipped+=("$dockerfile")
            continue
            ;;
    esac

    versioned="${REGISTRY}/allinone-${name}:${TAG}"
    latest="${REGISTRY}/allinone-${name}:latest"

    echo
    echo "============================================================"
    echo "Building:   ${versioned}"
    echo "Dockerfile: ${dockerfile}"
    echo "============================================================"

    if ! docker build \
        --label "org.opencontainers.image.source=${SOURCE_URL}" \
        -f "$dockerfile" \
        -t "$versioned" \
        -t "$latest" \
        .
    then
        echo
        echo "BUILD FAILED: ${versioned}"
        build_failed+=("$versioned")
        continue
    fi


    # --------------------------------------------------------
    # Validate ps
    # --------------------------------------------------------

    echo
    echo "Testing ps..."

    if ! docker run --rm \
        --entrypoint /bin/sh \
        "$versioned" \
        -c '
            command -v ps >/dev/null 2>&1 || exit 10
            ps -eo pid,ppid,comm >/dev/null 2>&1 || exit 11
        '
    then
        rc=$?

        echo "PS CHECK FAILED: ${versioned} (exit=${rc})"
        echo "IMAGE WILL NOT BE PUSHED."

        ps_failed+=("$versioned")
        continue
    fi

    echo "PASS: ps validation"


    # --------------------------------------------------------
    # Verify GitHub source metadata
    # --------------------------------------------------------

    source_label="$(
        docker image inspect "$versioned" \
            --format '{{ index .Config.Labels "org.opencontainers.image.source" }}'
    )"

    if [[ "$source_label" != "$SOURCE_URL" ]]; then
        echo "ERROR: Incorrect source label:"
        echo "  Found:    ${source_label}"
        echo "  Expected: ${SOURCE_URL}"
        echo "IMAGE WILL NOT BE PUSHED."

        build_failed+=("$versioned [source-label]")
        continue
    fi

    echo "PASS: source label"


    # --------------------------------------------------------
    # Push versioned release
    # --------------------------------------------------------

    echo
    echo "Pushing:"
    echo "  ${versioned}"

    if ! docker push "$versioned"; then
        echo "PUSH FAILED: ${versioned}"
        push_failed+=("$versioned")
        continue
    fi


    # --------------------------------------------------------
    # Push latest
    # --------------------------------------------------------

    echo
    echo "Pushing:"
    echo "  ${latest}"

    if ! docker push "$latest"; then
        echo "PUSH FAILED: ${latest}"
        push_failed+=("$latest")
        continue
    fi


    passed+=("$versioned")

    echo
    echo "SUCCESS:"
    echo "  ${versioned}"
    echo "  ${latest}"

done


# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

echo
echo
echo "============================================================"
echo "SUMMARY"
echo "============================================================"

echo
echo "Built, validated, and pushed: ${#passed[@]}"
for image in "${passed[@]}"; do
    echo "  ${image}"
done

echo
echo "Build failures: ${#build_failed[@]}"
for image in "${build_failed[@]}"; do
    echo "  ${image}"
done

echo
echo "ps failures: ${#ps_failed[@]}"
for image in "${ps_failed[@]}"; do
    echo "  ${image}"
done

echo
echo "Push failures: ${#push_failed[@]}"
for image in "${push_failed[@]}"; do
    echo "  ${image}"
done

echo
echo "Skipped Dockerfiles: ${#skipped[@]}"
for dockerfile in "${skipped[@]}"; do
    echo "  ${dockerfile}"
done

echo


if (( ${#build_failed[@]} > 0 ||
      ${#ps_failed[@]} > 0 ||
      ${#push_failed[@]} > 0 )); then

    echo "One or more images failed."
    exit 1
fi

echo "All images built, validated, and pushed successfully."
