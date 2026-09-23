#!/usr/bin/env bash

VERSION="1.0.0"

successful=()
failed=()
skipped=()

for dockerfile in Dockerfile.*; do

    name="${dockerfile#Dockerfile.}"

    case "$name" in
        base|conda-base|almalinux9.autoactivate)
            continue
            ;;
    esac

    image="ghcr.io/bdrd-genomics/allinone-${name}"

    if docker image inspect "${image}:${VERSION}" >/dev/null 2>&1; then
        echo
        echo "SKIPPING: ${image}:${VERSION} already exists"
        skipped+=("$image")
        continue
    fi

    echo
    echo "=================================================="
    echo "Building: $image:$VERSION"
    echo "Dockerfile: $dockerfile"
    echo "=================================================="

    if docker build \
        --no-cache \
        -f "$dockerfile" \
        -t "$image:$VERSION" \
        -t "$image:latest" \
        .
    then
        echo "SUCCESS: $image"
        successful+=("$image")
    else
        echo "FAILED: $image"
        failed+=("$image")
    fi

done

echo
echo "=================================================="
echo "BUILD SUMMARY"
echo "=================================================="

echo
echo "Skipped (${#skipped[@]}):"
printf '  %s\n' "${skipped[@]}"

echo
echo "Successful (${#successful[@]}):"
printf '  %s\n' "${successful[@]}"

echo
echo "Failed (${#failed[@]}):"
printf '  %s\n' "${failed[@]}"