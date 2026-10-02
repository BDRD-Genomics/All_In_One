#!/usr/bin/env bash

set -euo pipefail


###############################################################################
# All-In-One synthetic test data setup
#
# Creates:
#
#   - Synthetic reference genome
#   - 3,000 paired-end short reads
#   - 500 synthetic long reads
#   - Hybrid dataset containing both short and long reads
#   - Samplesheets for short, long, and hybrid test profiles
#
# Python execution is performed inside the All-In-One python_utils Docker image.
#
# Usage:
#
#   ./scripts/setup_test_data.sh
#
###############################################################################


###############################################################################
# Configuration
###############################################################################

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

TEST_DIR="${ROOT_DIR}/tests/data"

SHORT_DIR="${TEST_DIR}/short"
LONG_DIR="${TEST_DIR}/long"
HYBRID_DIR="${TEST_DIR}/hybrid"
REFERENCE_DIR="${TEST_DIR}/reference"

PYTHON_IMAGE="${AIO_PYTHON_IMAGE:-ghcr.io/bdrd-genomics/allinone-python_utils:latest}"


###############################################################################
# Header
###############################################################################

echo "============================================================"
echo " All-In-One synthetic test data setup"
echo "============================================================"
echo
echo "Repository : ${ROOT_DIR}"
echo "Test data  : ${TEST_DIR}"
echo "Image      : ${PYTHON_IMAGE}"
echo


###############################################################################
# Check Docker
###############################################################################

echo "============================================================"
echo " Checking Docker"
echo "============================================================"
echo

if ! command -v docker >/dev/null 2>&1; then

    echo "ERROR: Docker is required to generate the synthetic test data."
    echo
    echo "Install Docker and rerun:"
    echo
    echo "  ./scripts/setup_test_data.sh"
    echo

    exit 1

fi


if ! docker info >/dev/null 2>&1; then

    echo "ERROR: Docker is installed but is not currently usable."
    echo
    echo "Check that:"
    echo
    echo "  - the Docker daemon is running"
    echo "  - your user has permission to use Docker"
    echo
    echo "Then rerun:"
    echo
    echo "  ./scripts/setup_test_data.sh"
    echo

    exit 1

fi


echo "Docker is available."


###############################################################################
# Pull Python utility image
###############################################################################

echo
echo "============================================================"
echo " Pulling Python utility image"
echo "============================================================"
echo
echo "${PYTHON_IMAGE}"
echo


if ! docker pull "${PYTHON_IMAGE}"; then

    echo
    echo "ERROR: Unable to pull:"
    echo
    echo "  ${PYTHON_IMAGE}"
    echo
    echo "If the GHCR package is private, authenticate first:"
    echo
    echo "  docker login ghcr.io"
    echo

    exit 1

fi


###############################################################################
# Verify Python container
###############################################################################

echo
echo "============================================================"
echo " Checking Python environment"
echo "============================================================"
echo


if ! docker run \
    --rm \
    "${PYTHON_IMAGE}" \
    python3 -c 'import gzip, random, shutil; from pathlib import Path; print("Python environment OK")'
then

    echo
    echo "ERROR: The Python utility image does not provide the required"
    echo "Python environment."
    echo

    exit 1

fi


###############################################################################
# Create output directories
###############################################################################

mkdir -p \
    "${SHORT_DIR}" \
    "${LONG_DIR}" \
    "${HYBRID_DIR}" \
    "${REFERENCE_DIR}"


###############################################################################
# Generate synthetic test data
###############################################################################

echo
echo "============================================================"
echo " Creating All-In-One synthetic test data"
echo "============================================================"
echo


docker run \
    --rm \
    -i \
    --user "$(id -u):$(id -g)" \
    -v "${ROOT_DIR}:${ROOT_DIR}" \
    -w "${ROOT_DIR}" \
    "${PYTHON_IMAGE}" \
    python3 - <<PY

import gzip
import random
import shutil
from pathlib import Path


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

SEED = 42

GENOME_LENGTH = 100_000

SHORT_READ_PAIRS = 3_000
SHORT_READ_LENGTH = 150
SHORT_INSERT_SIZE = 350

LONG_READ_COUNT = 500
LONG_READ_MIN = 2_000
LONG_READ_MAX = 4_000


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

root = Path("${TEST_DIR}")

short_dir = root / "short"
long_dir = root / "long"
hybrid_dir = root / "hybrid"
reference_dir = root / "reference"


for directory in [
    short_dir,
    long_dir,
    hybrid_dir,
    reference_dir,
]:
    directory.mkdir(
        parents=True,
        exist_ok=True,
    )


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

random.seed(SEED)

BASES = "ACGT"


def random_dna(length):
    return "".join(
        random.choice(BASES)
        for _ in range(length)
    )


def reverse_complement(sequence):
    table = str.maketrans(
        "ACGT",
        "TGCA",
    )

    return sequence.translate(table)[::-1]


def circular_slice(sequence, start, length):
    """
    Extract sequence from a circular synthetic genome.
    """

    if start + length <= len(sequence):
        return sequence[start:start + length]

    first = sequence[start:]

    remaining = length - len(first)

    return first + sequence[:remaining]


def count_fastq_reads(path):
    """
    Count reads in a gzipped FASTQ.
    """

    lines = 0

    with gzip.open(path, "rt") as handle:
        for _ in handle:
            lines += 1

    return lines // 4


# ---------------------------------------------------------------------------
# Create synthetic reference genome
# ---------------------------------------------------------------------------

print()
print(
    f"Creating synthetic reference genome "
    f"({GENOME_LENGTH:,} bp)..."
)

genome = random_dna(
    GENOME_LENGTH
)

reference_path = (
    reference_dir /
    "test_reference.fasta"
)


with open(reference_path, "w") as handle:

    handle.write(
        ">aio_test_reference\\n"
    )

    for i in range(
        0,
        len(genome),
        80,
    ):
        handle.write(
            genome[i:i + 80] + "\\n"
        )


# ---------------------------------------------------------------------------
# Create paired-end short reads
# ---------------------------------------------------------------------------

print(
    f"Creating {SHORT_READ_PAIRS:,} paired-end "
    f"{SHORT_READ_LENGTH} bp reads..."
)


r1_path = (
    short_dir /
    "sample_R1.fastq.gz"
)

r2_path = (
    short_dir /
    "sample_R2.fastq.gz"
)


short_quality = (
    "I" * SHORT_READ_LENGTH
)


with gzip.open(r1_path, "wt") as r1, \
     gzip.open(r2_path, "wt") as r2:

    for read_number in range(
        1,
        SHORT_READ_PAIRS + 1,
    ):

        start = random.randrange(
            GENOME_LENGTH
        )

        fragment = circular_slice(
            genome,
            start,
            SHORT_INSERT_SIZE,
        )

        read1 = fragment[
            :SHORT_READ_LENGTH
        ]

        read2 = reverse_complement(
            fragment[
                -SHORT_READ_LENGTH:
            ]
        )

        r1.write(
            f"@test_short_{read_number}/1\\n"
            f"{read1}\\n"
            "+\\n"
            f"{short_quality}\\n"
        )

        r2.write(
            f"@test_short_{read_number}/2\\n"
            f"{read2}\\n"
            "+\\n"
            f"{short_quality}\\n"
        )


# ---------------------------------------------------------------------------
# Create long reads
# ---------------------------------------------------------------------------

print(
    f"Creating {LONG_READ_COUNT:,} "
    f"synthetic long reads..."
)


long_path = (
    long_dir /
    "sample_long.fastq.gz"
)


with gzip.open(
    long_path,
    "wt",
) as long_handle:

    for read_number in range(
        1,
        LONG_READ_COUNT + 1,
    ):

        read_length = random.randint(
            LONG_READ_MIN,
            LONG_READ_MAX,
        )

        start = random.randrange(
            GENOME_LENGTH
        )

        read = circular_slice(
            genome,
            start,
            read_length,
        )

        quality = (
            "I" * read_length
        )

        long_handle.write(
            f"@test_long_{read_number}\\n"
            f"{read}\\n"
            "+\\n"
            f"{quality}\\n"
        )


# ---------------------------------------------------------------------------
# Create hybrid dataset
# ---------------------------------------------------------------------------

print(
    "Creating hybrid test dataset..."
)


shutil.copy2(
    r1_path,
    hybrid_dir /
    "sample_R1.fastq.gz",
)

shutil.copy2(
    r2_path,
    hybrid_dir /
    "sample_R2.fastq.gz",
)

shutil.copy2(
    long_path,
    hybrid_dir /
    "sample_long.fastq.gz",
)


# ---------------------------------------------------------------------------
# Create samplesheets
# ---------------------------------------------------------------------------

print(
    "Creating samplesheets..."
)


# Short-read samplesheet

with open(
    short_dir / "samplesheet.csv",
    "w",
) as handle:

    handle.write(
        "sample_id,fastq_1,fastq_2,long_read\\n"
        "test_short,"
        "tests/data/short/sample_R1.fastq.gz,"
        "tests/data/short/sample_R2.fastq.gz,"
        "No_Read\\n"
    )


# Long-read samplesheet

with open(
    long_dir / "samplesheet.csv",
    "w",
) as handle:

    handle.write(
        "sample_id,fastq_1,fastq_2,long_read\\n"
        "test_long,"
        "No_Read,"
        "No_Read,"
        "tests/data/long/sample_long.fastq.gz\\n"
    )


# Hybrid samplesheet

with open(
    hybrid_dir / "samplesheet.csv",
    "w",
) as handle:

    handle.write(
        "sample_id,fastq_1,fastq_2,long_read\\n"
        "test_hybrid,"
        "tests/data/hybrid/sample_R1.fastq.gz,"
        "tests/data/hybrid/sample_R2.fastq.gz,"
        "tests/data/hybrid/sample_long.fastq.gz\\n"
    )


# ---------------------------------------------------------------------------
# Validate generated data
# ---------------------------------------------------------------------------

print()
print("Validating generated test data...")


r1_count = count_fastq_reads(
    r1_path
)

r2_count = count_fastq_reads(
    r2_path
)

long_count = count_fastq_reads(
    long_path
)


if r1_count != SHORT_READ_PAIRS:
    raise RuntimeError(
        f"Expected {SHORT_READ_PAIRS} R1 reads, "
        f"found {r1_count}"
    )


if r2_count != SHORT_READ_PAIRS:
    raise RuntimeError(
        f"Expected {SHORT_READ_PAIRS} R2 reads, "
        f"found {r2_count}"
    )


if long_count != LONG_READ_COUNT:
    raise RuntimeError(
        f"Expected {LONG_READ_COUNT} long reads, "
        f"found {long_count}"
    )


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

print()
print(
    "Synthetic test data created successfully."
)

print()
print(
    f"Reference genome : "
    f"{GENOME_LENGTH:,} bp"
)

print(
    f"Short read pairs : "
    f"{SHORT_READ_PAIRS:,}"
)

print(
    f"Short read size  : "
    f"{SHORT_READ_LENGTH} bp"
)

print(
    f"Long reads       : "
    f"{LONG_READ_COUNT:,}"
)

print(
    f"Long read range  : "
    f"{LONG_READ_MIN:,}-"
    f"{LONG_READ_MAX:,} bp"
)

print()
print(
    f"R1 validation   : "
    f"{r1_count:,} reads"
)

print(
    f"R2 validation   : "
    f"{r2_count:,} reads"
)

print(
    f"Long validation : "
    f"{long_count:,} reads"
)

PY


###############################################################################
# Verify expected files exist
###############################################################################

echo
echo "============================================================"
echo " Verifying generated files"
echo "============================================================"
echo


EXPECTED_FILES=(
    "${REFERENCE_DIR}/test_reference.fasta"

    "${SHORT_DIR}/sample_R1.fastq.gz"
    "${SHORT_DIR}/sample_R2.fastq.gz"
    "${SHORT_DIR}/samplesheet.csv"

    "${LONG_DIR}/sample_long.fastq.gz"
    "${LONG_DIR}/samplesheet.csv"

    "${HYBRID_DIR}/sample_R1.fastq.gz"
    "${HYBRID_DIR}/sample_R2.fastq.gz"
    "${HYBRID_DIR}/sample_long.fastq.gz"
    "${HYBRID_DIR}/samplesheet.csv"
)


for file in "${EXPECTED_FILES[@]}"; do

    if [[ ! -f "${file}" ]]; then

        echo "ERROR: Expected file was not created:"
        echo
        echo "  ${file}"
        echo

        exit 1

    fi

done


echo "All expected files were created."


###############################################################################
# Display generated files
###############################################################################

echo
echo "Generated files:"
echo

find "${TEST_DIR}" \
    -maxdepth 2 \
    -type f \
    -print \
    | sort


###############################################################################
# Finished
###############################################################################

echo
echo "============================================================"
echo " Test data setup complete"
echo "============================================================"
echo

echo "Reference genome:"
echo
echo "  tests/data/reference/test_reference.fasta"
echo

echo "Short-read test:"
echo
echo "  tests/data/short/samplesheet.csv"
echo

echo "Long-read test:"
echo
echo "  tests/data/long/samplesheet.csv"
echo

echo "Hybrid test:"
echo
echo "  tests/data/hybrid/samplesheet.csv"
echo

echo "The synthetic short and long reads originate from the same"
echo "reference genome, allowing the data to be used for QC,"
echo "mapping, assembly, and hybrid integration tests."
echo