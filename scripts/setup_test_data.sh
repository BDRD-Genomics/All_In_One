#!/usr/bin/env bash

set -euo pipefail

###############################################################################
# All-In-One synthetic test data setup
#
# Creates:
#   - Synthetic reference genome containing a known protein-coding region
#   - Synthetic host reference containing only Ns
#   - Protein FASTA matching the embedded coding region
#   - 3,000 paired-end short reads
#   - 500 synthetic long reads
#   - Hybrid dataset containing both short and long reads
#   - Samplesheets for short, long, and hybrid test profiles
#
# Python execution is performed inside the All-In-One python_utils Docker image.
#
# Usage:
#   ./scripts/setup_test_data.sh
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

if ! command -v docker >/dev/null 2>&1; then
    echo "ERROR: Docker is required to generate the synthetic test data."
    exit 1
fi

if ! docker info >/dev/null 2>&1; then
    echo "ERROR: Docker is installed but is not currently usable."
    echo "Check that the Docker daemon is running and your user can access it."
    exit 1
fi

echo "Docker is available."

###############################################################################
# Pull and verify Python utility image
###############################################################################

echo
echo "============================================================"
echo " Pulling Python utility image"
echo "============================================================"
echo

docker pull "${PYTHON_IMAGE}"

docker run --rm "${PYTHON_IMAGE}" \
    python3 -c 'import gzip, random, shutil; from pathlib import Path; print("Python environment OK")'

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
HOST_REFERENCE_LENGTH = 100_000

SHORT_READ_PAIRS = 3_000
SHORT_READ_LENGTH = 150
SHORT_INSERT_SIZE = 350

LONG_READ_COUNT = 500
LONG_READ_MIN = 2_000
LONG_READ_MAX = 4_000

# The sequence is synthetic. The real RefSeq accession is used only so that
# MEGAN's accession-to-taxonomy mapping is exercised during the functional test.
TEST_PROTEIN_ACCESSION = "NP_414543.1"
TEST_PROTEIN_HEADER = f"ref|{TEST_PROTEIN_ACCESSION}|"
TEST_PROTEIN_LENGTH = 300
TEST_CDS_START = 10_000

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

root = Path("${TEST_DIR}")
short_dir = root / "short"
long_dir = root / "long"
hybrid_dir = root / "hybrid"
reference_dir = root / "reference"

for directory in [short_dir, long_dir, hybrid_dir, reference_dir]:
    directory.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

random.seed(SEED)
BASES = "ACGT"
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"

CODONS = {
    "A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT",
    "G": "GGT", "H": "CAT", "I": "ATT", "K": "AAA", "L": "CTG",
    "M": "ATG", "N": "AAT", "P": "CCT", "Q": "CAA", "R": "CGT",
    "S": "TCT", "T": "ACT", "V": "GTT", "W": "TGG", "Y": "TAT",
}


def random_dna(length):
    return "".join(random.choice(BASES) for _ in range(length))


def reverse_complement(sequence):
    table = str.maketrans("ACGT", "TGCA")
    return sequence.translate(table)[::-1]


def circular_slice(sequence, start, length):
    if start + length <= len(sequence):
        return sequence[start:start + length]

    first = sequence[start:]
    remaining = length - len(first)
    return first + sequence[:remaining]


def count_fastq_reads(path):
    lines = 0
    with gzip.open(path, "rt") as handle:
        for _ in handle:
            lines += 1
    return lines // 4


def write_fasta(path, name, sequence):
    with open(path, "w") as handle:
        handle.write(f">{name}\n")
        for i in range(0, len(sequence), 80):
            handle.write(sequence[i:i + 80] + "\n")

# ---------------------------------------------------------------------------
# Create a synthetic protein and embed its coding sequence in the test genome
# ---------------------------------------------------------------------------

protein_rng = random.Random(SEED + 1)

test_protein = "M" + "".join(
    protein_rng.choice(AMINO_ACIDS)
    for _ in range(TEST_PROTEIN_LENGTH - 1)
)

test_cds = "".join(CODONS[aa] for aa in test_protein)

if TEST_CDS_START + len(test_cds) >= GENOME_LENGTH:
    raise RuntimeError("Embedded test CDS does not fit inside synthetic genome")

print()
print(f"Creating synthetic reference genome ({GENOME_LENGTH:,} bp)...")
print(
    f"Embedding {TEST_PROTEIN_LENGTH}-aa test protein at "
    f"genome position {TEST_CDS_START:,}"
)

genome = random_dna(GENOME_LENGTH)
genome = (
    genome[:TEST_CDS_START]
    + test_cds
    + genome[TEST_CDS_START + len(test_cds):]
)

reference_path = reference_dir / "test_reference.fasta"
write_fasta(reference_path, "aio_test_reference", genome)

protein_path = reference_dir / "test_protein.faa"
write_fasta(protein_path, TEST_PROTEIN_HEADER, test_protein)

# ---------------------------------------------------------------------------
# Create host-removal reference containing only Ns
# ---------------------------------------------------------------------------

print(f"Creating synthetic host reference ({HOST_REFERENCE_LENGTH:,} Ns)...")

host_reference_path = reference_dir / "test_host_reference.fasta"
host_reference = "N" * HOST_REFERENCE_LENGTH
write_fasta(host_reference_path, "aio_test_host_reference", host_reference)

# ---------------------------------------------------------------------------
# Create paired-end short reads
# ---------------------------------------------------------------------------

print(
    f"Creating {SHORT_READ_PAIRS:,} paired-end "
    f"{SHORT_READ_LENGTH} bp reads..."
)

r1_path = short_dir / "sample_R1.fastq.gz"
r2_path = short_dir / "sample_R2.fastq.gz"
short_quality = "I" * SHORT_READ_LENGTH

with gzip.open(r1_path, "wt") as r1, gzip.open(r2_path, "wt") as r2:
    for read_number in range(1, SHORT_READ_PAIRS + 1):
        start = random.randrange(GENOME_LENGTH)
        fragment = circular_slice(genome, start, SHORT_INSERT_SIZE)

        read1 = fragment[:SHORT_READ_LENGTH]
        read2 = reverse_complement(fragment[-SHORT_READ_LENGTH:])

        r1.write(
            f"@test_short_{read_number}/1\n"
            f"{read1}\n"
            "+\n"
            f"{short_quality}\n"
        )

        r2.write(
            f"@test_short_{read_number}/2\n"
            f"{read2}\n"
            "+\n"
            f"{short_quality}\n"
        )

# ---------------------------------------------------------------------------
# Create long reads
# ---------------------------------------------------------------------------

print(f"Creating {LONG_READ_COUNT:,} synthetic long reads...")

long_path = long_dir / "sample_long.fastq.gz"

with gzip.open(long_path, "wt") as long_handle:
    for read_number in range(1, LONG_READ_COUNT + 1):
        read_length = random.randint(LONG_READ_MIN, LONG_READ_MAX)
        start = random.randrange(GENOME_LENGTH)
        read = circular_slice(genome, start, read_length)
        quality = "I" * read_length

        long_handle.write(
            f"@test_long_{read_number}\n"
            f"{read}\n"
            "+\n"
            f"{quality}\n"
        )

# ---------------------------------------------------------------------------
# Create hybrid dataset
# ---------------------------------------------------------------------------

print("Creating hybrid test dataset...")

shutil.copy2(r1_path, hybrid_dir / "sample_R1.fastq.gz")
shutil.copy2(r2_path, hybrid_dir / "sample_R2.fastq.gz")
shutil.copy2(long_path, hybrid_dir / "sample_long.fastq.gz")

# ---------------------------------------------------------------------------
# Create samplesheets
# ---------------------------------------------------------------------------

print("Creating samplesheets...")

with open(short_dir / "samplesheet.csv", "w") as handle:
    handle.write(
        "sample_id,fastq_1,fastq_2,long_read\n"
        "test_short,tests/data/short/sample_R1.fastq.gz,"
        "tests/data/short/sample_R2.fastq.gz,No_Read\n"
    )

with open(long_dir / "samplesheet.csv", "w") as handle:
    handle.write(
        "sample_id,fastq_1,fastq_2,long_read\n"
        "test_long,No_Read,No_Read,tests/data/long/sample_long.fastq.gz\n"
    )

with open(hybrid_dir / "samplesheet.csv", "w") as handle:
    handle.write(
        "sample_id,fastq_1,fastq_2,long_read\n"
        "test_hybrid,tests/data/hybrid/sample_R1.fastq.gz,"
        "tests/data/hybrid/sample_R2.fastq.gz,"
        "tests/data/hybrid/sample_long.fastq.gz\n"
    )

# ---------------------------------------------------------------------------
# Validate generated data
# ---------------------------------------------------------------------------

print()
print("Validating generated test data...")

r1_count = count_fastq_reads(r1_path)
r2_count = count_fastq_reads(r2_path)
long_count = count_fastq_reads(long_path)

if r1_count != SHORT_READ_PAIRS:
    raise RuntimeError(
        f"Expected {SHORT_READ_PAIRS} R1 reads, found {r1_count}"
    )

if r2_count != SHORT_READ_PAIRS:
    raise RuntimeError(
        f"Expected {SHORT_READ_PAIRS} R2 reads, found {r2_count}"
    )

if long_count != LONG_READ_COUNT:
    raise RuntimeError(
        f"Expected {LONG_READ_COUNT} long reads, found {long_count}"
    )

for required in [reference_path, host_reference_path, protein_path]:
    if not required.exists() or required.stat().st_size == 0:
        raise RuntimeError(f"Required file was not created: {required}")

print()
print("Synthetic test data created successfully.")
print()
print(f"Reference genome : {GENOME_LENGTH:,} bp")
print(f"Host reference   : {HOST_REFERENCE_LENGTH:,} Ns")
print(f"Test protein     : {TEST_PROTEIN_ACCESSION} ({TEST_PROTEIN_LENGTH} aa)")
print(f"Short read pairs : {SHORT_READ_PAIRS:,}")
print(f"Short read size  : {SHORT_READ_LENGTH} bp")
print(f"Long reads       : {LONG_READ_COUNT:,}")
print(f"Long read range  : {LONG_READ_MIN:,}-{LONG_READ_MAX:,} bp")
print()
print(f"R1 validation    : {r1_count:,} reads")
print(f"R2 validation    : {r2_count:,} reads")
print(f"Long validation  : {long_count:,} reads")
PY

###############################################################################
# Verify expected files exist
###############################################################################

EXPECTED_FILES=(
    "${REFERENCE_DIR}/test_reference.fasta"
    "${REFERENCE_DIR}/test_host_reference.fasta"
    "${REFERENCE_DIR}/test_protein.faa"
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
    if [[ ! -s "${file}" ]]; then
        echo "ERROR: Expected file was not created or is empty:"
        echo "  ${file}"
        exit 1
    fi
done

echo
echo "All expected files were created."
echo
echo "Generated files:"
find "${TEST_DIR}" -maxdepth 2 -type f -print | sort

echo
echo "============================================================"
echo " Test data setup complete"
echo "============================================================"
echo
echo "Synthetic read reference:"
echo "  tests/data/reference/test_reference.fasta"
echo
echo "Host-removal reference:"
echo "  tests/data/reference/test_host_reference.fasta"
echo
echo "Protein reference used by the toy protein databases:"
echo "  tests/data/reference/test_protein.faa"
echo
echo "Short-read test:"
echo "  tests/data/short/samplesheet.csv"
echo
echo "Long-read test:"
echo "  tests/data/long/samplesheet.csv"
echo
echo "Hybrid test:"
echo "  tests/data/hybrid/samplesheet.csv"
echo
