#!/usr/bin/env bash

set -euo pipefail


###############################################################################
# All-In-One full-test database setup
#
# This script prepares databases used by the full integration test.
#
# Database strategy:
#
#   REAL databases:
#       CheckM
#       CheckM2
#       CheckV
#       BUSCO
#       AMRFinderPlus
#
#   SYNTHETIC / TOY databases:
#       BLAST nucleotide
#       BLAST protein
#       DIAMOND
#       MMseqs2 nucleotide
#       MMseqs2 protein
#       Kraken2
#
# Software strategy:
#
#   1. Prefer Docker if Docker is installed and usable.
#      Tool-specific All-In-One containers are pulled from GHCR.
#
#   2. If Docker is unavailable, use Conda/Mamba.
#
#   3. If neither Docker nor Conda is available, exit with an error.
#
# Usage:
#
#   ./scripts/setup_test_databases.sh
#
# Optional custom database root:
#
#   ./scripts/setup_test_databases.sh /path/to/test_databases
#
###############################################################################


###############################################################################
# Paths
###############################################################################

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

DB_ROOT="${1:-${ROOT_DIR}/tests/databases}"

TOY_DIR="${DB_ROOT}/toy"


###############################################################################
# GHCR images
###############################################################################

GHCR_ROOT="ghcr.io/bdrd-genomics"

BLAST_IMAGE="${AIO_BLAST_IMAGE:-${GHCR_ROOT}/allinone-blast:latest}"
DIAMOND_IMAGE="${AIO_DIAMOND_IMAGE:-${GHCR_ROOT}/allinone-diamond:latest}"
MMSEQS_IMAGE="${AIO_MMSEQS_IMAGE:-${GHCR_ROOT}/allinone-mmseqs2:latest}"
KRAKEN_IMAGE="${AIO_KRAKEN_IMAGE:-${GHCR_ROOT}/allinone-kraken2:latest}"


###############################################################################
# Conda fallback
###############################################################################

TEST_DB_ENV="aio_test_databases"
TEST_DB_ENV_FILE="${ROOT_DIR}/env/test_databases.yml"


###############################################################################
# Initial setup
###############################################################################

mkdir -p "${DB_ROOT}"
mkdir -p "${TOY_DIR}"


echo "============================================================"
echo " All-In-One test database setup"
echo "============================================================"
echo
echo "Repository : ${ROOT_DIR}"
echo "Database   : ${DB_ROOT}"
echo


###############################################################################
# Helper functions
###############################################################################

have_command() {
    command -v "$1" >/dev/null 2>&1
}


section() {
    echo
    echo "============================================================"
    echo " $1"
    echo "============================================================"
}


###############################################################################
# Download helper
###############################################################################

download_file() {

    local url="$1"
    local output="$2"

    if have_command wget; then

        wget \
            -O "${output}" \
            "${url}"

    elif have_command curl; then

        curl \
            -L \
            --fail \
            -o "${output}" \
            "${url}"

    else

        echo
        echo "ERROR: Neither wget nor curl is available."
        echo "One is required for database downloads."
        exit 1

    fi
}


###############################################################################
# Detect Docker / Conda
###############################################################################

BUILD_BACKEND=""

section "Detecting test database build environment"


###############################################################################
# Prefer Docker
###############################################################################

if have_command docker && docker info >/dev/null 2>&1; then

    BUILD_BACKEND="docker"

    echo
    echo "Docker detected and usable."
    echo
    echo "The following All-In-One GHCR images will be used:"
    echo
    echo "  BLAST   : ${BLAST_IMAGE}"
    echo "  DIAMOND : ${DIAMOND_IMAGE}"
    echo "  MMseqs2 : ${MMSEQS_IMAGE}"
    echo "  Kraken2 : ${KRAKEN_IMAGE}"
    echo

    section "Pulling All-In-One database builder images"

    docker pull "${BLAST_IMAGE}"
    docker pull "${DIAMOND_IMAGE}"
    docker pull "${MMSEQS_IMAGE}"
    docker pull "${KRAKEN_IMAGE}"


###############################################################################
# Conda fallback
###############################################################################

elif have_command conda; then

    BUILD_BACKEND="conda"

    echo
    echo "Docker is unavailable."
    echo "Conda detected."
    echo
    echo "Using Conda environment:"
    echo
    echo "  ${TEST_DB_ENV}"
    echo

    if [[ ! -f "${TEST_DB_ENV_FILE}" ]]; then

        echo "ERROR: Conda environment definition does not exist:"
        echo
        echo "  ${TEST_DB_ENV_FILE}"
        echo
        exit 1

    fi


    ###########################################################################
    # Determine whether environment already exists
    ###########################################################################

    if conda env list | awk '{print $1}' | grep -qx "${TEST_DB_ENV}"; then

        echo "Conda environment already exists."
        echo "Skipping environment creation."

    else

        echo "Creating Conda environment..."
        echo

        if have_command mamba; then

            mamba env create \
                -f "${TEST_DB_ENV_FILE}"

        else

            conda env create \
                -f "${TEST_DB_ENV_FILE}"

        fi

    fi


###############################################################################
# Neither available
###############################################################################

else

    echo
    echo "ERROR: No supported test database build backend was found."
    echo
    echo "The All-In-One full test requires either:"
    echo
    echo "  Docker"
    echo
    echo "or:"
    echo
    echo "  Conda/Mamba"
    echo
    echo "Docker is recommended because the test uses the same"
    echo "All-In-One container images used by the pipeline."
    echo
    echo "After installing Docker or Conda, rerun:"
    echo
    echo "  ./scripts/setup_test_databases.sh"
    echo

    exit 1

fi


###############################################################################
# Generic database tool runner
#
# Arguments:
#
#   run_db_tool IMAGE COMMAND [ARGS ...]
#
# Docker:
#   Executes COMMAND inside IMAGE.
#
# Conda:
#   IMAGE is ignored and COMMAND executes through conda run.
###############################################################################

run_db_tool() {

    local image="$1"

    shift

    if [[ "${BUILD_BACKEND}" == "docker" ]]; then

        docker run \
            --rm \
            --user "$(id -u):$(id -g)" \
            -v "${ROOT_DIR}:${ROOT_DIR}" \
            -v "${DB_ROOT}:${DB_ROOT}" \
            -w "${ROOT_DIR}" \
            "${image}" \
            "$@"

    elif [[ "${BUILD_BACKEND}" == "conda" ]]; then

        conda run \
            -n "${TEST_DB_ENV}" \
            "$@"

    else

        echo "ERROR: Database build backend is not configured."
        exit 1

    fi
}


###############################################################################
# Verify tool availability
###############################################################################

section "Checking database build software"

echo "Checking BLAST..."
run_db_tool "${BLAST_IMAGE}" makeblastdb -version >/dev/null

echo "Checking DIAMOND..."
run_db_tool "${DIAMOND_IMAGE}" diamond version >/dev/null

echo "Checking MMseqs2..."
run_db_tool "${MMSEQS_IMAGE}" mmseqs version >/dev/null

echo "Checking Kraken2..."
run_db_tool "${KRAKEN_IMAGE}" kraken2-build --help >/dev/null 2>&1 || true

echo
echo "Database build environment is ready."


###############################################################################
# REAL DATABASES
###############################################################################

section "Installing real test databases"


if [[ ! -x "${ROOT_DIR}/scripts/download_databases.sh" ]]; then

    echo
    echo "ERROR: Database installer was not found or is not executable:"
    echo
    echo "  ${ROOT_DIR}/scripts/download_databases.sh"
    echo

    exit 1

fi


###############################################################################
# CheckM
###############################################################################

echo
echo "[1/5] Installing CheckM database"

"${ROOT_DIR}/scripts/download_databases.sh" \
    --db-root "${DB_ROOT}" \
    --db checkm


###############################################################################
# CheckM2
###############################################################################

echo
echo "[2/5] Installing CheckM2 database"

"${ROOT_DIR}/scripts/download_databases.sh" \
    --db-root "${DB_ROOT}" \
    --db checkm2


###############################################################################
# CheckV
###############################################################################

echo
echo "[3/5] Installing CheckV database"

"${ROOT_DIR}/scripts/download_databases.sh" \
    --db-root "${DB_ROOT}" \
    --db checkv


###############################################################################
# BUSCO
#
# Archive:
#
#   bacteria_odb12.2.2026-05-22.tar.gz
#
# Extracted directory:
#
#   bacteria_odb12.2/
###############################################################################

echo
echo "[4/5] Installing BUSCO database"


BUSCO_ROOT="${DB_ROOT}/busco"

BUSCO_LINEAGE="bacteria_odb12.2"

BUSCO_ARCHIVE="bacteria_odb12.2.2026-05-22.tar.gz"

BUSCO_URL="https://busco-data.ezlab.org/v6/data/lineages/${BUSCO_ARCHIVE}"

BUSCO_DIR="${BUSCO_ROOT}/${BUSCO_LINEAGE}"

BUSCO_ARCHIVE_PATH="${BUSCO_ROOT}/${BUSCO_ARCHIVE}"


mkdir -p "${BUSCO_ROOT}"


if [[ -d "${BUSCO_DIR}" ]]; then

    echo "BUSCO database already installed:"
    echo
    echo "  ${BUSCO_DIR}"
    echo
    echo "Skipping BUSCO download and extraction."

else

    ###########################################################################
    # Download BUSCO
    ###########################################################################

    if [[ -f "${BUSCO_ARCHIVE_PATH}" ]]; then

        echo "BUSCO archive already exists:"
        echo
        echo "  ${BUSCO_ARCHIVE_PATH}"
        echo
        echo "Skipping download."

    else

        echo "Downloading BUSCO:"
        echo
        echo "  ${BUSCO_URL}"
        echo

        download_file \
            "${BUSCO_URL}" \
            "${BUSCO_ARCHIVE_PATH}"

    fi


    ###########################################################################
    # Validate archive
    ###########################################################################

    echo
    echo "Checking BUSCO archive..."

    if ! tar -tzf "${BUSCO_ARCHIVE_PATH}" >/dev/null 2>&1; then

        echo
        echo "ERROR: BUSCO archive is invalid or corrupted:"
        echo
        echo "  ${BUSCO_ARCHIVE_PATH}"
        echo

        exit 1

    fi


    echo "BUSCO archive is valid."


    ###########################################################################
    # Extract archive
    ###########################################################################

    echo
    echo "Extracting BUSCO database..."

    tar \
        -xzf "${BUSCO_ARCHIVE_PATH}" \
        -C "${BUSCO_ROOT}"


    ###########################################################################
    # Verify extraction
    ###########################################################################

    if [[ ! -d "${BUSCO_DIR}" ]]; then

        echo
        echo "ERROR: BUSCO extraction completed but expected directory"
        echo "was not found."
        echo
        echo "Expected:"
        echo
        echo "  ${BUSCO_DIR}"
        echo
        echo "Directory contents:"
        echo

        ls -lah "${BUSCO_ROOT}"

        exit 1

    fi


    echo
    echo "BUSCO database installed successfully:"
    echo
    echo "  ${BUSCO_DIR}"

fi


###############################################################################
# AMRFinderPlus
###############################################################################

echo
echo "[5/5] Installing AMRFinderPlus database"

"${ROOT_DIR}/scripts/download_databases.sh" \
    --db-root "${DB_ROOT}" \
    --db amrfinderplus


###############################################################################
# CREATE TOY NUCLEOTIDE REFERENCE
###############################################################################

section "Creating toy nucleotide reference"


mkdir -p "${TOY_DIR}/fasta"


cat > "${TOY_DIR}/fasta/toy_nt.fna" <<'EOF'
>toy_bacterium_1
ATGCGTACGTTAGCTAGCTACGATCGATCGTAGCTAGCTACGATCGATCGTAGCTAGCATCGATCG
ATCGTAGCTAGCTACGATCGATCGTAGCTAGCTACGATCGATCGTAGCTAGCATCGATCGATCGTA
GCTAGCTACGATCGATCGTAGCTAGCTACGATCGATCGTAGCTAGCATCGATCGATCGTAGCTAGC
>toy_bacterium_2
ATGGCTGCTGCTGATCGATCGTAGCGATCGATGCTAGCTAGCATCGATCGTAGCTAGCATCGATGC
TAGCTAGCATCGATCGTAGCTAGCTAGCATCGATCGATGCTAGCTAGCATCGATCGTAGCTAGCAT
CGATGCTAGCTAGCATCGATCGTAGCTAGCATCGATGCTAGCTAGCATCGATCGTAGCTAGCATCG
>toy_virus_1
ATGAAAGGGTTTCCCAAAGGGTTTCCCAAAGGGTTTCCCAAAGGGTTTCCCAAAGGGTTTCCCAAA
GGGTTTCCCAAAGGGTTTCCCAAAGGGTTTCCCAAAGGGTTTCCCAAAGGGTTTCCCAAAGGGTTT
CCCAAAGGGTTTCCCAAAGGGTTTCCCAAAGGGTTTCCCAAAGGGTTTCCCAAAGGGTTTCCCAAA
EOF


###############################################################################
# CREATE TOY PROTEIN REFERENCE
###############################################################################

cat > "${TOY_DIR}/fasta/toy_nr.faa" <<'EOF'
>toy_protein_1
MKKIGYSAPRQTREAIAQLGADVVVVATGGTDSAALAAAGADVVVVNNAGITPEQARR
>toy_protein_2
MNNIRRVAILGLGLAALATASAAADKPEPTAEEAKKAGVEVVVVDDTPGHADYVAAK
>toy_viral_protein_1
MSTNPKPQRKTKRNTNRRPQDVKFPGGGQIVGGVYLLPRRGPRLGVRAPR
EOF


###############################################################################
# BUILD BLAST DATABASES
###############################################################################

section "Building toy BLAST databases"


mkdir -p "${DB_ROOT}/blast/nt"
mkdir -p "${DB_ROOT}/blast/nr"


###############################################################################
# BLAST nucleotide database
###############################################################################

if [[ -f "${DB_ROOT}/blast/nt/nt.nsq" ]] || \
   [[ -f "${DB_ROOT}/blast/nt/nt.ndb" ]]; then

    echo "Toy BLAST nucleotide database already exists."
    echo "Skipping."

else

    echo "Building toy nucleotide BLAST database using:"
    echo
    echo "  ${BLAST_IMAGE}"
    echo

    run_db_tool \
        "${BLAST_IMAGE}" \
        makeblastdb \
        -in "${TOY_DIR}/fasta/toy_nt.fna" \
        -dbtype nucl \
        -parse_seqids \
        -out "${DB_ROOT}/blast/nt/nt"

fi


###############################################################################
# BLAST protein database
###############################################################################

if [[ -f "${DB_ROOT}/blast/nr/nr.psq" ]] || \
   [[ -f "${DB_ROOT}/blast/nr/nr.pdb" ]]; then

    echo "Toy BLAST protein database already exists."
    echo "Skipping."

else

    echo
    echo "Building toy protein BLAST database using:"
    echo
    echo "  ${BLAST_IMAGE}"
    echo

    run_db_tool \
        "${BLAST_IMAGE}" \
        makeblastdb \
        -in "${TOY_DIR}/fasta/toy_nr.faa" \
        -dbtype prot \
        -parse_seqids \
        -out "${DB_ROOT}/blast/nr/nr"

fi


###############################################################################
# BUILD DIAMOND DATABASE
###############################################################################

section "Building toy DIAMOND database"


mkdir -p "${DB_ROOT}/diamond"


if [[ -f "${DB_ROOT}/diamond/nr.dmnd" ]]; then

    echo "Toy DIAMOND database already exists."
    echo "Skipping."

else

    echo "Building DIAMOND database using:"
    echo
    echo "  ${DIAMOND_IMAGE}"
    echo

    run_db_tool \
        "${DIAMOND_IMAGE}" \
        diamond makedb \
        --in "${TOY_DIR}/fasta/toy_nr.faa" \
        --db "${DB_ROOT}/diamond/nr"

fi


###############################################################################
# BUILD MMSEQS2 DATABASES
###############################################################################

section "Building toy MMseqs2 databases"


mkdir -p "${DB_ROOT}/mmseqs"


###############################################################################
# MMseqs nucleotide database
###############################################################################

if [[ -f "${DB_ROOT}/mmseqs/nt" ]] && \
   [[ -f "${DB_ROOT}/mmseqs/nt.dbtype" ]]; then

    echo "Toy MMseqs nucleotide database already exists."
    echo "Skipping."

else

    rm -f \
        "${DB_ROOT}/mmseqs/nt" \
        "${DB_ROOT}/mmseqs/nt.dbtype" \
        "${DB_ROOT}/mmseqs/nt.index" \
        "${DB_ROOT}/mmseqs/nt.lookup" \
        "${DB_ROOT}/mmseqs/nt.source" \
        2>/dev/null || true


    echo "Building MMseqs nucleotide database using:"
    echo
    echo "  ${MMSEQS_IMAGE}"
    echo


    run_db_tool \
        "${MMSEQS_IMAGE}" \
        mmseqs createdb \
        "${TOY_DIR}/fasta/toy_nt.fna" \
        "${DB_ROOT}/mmseqs/nt"

fi


###############################################################################
# MMseqs protein database
###############################################################################

if [[ -f "${DB_ROOT}/mmseqs/nr" ]] && \
   [[ -f "${DB_ROOT}/mmseqs/nr.dbtype" ]]; then

    echo "Toy MMseqs protein database already exists."
    echo "Skipping."

else

    rm -f \
        "${DB_ROOT}/mmseqs/nr" \
        "${DB_ROOT}/mmseqs/nr.dbtype" \
        "${DB_ROOT}/mmseqs/nr.index" \
        "${DB_ROOT}/mmseqs/nr.lookup" \
        "${DB_ROOT}/mmseqs/nr.source" \
        2>/dev/null || true


    echo
    echo "Building MMseqs protein database using:"
    echo
    echo "  ${MMSEQS_IMAGE}"
    echo


    run_db_tool \
        "${MMSEQS_IMAGE}" \
        mmseqs createdb \
        "${TOY_DIR}/fasta/toy_nr.faa" \
        "${DB_ROOT}/mmseqs/nr"

fi


###############################################################################
# BUILD TINY KRAKEN2 DATABASE
###############################################################################

section "Building tiny Kraken2 database"


KRAKEN_DB="${DB_ROOT}/kraken2/test"

mkdir -p "${KRAKEN_DB}"


if [[ -f "${KRAKEN_DB}/hash.k2d" ]] && \
   [[ -f "${KRAKEN_DB}/opts.k2d" ]] && \
   [[ -f "${KRAKEN_DB}/taxo.k2d" ]]; then

    echo "Toy Kraken2 database already exists."
    echo "Skipping."

else

    ###########################################################################
    # Minimal taxonomy
    ###########################################################################

    mkdir -p "${KRAKEN_DB}/taxonomy"


    cat > "${KRAKEN_DB}/taxonomy/nodes.dmp" <<'EOF'
1	|	1	|	no rank	|		|	0	|	0	|	11	|	1	|	0	|	1	|	0	|		|
2	|	1	|	superkingdom	|		|	0	|	0	|	11	|	1	|	0	|	1	|	0	|		|
562	|	2	|	species	|		|	0	|	0	|	11	|	1	|	0	|	1	|	0	|		|
EOF


    cat > "${KRAKEN_DB}/taxonomy/names.dmp" <<'EOF'
1	|	root	|		|	scientific name	|
2	|	Bacteria	|		|	scientific name	|
562	|	Escherichia coli	|		|	scientific name	|
EOF


    ###########################################################################
    # Small Kraken test sequence
    ###########################################################################

    cat > "${TOY_DIR}/fasta/kraken_test.fna" <<'EOF'
>kraken:taxid|562|toy_ecoli
ATGCGTACGTTAGCTAGCTACGATCGATCGTAGCTAGCTACGATCGATCGTAGCTAGCATCGATCG
ATCGTAGCTAGCTACGATCGATCGTAGCTAGCTACGATCGATCGTAGCTAGCATCGATCGATCGTA
GCTAGCTACGATCGATCGTAGCTAGCTACGATCGATCGTAGCTAGCATCGATCGATCGTAGCTAGC
EOF


    echo "Building Kraken2 database using:"
    echo
    echo "  ${KRAKEN_IMAGE}"
    echo


    ###########################################################################
    # Add sequence
    ###########################################################################

    run_db_tool \
        "${KRAKEN_IMAGE}" \
        kraken2-build \
        --add-to-library "${TOY_DIR}/fasta/kraken_test.fna" \
        --db "${KRAKEN_DB}"


    ###########################################################################
    # Build DB
    ###########################################################################

    run_db_tool \
        "${KRAKEN_IMAGE}" \
        kraken2-build \
        --build \
        --db "${KRAKEN_DB}" \
        --threads 2


    ###########################################################################
    # Clean intermediate files
    ###########################################################################

    run_db_tool \
        "${KRAKEN_IMAGE}" \
        kraken2-build \
        --clean \
        --db "${KRAKEN_DB}"


    ###########################################################################
    # Verify
    ###########################################################################

    if [[ ! -f "${KRAKEN_DB}/hash.k2d" ]]; then

        echo
        echo "ERROR: Kraken2 database build failed."
        echo

        exit 1

    fi

fi


###############################################################################
# SIMPLE TEST TAXONOMY
###############################################################################

section "Creating simple taxonomy resources"


mkdir -p "${DB_ROOT}/taxdb"


cat > "${DB_ROOT}/taxdb/test_accession_taxid.tsv" <<'EOF'
accession	taxid
toy_bacterium_1	562
toy_bacterium_2	562
toy_virus_1	10239
EOF


###############################################################################
# DATABASE INFORMATION MARKER
###############################################################################

cat > "${DB_ROOT}/TEST_DATABASE_INFO.txt" <<EOF
All-In-One lightweight/full functional test database

Generated:
$(date)

Database root:
${DB_ROOT}

Build backend:
${BUILD_BACKEND}

Real databases:
  - CheckM
  - CheckM2
  - CheckV
  - BUSCO
  - AMRFinderPlus

Synthetic databases:
  - BLAST nucleotide
  - BLAST protein
  - DIAMOND protein
  - MMseqs2 nucleotide
  - MMseqs2 protein
  - Kraken2
  - lightweight taxonomy resources

Docker images used when Docker backend is selected:

  BLAST:
    ${BLAST_IMAGE}

  DIAMOND:
    ${DIAMOND_IMAGE}

  MMseqs2:
    ${MMSEQS_IMAGE}

  Kraken2:
    ${KRAKEN_IMAGE}

These databases are intended only for pipeline functional testing.

They MUST NOT be used for biological interpretation.
EOF


###############################################################################
# Finished
###############################################################################

echo
echo "============================================================"
echo " Test database setup complete"
echo "============================================================"
echo
echo "Database root:"
echo
echo "  ${DB_ROOT}"
echo
echo "Build backend:"
echo
echo "  ${BUILD_BACKEND}"
echo
echo "The databases can now be used by the full test."
echo
echo "Example:"
echo
echo "  nextflow run . \\"
echo "      -profile test_full,apptainer \\"
echo "      --database_dir \"${DB_ROOT}\""
echo