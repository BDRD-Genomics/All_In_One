#!/usr/bin/env bash
set -euo pipefail

DB_ROOT="${AIO_DATABASE_DIR:-}"
PRESET=""
DB_LIST=""
FORCE=0
KEEP_ARCHIVES=0

usage() {
    cat <<'EOF'
Usage:
  download_databases.sh --db-root PATH --preset standard
  download_databases.sh --db-root PATH --db checkm,checkm2,checkv
  download_databases.sh --list

Options:
  --db-root PATH
  --preset NAME        minimal | standard | full
  --db LIST            Comma-separated database names
  --force
  --keep-archives
  --list
  -h, --help

Presets:

  minimal:
    checkm2,checkv,busco

  standard:
    checkm,checkm2,checkv,busco,kraken2,metaphlan4,
    amrfinderplus,mobsuite,plasme

  full:
    standard plus nt,nr

Notes:
  * Downloading does not require the corresponding bioinformatics
    software to be installed.
  * GOTTCHA2 is intentionally excluded.
  * AMRFinderPlus and MOB-suite require a later indexing/setup step.
EOF
}

die() {
    echo "ERROR: $*" >&2
    exit 1
}

log() {
    echo "[$(date '+%F %T')] $*"
}

need() {
    command -v "$1" >/dev/null 2>&1 ||
        die "Required system command '$1' was not found."
}

download() {
    local url="$1"
    local output="$2"

    if [[ -s "$output" && "$FORCE" -eq 0 ]]; then
        log "Already downloaded: $output"
        return
    fi

    wget -c -O "$output" "$url"
}


while [[ $# -gt 0 ]]; do
    case "$1" in
        --db-root)
            DB_ROOT="$2"
            shift 2
            ;;
        --preset)
            PRESET="$2"
            shift 2
            ;;
        --db)
            DB_LIST="$2"
            shift 2
            ;;
        --force)
            FORCE=1
            shift
            ;;
        --keep-archives)
            KEEP_ARCHIVES=1
            shift
            ;;
        --list)
            echo "checkm"
            echo "checkm2"
            echo "checkv"
            echo "busco"
            echo "kraken2"
            echo "metaphlan4"
            echo "amrfinderplus"
            echo "mobsuite"
            echo "plasme"
            echo "nt"
            echo "nr"
            echo "taxonomy"
            exit 0
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            die "Unknown argument: $1"
            ;;
    esac
done


[[ -n "$DB_ROOT" ]] ||
    die "Set --db-root or AIO_DATABASE_DIR"

mkdir -p "$DB_ROOT"

DB_ROOT="$(cd "$DB_ROOT" && pwd)"

export AIO_DATABASE_DIR="$DB_ROOT"

need wget
need tar
need curl

if [[ -n "$DB_LIST" && -n "$PRESET" ]]; then
    die "Use either --db or --preset, not both"
fi


if [[ -z "$DB_LIST" ]]; then

    case "${PRESET:-standard}" in

        minimal)
            DB_LIST="checkm2,checkv,busco"
            ;;

        standard)
            DB_LIST="checkm,checkm2,checkv,busco,kraken2,metaphlan4,amrfinderplus,mobsuite,plasme"
            ;;

        full)
            DB_LIST="checkm,checkm2,checkv,busco,kraken2,metaphlan4,amrfinderplus,mobsuite,plasme,nt,nr"
            ;;

        *)
            die "Unknown preset: $PRESET"
            ;;

    esac

fi


IFS=',' read -r -a DBS <<< "$DB_LIST"


# ============================================================
# CheckM
# ============================================================

download_checkm() {

    local dest="$DB_ROOT/checkm"
    local archive="$dest/checkm_data_2015_01_16.tar.gz"

    local url="https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz"

    mkdir -p "$dest"

    if [[ -f "$dest/hmms/phylo.hmm" && "$FORCE" -eq 0 ]]; then
        log "CheckM already present"
        return
    fi

    log "Downloading CheckM"

    download "$url" "$archive"

    tar -xzf "$archive" -C "$dest"

    [[ -f "$dest/hmms/phylo.hmm" ]] ||
        die "CheckM extraction failed: hmms/phylo.hmm missing"

    [[ "$KEEP_ARCHIVES" -eq 1 ]] ||
        rm -f "$archive"
}


# ============================================================
# CheckM2
# ============================================================

download_checkm2() {

    local dest="$DB_ROOT/checkm2"
    local archive="$dest/checkm2_database.tar.gz"

    local url="https://zenodo.org/records/14897628/files/checkm2_database.tar.gz?download=1"

    mkdir -p "$dest"

    if [[ -f "$dest/CheckM2_database/uniref100.KO.1.dmnd" &&
          "$FORCE" -eq 0 ]]
    then
        log "CheckM2 already present"
        return
    fi

    log "Downloading CheckM2"

    download "$url" "$archive"

    tar -xzf "$archive" -C "$dest"

    [[ -f "$dest/CheckM2_database/uniref100.KO.1.dmnd" ]] ||
        die "CheckM2 extraction did not produce expected database"

    [[ "$KEEP_ARCHIVES" -eq 1 ]] ||
        rm -f "$archive"
}


# ============================================================
# CheckV
# ============================================================

download_checkv() {

    local dest="$DB_ROOT/checkv"
    local archive="$dest/checkv-db-v1.5.tar.gz"

    local url="https://portal.nersc.gov/CheckV/checkv-db-v1.5.tar.gz"

    mkdir -p "$dest"

    if [[ -d "$dest/checkv-db-v1.5" &&
          "$FORCE" -eq 0 ]]
    then
        log "CheckV already present"
        return
    fi

    log "Downloading CheckV"

    download "$url" "$archive"

    tar -xzf "$archive" -C "$dest"

    [[ -d "$dest/checkv-db-v1.5" ]] ||
        die "CheckV extraction failed"

    [[ "$KEEP_ARCHIVES" -eq 1 ]] ||
        rm -f "$archive"
}


# ============================================================
# BUSCO
# ============================================================

download_busco() {

    local dest="$DB_ROOT/busco"

    local dataset="bacteria_odb12.2.2026-05-22"

    local archive="$dest/${dataset}.tar.gz"

    local url="https://busco-data.ezlab.org/v6/data/lineages/${dataset}.tar.gz"

    mkdir -p "$dest"

    if [[ -d "$dest/$dataset" &&
          "$FORCE" -eq 0 ]]
    then
        log "BUSCO dataset already present"
        return
    fi

    log "Downloading BUSCO $dataset"

    download "$url" "$archive"

    tar -xzf "$archive" -C "$dest"

    [[ -d "$dest/$dataset" ]] ||
        die "BUSCO extraction failed"

    [[ "$KEEP_ARCHIVES" -eq 1 ]] ||
        rm -f "$archive"
}


# ============================================================
# Kraken2 / Bracken Standard-8
# ============================================================

download_kraken2() {

    local dest="$DB_ROOT/kraken2/standard_8"

    local archive="$DB_ROOT/kraken2/k2_standard_08_GB_20260626.tar.gz"

    local url="https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08_GB_20260626.tar.gz"

    mkdir -p "$dest"

    if [[ -f "$dest/hash.k2d" &&
          -f "$dest/opts.k2d" &&
          -f "$dest/taxo.k2d" &&
          "$FORCE" -eq 0 ]]
    then
        log "Kraken2 Standard-8 already present"
        return
    fi

    log "Downloading Kraken2 Standard-8"

    download "$url" "$archive"

    tar -xzf "$archive" -C "$dest"

    [[ -f "$dest/hash.k2d" ]] ||
        die "Kraken2 hash.k2d missing"

    [[ -f "$dest/opts.k2d" ]] ||
        die "Kraken2 opts.k2d missing"

    [[ -f "$dest/taxo.k2d" ]] ||
        die "Kraken2 taxo.k2d missing"

    [[ "$KEEP_ARCHIVES" -eq 1 ]] ||
        rm -f "$archive"
}


# ============================================================
# MetaPhlAn
# ============================================================

download_metaphlan4() {

    local dest="$DB_ROOT/metaphlan4"

    local base="https://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases"

    mkdir -p "$dest"

    log "Getting current MetaPhlAn database name"

    wget -q \
        -O "$dest/mpa_latest" \
        "$base/mpa_latest"

    local index

    index="$(
        grep -v '^#' "$dest/mpa_latest" |
        sed '/^[[:space:]]*$/d' |
        head -1
    )"

    [[ -n "$index" ]] ||
        die "Could not determine MetaPhlAn database version"

    local archive="$dest/${index}.tar"

    log "MetaPhlAn database: $index"

    if [[ -f "$dest/${index}.pkl" &&
          "$FORCE" -eq 0 ]]
    then
        log "MetaPhlAn already present"
        return
    fi

    download \
        "$base/${index}.tar" \
        "$archive"

    tar -xf "$archive" -C "$dest"

    # Auxiliary current-database files.
    wget -c \
        -P "$dest" \
        "$base/${index}.nwk" \
        "$base/${index}_marker_info.txt.bz2" \
        "$base/${index}_species.txt.bz2" \
        || true

    [[ -f "$dest/${index}.pkl" ]] ||
        die "MetaPhlAn extraction failed"

    [[ "$KEEP_ARCHIVES" -eq 1 ]] ||
        rm -f "$archive"
}


# ============================================================
# PLASMe
# ============================================================

download_plasme() {

    need unzip

    local root="$DB_ROOT/PLASMe"
    local dest="$root/DB"

    local archive="$root/DB.zip"

    local url="https://zenodo.org/records/8046934/files/DB.zip?download=1"

    mkdir -p "$root"

    if [[ -d "$dest" &&
          -n "$(find "$dest" -mindepth 1 -print -quit 2>/dev/null)" &&
          "$FORCE" -eq 0 ]]
    then
        log "PLASMe already present"
        return
    fi

    log "Downloading PLASMe database"

    download "$url" "$archive"

    (
        cd "$root"
        unzip -o DB.zip
    )

    [[ -d "$dest" ]] ||
        die "PLASMe extraction failed"

    [[ "$KEEP_ARCHIVES" -eq 1 ]] ||
        rm -f "$archive"
}

# ============================================================
# NCBI Taxonomy
# ============================================================

install_ncbi_taxonomy() {
    local taxonomy_dir="${DB_ROOT}/taxonomy"
    local taxdump="${taxonomy_dir}/taxdump.tar.gz"
    local taxa_sqlite="${taxonomy_dir}/taxa.sqlite"
    local traverse="${taxonomy_dir}/taxa.sqlite.traverse.pkl"

    echo "============================================================"
    echo "Installing NCBI taxonomy / ETE3 taxa.sqlite"
    echo "============================================================"

    mkdir -p "${taxonomy_dir}"

    # Download the official NCBI taxonomy archive
    if [[ ! -s "${taxdump}" ]]; then
        curl -L \
            --retry 3 \
            -o "${taxdump}" \
            "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
    else
        echo "NCBI taxdump already exists; skipping download."
    fi

    # Build the ETE3 SQLite database
    if [[ ! -s "${taxa_sqlite}" ]]; then

        # Remove any partial files left from a failed build
        rm -f "${taxa_sqlite}" "${traverse}"

        docker run --rm -i \
            --user "$(id -u):$(id -g)" \
            -v "${DB_ROOT}:${DB_ROOT}" \
            -w /tmp \
            ghcr.io/bdrd-genomics/allinone-python_utils:latest \
            python - "${taxdump}" "${taxa_sqlite}" <<'PY'
import sys

from ete3 import NCBITaxa

taxdump = sys.argv[1]
dbfile = sys.argv[2]

print("Building ETE3 taxonomy database:")
print(f"  taxdump: {taxdump}")
print(f"  output : {dbfile}")

try:
    ncbi = NCBITaxa(
        dbfile=dbfile,
        taxdump_file=taxdump
    )
except TypeError:
    # Compatibility with older ETE3 releases
    from ete3.ncbi_taxonomy.ncbiquery import update_db

    update_db(dbfile, taxdump)
    ncbi = NCBITaxa(dbfile=dbfile)

# Basic validation
names = ncbi.get_taxid_translator([9606])

if 9606 not in names:
    raise RuntimeError(
        "ETE3 taxonomy validation failed: taxid 9606 could not be resolved"
    )

print(f"Validation successful: 9606 -> {names[9606]}")
print(f"Created: {dbfile}")
PY

    else
        echo "ETE3 taxa.sqlite already exists; skipping build."
    fi

    if [[ ! -s "${taxa_sqlite}" ]]; then
        echo "ERROR: Failed to create ${taxa_sqlite}" >&2
        exit 1
    fi

    echo "NCBI taxonomy database ready:"
    echo "  ${taxa_sqlite}"
}

# ============================================================
# MOB-suite raw database
# ============================================================

download_mobsuite() {

    local dest="$DB_ROOT/mob_suite"
    local archive="$dest/data.tar.gz"

    local url="https://zenodo.org/records/10304948/files/data.tar.gz?download=1"

    mkdir -p "$dest"

    log "Downloading MOB-suite database payload"

    download "$url" "$archive"

    tar -xzf "$archive" -C "$dest"

    [[ "$KEEP_ARCHIVES" -eq 1 ]] ||
        rm -f "$archive"

    log "NOTE: MOB-suite payload downloaded."
    log "MOB-suite still requires database indexing before use."
}


# ============================================================
# AMRFinderPlus
# ============================================================

download_amrfinderplus() {

    local dest="$DB_ROOT/amrfinderplus"

    mkdir -p "$dest"

    log "AMRFinderPlus cannot be made completely ready with wget alone."
    log "NCBI distributes the source database files and amrfinder_index"
    log "builds the final searchable indices."

    touch "$dest/REQUIRES_INDEXING"

    log "Skipping AMRFinderPlus initialization for download-only mode."
}


# ============================================================
# BLAST nt / nr
# ============================================================

download_blast_db() {

    local name="$1"

    local dest="$DB_ROOT/blastdb/$name"

    mkdir -p "$dest"

    log "$name is intentionally not downloaded by this wget-only"
    log "test installer."

    log "Use the production database installation workflow for $name."
}


install_one() {

    case "$1" in

        checkm)
            download_checkm
            ;;

        checkm2)
            download_checkm2
            ;;

        checkv)
            download_checkv
            ;;

        busco)
            download_busco
            ;;

        kraken2)
            download_kraken2
            ;;

        metaphlan4)
            download_metaphlan4
            ;;

        plasme)
            download_plasme
            ;;

        mobsuite)
            download_mobsuite
            ;;

        amrfinderplus)
            download_amrfinderplus
            ;;

        nt|nr)
            download_blast_db "$1"
            ;;
        taxonomy)
            install_ncbi_taxonomy
            ;;
        *)
            die "Unknown database: $1"
            ;;

    esac
}


log "Database root: $DB_ROOT"
log "Requested: ${DBS[*]}"

for db in "${DBS[@]}"; do

    db="${db//[[:space:]]/}"

    [[ -n "$db" ]] || continue

    log "============================================================"
    log "Installing: $db"
    log "============================================================"

    install_one "$db"

done

log "============================================================"
log "Download complete"
log "============================================================"

log "AIO_DATABASE_DIR=$DB_ROOT"