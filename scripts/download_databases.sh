#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
MANIFEST="${REPO_ROOT}/databases/manifest.tsv"

DB_ROOT="${AIO_DATABASE_DIR:-}"
PRESET=""
DB_LIST=""
THREADS="${AIO_DB_THREADS:-8}"
FORCE=0
KEEP_ARCHIVES=0

usage() {
    cat <<'EOF'
Usage:
  download_databases.sh --db-root PATH --preset standard
  download_databases.sh --db-root PATH --preset full
  download_databases.sh --db-root PATH --db checkm2,checkv,taxdb
  download_databases.sh --list

Options:
  --db-root PATH       Database root. Can also use AIO_DATABASE_DIR.
  --preset NAME        minimal | standard | full
  --db LIST            Comma-separated database names.
  --threads N          Build/download threads where supported. Default: 8
  --force              Re-download/rebuild where supported.
  --keep-archives      Keep downloaded archives.
  --list               Show available database names.
  -h, --help           Show this help.

Presets:
  minimal   taxdb,checkm2,checkv
  standard  minimal + checkm,busco,kraken2,metaphlan4,amrfinderplus
  full      standard + nt,nr

Notes:
  * nt and nr are never downloaded by default.
  * Database downloads are intentionally separate from Nextflow execution.
  * Upstream tool commands must be available in PATH for tool-managed databases.
EOF
}

die() { echo "ERROR: $*" >&2; exit 1; }
log() { echo "[$(date '+%F %T')] $*"; }
have() { command -v "$1" >/dev/null 2>&1; }

list_dbs() {
    awk -F'\t' 'NR>1 {printf "%-16s %-14s %-10s %s\n",$1,$2,$3,$8}' "$MANIFEST"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --db-root) DB_ROOT="$2"; shift 2 ;;
        --preset) PRESET="$2"; shift 2 ;;
        --db) DB_LIST="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        --force) FORCE=1; shift ;;
        --keep-archives) KEEP_ARCHIVES=1; shift ;;
        --list) list_dbs; exit 0 ;;
        -h|--help) usage; exit 0 ;;
        *) die "Unknown argument: $1" ;;
    esac
done

[[ -n "$DB_ROOT" ]] || die "Set --db-root or AIO_DATABASE_DIR."
mkdir -p "$DB_ROOT"
DB_ROOT="$(cd "$DB_ROOT" && pwd)"
export AIO_DATABASE_DIR="$DB_ROOT"

if [[ -n "$DB_LIST" && -n "$PRESET" ]]; then
    die "Use either --db or --preset, not both."
fi

if [[ -z "$DB_LIST" ]]; then
    case "${PRESET:-standard}" in
        minimal)
            DB_LIST="taxdb,checkm2,checkv"
            ;;
        standard)
            DB_LIST="taxdb,checkm,checkm2,checkv,busco,kraken2,metaphlan4,amrfinderplus"
            ;;
        full)
            DB_LIST="taxdb,checkm,checkm2,checkv,busco,kraken2,metaphlan4,amrfinderplus,nt,nr"
            ;;
        *)
            die "Unknown preset: ${PRESET}"
            ;;
    esac
fi

IFS=',' read -r -a DBS <<< "$DB_LIST"

require_cmd() {
    have "$1" || die "Required command '$1' was not found in PATH."
}

download_blast_db() {
    local name="$1"
    require_cmd update_blastdb.pl
    local dest="$DB_ROOT/blast/$name"
    mkdir -p "$dest"
    log "Installing NCBI BLAST database: $name -> $dest"
    (
        cd "$dest"
        update_blastdb.pl --decompress "$name"
    )
}

download_checkm() {
    require_cmd wget
    require_cmd tar
    local dest="$DB_ROOT/checkm"
    local url="https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz"
    local archive="$dest/checkm_data_2015_01_16.tar.gz"
    mkdir -p "$dest"

    if [[ -f "$dest/marker_sets/taxon_marker_sets.tsv" && "$FORCE" -eq 0 ]]; then
        log "CheckM data appears present; skipping. Use --force to reinstall."
        return
    fi

    log "Downloading CheckM v1 database"
    wget -c -O "$archive" "$url"
    tar -xzf "$archive" -C "$dest"
    [[ "$KEEP_ARCHIVES" -eq 1 ]] || rm -f "$archive"

    if have checkm; then
        log "Registering CheckM data root"
        checkm data setRoot "$dest" || true
    else
        log "checkm executable not found; data downloaded but setRoot was not run."
    fi
}

download_checkm2() {
    require_cmd checkm2
    local dest="$DB_ROOT/checkm2"
    mkdir -p "$dest"
    log "Downloading CheckM2 database -> $dest"
    checkm2 database --download --path "$dest"
}

download_checkv() {
    require_cmd checkv
    local dest="$DB_ROOT/checkv"
    mkdir -p "$dest"
    log "Downloading CheckV database -> $dest"
    checkv download_database "$dest"
}

download_busco() {
    require_cmd busco
    local dest="$DB_ROOT/busco"
    mkdir -p "$dest"
    log "Downloading BUSCO prokaryota datasets -> $dest"
    busco --download prokaryota --download_path "$dest"
}

download_kraken2_core_nt() {
    local dest="$DB_ROOT/kraken2/core_nt"
    mkdir -p "$dest"

    local url="CURRENT_CORE_NT_TAR_GZ_URL"
    local archive="$dest/core_nt.tar.gz"

    wget -c -O "$archive" "$url"
    tar -xzf "$archive" -C "$dest"

    [[ "$KEEP_ARCHIVES" -eq 1 ]] || rm -f "$archive"
}

download_metaphlan4() {
    require_cmd metaphlan
    local dest="$DB_ROOT/metaphlan4"
    mkdir -p "$dest"
    log "Installing MetaPhlAn database -> $dest"

    # MetaPhlAn 4.1+ uses --db_dir. Older releases used --bowtie2db.
    if metaphlan --help 2>&1 | grep -q -- '--db_dir'; then
        metaphlan --install --db_dir "$dest"
    else
        metaphlan --install --bowtie2db "$dest"
    fi
}

download_amrfinderplus() {
    require_cmd amrfinder
    local dest="$DB_ROOT/amrfinderplus"
    mkdir -p "$dest"

    log "Updating AMRFinderPlus database using upstream updater"
    if [[ "$FORCE" -eq 1 ]]; then
        amrfinder -U
    else
        amrfinder -u
    fi

    # Ask AMRFinder where its compatible DB was installed.
    local dbdir
    dbdir="$(amrfinder --database_version 2>&1 | sed -n 's/^Database directory:[[:space:]]*//p' | head -1 || true)"

    if [[ -z "$dbdir" || ! -d "$dbdir" ]]; then
        log "AMRFinder updated successfully, but its database location could not be auto-detected."
        log "Run: amrfinder --database_version"
        log "Then copy/symlink that database into: $dest"
        return
    fi

    log "Copying AMRFinderPlus database from $dbdir -> $dest"
    rm -rf "$dest"/latest "$dest"/db
    mkdir -p "$dest"
    cp -a "$dbdir" "$dest/db"
    ln -sfn "$dest/db" "$dest/latest"
}

install_one() {
    case "$1" in
        taxdb|nt|nr) download_blast_db "$1" ;;
        checkm) download_checkm ;;
        checkm2) download_checkm2 ;;
        checkv) download_checkv ;;
        busco) download_busco ;;
        kraken2) download_kraken2 ;;
        metaphlan4) download_metaphlan4 ;;
        amrfinderplus) download_amrfinderplus ;;
        *) die "Unknown database '$1'. Use --list." ;;
    esac
}

log "Database root: $DB_ROOT"
log "Requested: ${DBS[*]}"
log "Threads: $THREADS"

for db in "${DBS[@]}"; do
    db="${db//[[:space:]]/}"
    [[ -n "$db" ]] || continue
    install_one "$db"
done

log "Database installation complete."
log "Run scripts/verify_databases.sh --db-root '$DB_ROOT' to verify."
