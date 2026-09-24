#!/usr/bin/env bash
set -euo pipefail

DB_ROOT="${AIO_DATABASE_DIR:-}"

usage() {
    echo "Usage: verify_databases.sh --db-root PATH"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --db-root) DB_ROOT="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

[[ -n "$DB_ROOT" ]] || { echo "Set --db-root or AIO_DATABASE_DIR." >&2; exit 1; }

ok=0
missing=0

check_path() {
    local name="$1"
    local path="$2"
    if compgen -G "$path" >/dev/null 2>&1; then
        printf "OK      %-16s %s\n" "$name" "$path"
        ok=$((ok+1))
    else
        printf "MISSING %-16s %s\n" "$name" "$path"
        missing=$((missing+1))
    fi
}

check_path taxdb          "$DB_ROOT/blast/taxdb/taxdb.*"
check_path nt             "$DB_ROOT/blast/nt/nt.*"
check_path nr             "$DB_ROOT/blast/nr/nr.*"
check_path checkm         "$DB_ROOT/checkm/marker_sets/*"
if find "$DB_ROOT/checkm2" -type f -name 'uniref100.KO.1.dmnd' -print -quit 2>/dev/null | grep -q .; then
    printf "OK      %-16s %s\n" "checkm2" "$DB_ROOT/checkm2/.../uniref100.KO.1.dmnd"
    ok=$((ok+1))
else
    printf "MISSING %-16s %s\n" "checkm2" "$DB_ROOT/checkm2/.../uniref100.KO.1.dmnd"
    missing=$((missing+1))
fi
check_path checkv         "$DB_ROOT/checkv/*"
check_path busco          "$DB_ROOT/busco/*"
check_path kraken2        "$DB_ROOT/kraken2/standard/hash.k2d"
check_path metaphlan4     "$DB_ROOT/metaphlan4/*"
check_path amrfinderplus  "$DB_ROOT/amrfinderplus/latest/*"

echo
echo "Present checks: $ok"
echo "Missing checks: $missing"
echo
echo "Note: missing nt/nr is normal unless you intentionally installed the full preset."

[[ "$missing" -eq 0 ]] || exit 2
