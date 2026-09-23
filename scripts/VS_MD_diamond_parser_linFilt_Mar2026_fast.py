#!/usr/bin/env python3
"""
Optimized VS/MD DIAMOND/MMseqs parser

Changes vs original:
- Avoids per-row SQLite lookups by batching unique accession queries
"""

from __future__ import annotations

import argparse
import os
import re
import sqlite3
from functools import lru_cache
from typing import Dict, Iterable, List, Optional, Tuple

import pandas as pd
from ete3 import NCBITaxa


CUSTOM_FIELDS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "sallgi", "qlen", "slen"
]

NUMERIC_COLS = [
    "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send",
    "evalue", "bitscore", "qlen", "slen"
]

NCBI: Optional[NCBITaxa] = None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="VS_MD_parser",
        description="Parse blast table",
        usage="python3 VS_MD_diamond_parser_linFilt_Mar2026_fast.py -i <blast_input> -n <ncbi_database> -v <vhunter_database>",
        epilog="Parse the blast table from sequence comparative algorithm (MMSeqs, blastn, blastx).",
    )
    parser.add_argument("-i", "--input", required=False, help="Input BLAST/MMseqs tabular file")
    parser.add_argument("-s", "--sample_id", required=False, help="Optional sample ID override")
    parser.add_argument("-t", "--type", required=False, choices=["mmseqs", "blastx"])
    parser.add_argument("-r", "--rank", required=False, help="Kept for CLI compatibility")
    parser.add_argument("-o", "--outdir", required=False, default=".")
    parser.add_argument("-n", "--ncbidb", required=False, help="NCBITaxa sqlite dbfile path")
    parser.add_argument("-v", "--vhunter", required=False, help="VHunter sqlite mapping DB path")
    parser.add_argument("--no-sort", action="store_true", help="Skip full DataFrame sort for faster runtime")
    parser.add_argument("--batch-size", type=int, default=50000, help="Batch size for SQLite accession lookups")
    return parser.parse_args()


def chunked(seq: List[str], size: int) -> Iterable[List[str]]:
    for i in range(0, len(seq), size):
        yield seq[i:i + size]


def detect_taxid_column(cursor: sqlite3.Cursor, table: str) -> Optional[str]:
    """Try to find an explicit taxid column. If not found, caller can fall back to row[2]."""
    try:
        rows = cursor.execute(f"PRAGMA table_info({table})").fetchall()
    except Exception:
        return None

    cols = [r[1] for r in rows]
    for candidate in ("taxid", "TaxID", "tax_id"):
        if candidate in cols:
            return candidate
    return None


def fetch_accession_taxids(
    cursor: sqlite3.Cursor,
    accessions: List[str],
    table: str,
    batch_size: int
) -> Dict[str, int]:
    """
    Batch accession -> taxid lookup.
    If no explicit taxid column is present, fall back to original-script behavior:
    SELECT * ... and use row[2] as taxid.
    """
    mapping: Dict[str, int] = {}

    if not accessions:
        return mapping

    taxid_col = detect_taxid_column(cursor, table)

    for chunk in chunked(accessions, batch_size):
        placeholders = ",".join(["?"] * len(chunk))

        if taxid_col is not None:
            query = (
                f"SELECT ACCESSION, {taxid_col} "
                f"FROM {table} WHERE Accession IN ({placeholders})"
            )
            rows = cursor.execute(query, chunk)
            for Accession, taxid in rows:
                try:
                    if Accession is not None and taxid is not None:
                        mapping[str(Accession)] = int(taxid)
                except Exception:
                    continue
        else:
            query = f"SELECT * FROM {table} WHERE Accession IN ({placeholders})"
            rows = cursor.execute(query, chunk)
            for row in rows:
                try:
                    Accession = row[0]
                    taxid = row[2]
                    if Accession is not None and taxid is not None:
                        mapping[str(Accession)] = int(taxid)
                except Exception:
                    continue

    return mapping


@lru_cache(maxsize=None)
def taxid_to_lineage_cached(taxid: int) -> Tuple[int, ...]:
    global NCBI
    if NCBI is None:
        return tuple()
    try:
        lineage = tuple(NCBI.get_lineage(int(taxid)) or [])
        return lineage
    except Exception:
        return tuple()


def lineage_to_strings(lineage: Tuple[int, ...], name_map: Dict[int, str]) -> Tuple[str, str]:
    if not lineage:
        return ("", "")
    names = [name_map.get(node_id, str(node_id)) for node_id in lineage[1:]]
    if not names:
        return ("", "")
    return (";".join(names) + ";", names[-1])


def derive_sample_id(input_file: str, sample_id_override: Optional[str]) -> str:
    if sample_id_override:
        return sample_id_override

    sample_id = os.path.basename(input_file)

    # Preserve your earlier "__sample." extraction if present
    match = re.search(r"__(.*?)\.", sample_id)
    if match:
        return match.group(1)

    # Otherwise strip common suffixes for cleaner output names
    sample_id = re.sub(r"\.(out|tsv|txt|tab)$", "", sample_id)
    return sample_id


def read_input_table(input_file: str, skip_header: int) -> pd.DataFrame:
    # Try pyarrow first for speed, then fall back
    try:
        return pd.read_csv(
            input_file,
            sep="\t",
            names=CUSTOM_FIELDS,
            header=None,
            skiprows=skip_header,
            dtype=str,
            engine="pyarrow",
        )
    except Exception:
        return pd.read_csv(
            input_file,
            sep="\t",
            names=CUSTOM_FIELDS,
            header=None,
            skiprows=skip_header,
            dtype=str,
            engine="c",
        )


def main() -> int:
    global NCBI

    args = parse_args()
    print(args.input, args.type, args.rank, args.outdir)

    if args.input is None or args.ncbidb is None or args.vhunter is None or args.type is None:
        print("Missing required args.")
        print("  Required: -i/--input, -t/--type, -n/--ncbidb, -v/--vhunter")
        return 1

    input_file = args.input
    outdir = args.outdir or "."
    sample_id = derive_sample_id(input_file, args.sample_id)

    if args.type == "mmseqs":
        out_file = os.path.join(outdir, f"{sample_id}.mmseqs.parsed")
        skip_header=1
        bx = False
    else:
        out_file = os.path.join(outdir, f"{sample_id}.blastx.parsed")
        skip_header=0
        bx = True

    try:
        connector = sqlite3.connect(args.vhunter)
        cursor = connector.cursor()
    except Exception as e:
        print(f"Cannot connect to DB: {args.vhunter} ({e})")
        return 2

    try:
        NCBI = NCBITaxa(dbfile=args.ncbidb)
    except Exception as e:
        print(f"Cannot open NCBI taxonomy DB: {args.ncbidb} ({e})")
        connector.close()
        return 2

    try:
        df = read_input_table(input_file, skip_header)
    except Exception as e:
        print(f"Failed to read input file as tabular data: {e}")
        connector.close()
        return 3

    if df.empty:
        os.makedirs(outdir, exist_ok=True)
        with open(out_file, "w") as out:
            out.write("empty diamond file")
        print(f"Wrote empty output: {out_file}")
        connector.close()
        return 0

    for c in NUMERIC_COLS:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    if "qseqid" not in df.columns or "sseqid" not in df.columns:
        print("Input file missing required columns; expected qseqid and sseqid.")
        connector.close()
        return 4

    if not args.no_sort:
        df = df.sort_values(
            ["qseqid", "sseqid", "bitscore", "evalue"],
            ascending=[True, True, False, True],
            kind="mergesort",
        )

    os.makedirs(outdir, exist_ok=True)

    print("parsing blast output files...\n")

    table = "mappings" if bx else "mappings"

    # Unique accession -> taxid
    unique_accessions = df["sseqid"].dropna().astype(str).unique().tolist()
    acc2tax = fetch_accession_taxids(cursor, unique_accessions, table, args.batch_size)
    df["txid"] = df["sseqid"].map(acc2tax)

    # Unique taxid -> lineage
    unique_taxids_series = pd.Series(df["txid"].dropna().unique())
    if len(unique_taxids_series) > 0:
        unique_taxids = unique_taxids_series.astype(int).tolist()
    else:
        unique_taxids = []

    taxid2lineage = {txid: taxid_to_lineage_cached(txid) for txid in unique_taxids}

    # Bulk translate all lineage node ids once
    all_nodes = set()
    for lineage in taxid2lineage.values():
        all_nodes.update(lineage)

    try:
        name_map = NCBI.get_taxid_translator(list(all_nodes)) if all_nodes else {}
    except Exception:
        name_map = {}

    # Precompute lineage strings once per unique taxid
    taxid2strings = {
        txid: lineage_to_strings(lineage, name_map)
        for txid, lineage in taxid2lineage.items()
    }

    mapped = df["txid"].map(taxid2strings)
    mapped_filled = [x if isinstance(x, tuple) else ("", "") for x in mapped.tolist()]
    df[["lineage", "lowest_classification"]] = pd.DataFrame(
        mapped_filled,
        index=df.index,
        columns=["lineage", "lowest_classification"],
    )

    df["alignment_method"] = "diamond" if bx else "mmseqs"
    df.insert(0, "Sample",sample_id)

    print(f"Writing output to: {out_file}")
    try:
        df.to_csv(out_file, index=False)
    except Exception as e:
        print(f"Failed while writing output file {out_file}: {e}")
        connector.close()
        return 5

    print(f"Done. Output written: {out_file}")
    connector.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
