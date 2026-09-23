#!/usr/bin/env python3
"""
Enhanced ge_screen HTML dashboard.

Improvements in this version:
  - embeds Plotly JS by default for offline/HPC use
  - horizontal scrolling for wide tables and sample-heavy plots
  - safer Plotly marker sizes
  - overview charts split into readable groups
  - top-N charts for high-cardinality metrics
  - all-sample heatmaps remain available
  - sample navigation filter
  - summary table links to sample sections
  - coverage parser supports samtools '#rname'
  - stricter mobile/engineering evidence scoring than the early prototype

Usage:
  python ge_screen_dashboard_scrollfix.py
  python ge_screen_dashboard_scrollfix.py --input-dir ge_screen_results --output-html ge_screen_dashboard.html
"""

from __future__ import annotations

import argparse
import datetime as _dt
import html
import re
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import pandas as pd

try:
    import plotly.express as px
    from plotly.offline import plot, get_plotlyjs
except ImportError as e:
    raise SystemExit(
        "This script requires pandas and plotly.\n"
        "Install with:\n"
        "  mamba install -c conda-forge pandas plotly"
    ) from e


NT_COLS = ["qseqid", "sacc", "stitle", "pident", "length", "qcovs", "evalue", "bitscore", "staxids"]
UNIVEC_COLS = ["qseqid", "sacc", "stitle", "pident", "length", "qcovs", "evalue", "bitscore"]
COVERAGE_COLS = ["rname", "startpos", "endpos", "numreads", "covbases", "coverage", "meandepth", "meanbaseq", "meanmapq"]
ABRICATE_SUMMARY_COLS = ["#FILE", "SEQUENCE", "DB", "NUM_FOUND"]
ABRICATE_HIT_COLS = [
    "#FILE", "SEQUENCE", "START", "END", "STRAND", "GENE", "COVERAGE", "COVERAGE_MAP", "GAPS",
    "%COVERAGE", "%IDENTITY", "DATABASE", "ACCESSION", "PRODUCT", "RESISTANCE"
]

AMRFINDER_COLS = [
    "Protein identifier", "Contig id", "Start", "Stop", "Strand", "Gene symbol",
    "Sequence name", "Scope", "Element type", "Element subtype", "Class", "Subclass",
    "Method", "Target length", "Reference sequence length", "% Coverage of reference sequence",
    "% Identity to reference sequence", "Alignment length", "Accession of closest sequence",
    "Name of closest sequence", "HMM id", "HMM description"
]

# Compact HTML defaults. Full BLAST and coverage TSVs can be enormous; embedding
# all rows from all samples can create multi-GB HTML files.
DEFAULT_MAX_NT_ROWS = 100
DEFAULT_MAX_UNIVEC_ROWS = 100
DEFAULT_MAX_COVERAGE_ROWS = 250
DEFAULT_MAX_ABRICATE_ROWS = 500
DEFAULT_MAX_AMRFINDER_ROWS = 500

MOBILE_PATTERNS: Dict[str, List[str]] = {
    "transposase": [
        r"\btransposase\b", r"\btnp[ab]?\b", r"\btransposable element\b",
    ],
    "insertion_sequence": [
        r"\binsertion sequence\b", r"\bIS\d+[A-Za-z0-9_-]*\b", r"\bIS[0-9A-Za-z_-]+\s+family\b",
    ],
    "integrase_recombinase": [
        r"\bintegrase\b", r"\bsite[- ]specific recombinase\b", r"\btyrosine recombinase\b",
        r"\bserine recombinase\b", r"\bresolvase\b", r"\binvertase\b", r"\bintI\d*\b",
    ],
    "conjugation_mobilization": [
        r"\brelaxase\b", r"\bmobilization protein\b", r"\bMob[A-Z0-9_-]*\b",
        r"\bconjugative transfer\b", r"\btype iv secretion\b", r"\bTra[A-Z0-9_-]*\b",
        r"\bTrb[A-Z0-9_-]*\b",
    ],
    "plasmid_rep_partition": [
        r"\breplication initiation protein\b", r"\bRep[A-Z0-9_-]*\b",
        r"\bplasmid partition\b", r"\bParA\b", r"\bParB\b", r"\bstability protein\b",
    ],
    "phage_integration_packaging": [
        r"\bprophage\b", r"\bphage integrase\b", r"\bterminase\b", r"\bportal protein\b",
        r"\bhead[- ]tail connector\b",
    ],
    "phage_structural": [
        r"\bcapsid protein\b", r"\btail fiber protein\b", r"\btail sheath protein\b",
        r"\btail tape measure protein\b", r"\bbaseplate protein\b", r"\bholin\b", r"\bendolysin\b",
    ],
}

ENGINEERING_PATTERNS: Dict[str, List[str]] = {
    "vector_backbone": [
        r"\bcloning vector\b", r"\bexpression vector\b", r"\bplasmid vector\b", r"\bbinary vector\b",
        r"\bsynthetic construct\b", r"\bmultiple cloning site\b", r"\bMCS\b",
        r"\bpUC\d*\b", r"\bpBR322\b", r"\bpET[-_ ]?\d*\b", r"\bpSMART\b", r"\bpSEVA\b",
        r"\bpBlueSTAR\b",
    ],
    "reporter_marker": [
        r"\bGFP\b", r"\bEGFP\b", r"\bRFP\b", r"\bmCherry\b", r"\bluciferase\b",
        r"\blacZ\b", r"\bbeta[- ]galactosidase alpha\b", r"\blacI\b", r"\blacY\b", r"\blacA\b",
    ],
    "promoter_recombination": [
        r"\bT7 promoter\b", r"\bSP6 promoter\b", r"\bCMV promoter\b", r"\bSV40 promoter\b",
        r"\bpromoter probe vector\b", r"\bloxP\b", r"\bFRT\b", r"\bGateway\b",
    ],
    "engineering_nuclease": [
        r"\bCas9\b", r"\bCRISPR\b", r"\bsgRNA\b", r"\bgRNA\b",
    ],
    "selectable_marker": [
        r"\bkanamycin resistance\b", r"\bneomycin resistance\b", r"\bampicillin resistance\b",
        r"\bchloramphenicol resistance\b", r"\bhygromycin resistance\b", r"\bbleomycin resistance\b",
        r"\bsacB\b", r"\bpositive selection\b",
    ],
}


def slugify(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_-]+", "_", str(value)).strip("_") or "sample"


def file_uri(path: str) -> str:
    if not path:
        return ""
    try:
        return Path(path).resolve().as_uri()
    except Exception:
        return ""


def normalize_output_columns(df: pd.DataFrame) -> pd.DataFrame:
    if df is None or df.empty:
        return df
    df = df.copy()
    rename = {}
    if "#rname" in df.columns and "rname" not in df.columns:
        rename["#rname"] = "rname"
    if "FILE" in df.columns and "#FILE" not in df.columns:
        rename["FILE"] = "#FILE"
    if rename:
        df = df.rename(columns=rename)
    return df


def safe_read_tsv(path: Path | None, cols: List[str], headerless_default: bool = False) -> pd.DataFrame:
    if path is None or not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame(columns=cols)
    try:
        first = ""
        with path.open("rt", errors="replace") as handle:
            for line in handle:
                if line.strip():
                    first = line.rstrip("\n")
                    break
        fields = first.split("\t") if first else []
        normalized_cols = [c.lstrip("#") for c in cols]
        has_header = bool(fields and (fields[0] in cols or fields[0].lstrip("#") in normalized_cols))

        if headerless_default and not has_header:
            # Infer common headerless BLAST outfmt 6 variants.
            # Legacy ge_screen NT:
            #   qseqid sacc stitle pident length qcovs evalue bitscore staxids
            # Extended supported:
            #   qseqid qlen qstart qend sacc stitle sstart send length slen
            #   pident qcovs sstrand gaps evalue bitscore score [staxids]
            if len(fields) >= 17:
                ext_cols = [
                    "qseqid", "qlen", "qstart", "qend", "sacc", "stitle",
                    "sstart", "send", "length", "slen", "pident", "qcovs",
                    "sstrand", "gaps", "evalue", "bitscore", "score"
                ]
                if len(fields) >= 18:
                    ext_cols.append("staxids")
                df = pd.read_csv(path, sep="\t", header=None, names=ext_cols[:len(fields)], dtype=str)
            elif len(fields) == 10 and len(cols) == 9:
                ext_cols = ["qseqid", "sacc", "stitle", "pident", "length", "slen", "qcovs", "evalue", "bitscore", "staxids"]
                df = pd.read_csv(path, sep="\t", header=None, names=ext_cols, dtype=str)
            elif len(fields) == 9 and len(cols) == 8:
                ext_cols = ["qseqid", "sacc", "stitle", "pident", "length", "slen", "qcovs", "evalue", "bitscore"]
                df = pd.read_csv(path, sep="\t", header=None, names=ext_cols, dtype=str)
            else:
                df = pd.read_csv(path, sep="\t", header=None, names=cols, dtype=str)
        else:
            df = pd.read_csv(path, sep="\t", dtype=str)
            if not has_header and df.shape[1] == len(cols):
                df = pd.read_csv(path, sep="\t", header=None, names=cols, dtype=str)

        if df.shape[1] == len(cols) and headerless_default and list(df.columns) != cols:
            df.columns = cols
        return normalize_output_columns(df)
    except pd.errors.EmptyDataError:
        return pd.DataFrame(columns=cols)
    except Exception:
        try:
            return normalize_output_columns(pd.read_csv(path, sep="\t", header=None, names=cols, dtype=str))
        except Exception:
            return pd.DataFrame(columns=cols)


def coerce_numeric(df: pd.DataFrame, cols: Iterable[str]) -> pd.DataFrame:
    df = df.copy()
    for col in cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def add_subject_coverage(df: pd.DataFrame) -> pd.DataFrame:
    """Add scovs, the percentage of the BLAST subject covered by the alignment.

    Requires subject length (slen). If sstart/send are present, use the subject
    span: abs(send - sstart) + 1. Otherwise fall back to alignment length/slen.
    """
    if df is None or df.empty:
        return df

    df = df.copy()

    # Normalize a few common subject-length aliases.
    if "slen" not in df.columns:
        for alt in ["subject_length", "s_length", "sseq_len", "sseq_length"]:
            if alt in df.columns:
                df = df.rename(columns={alt: "slen"})
                break

    numeric_cols = [c for c in ["sstart", "send", "slen", "length"] if c in df.columns]
    df = coerce_numeric(df, numeric_cols)

    if "scovs" not in df.columns:
        df["scovs"] = pd.NA

    if "slen" not in df.columns:
        return df

    valid_slen = df["slen"].notna() & (df["slen"] > 0)

    if {"sstart", "send"}.issubset(df.columns):
        span = (df["send"] - df["sstart"]).abs() + 1
        mask = valid_slen & span.notna()
        df.loc[mask, "scovs"] = (100.0 * span[mask] / df.loc[mask, "slen"]).clip(upper=100.0)

    if "length" in df.columns:
        missing = df["scovs"].isna()
        mask = missing & valid_slen & df["length"].notna()
        df.loc[mask, "scovs"] = (100.0 * df.loc[mask, "length"] / df.loc[mask, "slen"]).clip(upper=100.0)

    return df


def ensure_numeric_columns(df: pd.DataFrame, cols: Iterable[str], default: float = 0.0) -> pd.DataFrame:
    df = df.copy()
    for col in cols:
        if col not in df.columns:
            df[col] = default
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(default)
        df[col] = df[col].replace([float("inf"), float("-inf")], default)
    return df


def positive_size(series: pd.Series, minimum: float = 1.0) -> pd.Series:
    return pd.to_numeric(series, errors="coerce").fillna(0).clip(lower=minimum)


def find_first(directory: Path, patterns: Iterable[str]) -> Path | None:
    for patt in patterns:
        hits = sorted(directory.glob(patt))
        if hits:
            return hits[0]
    return None


def has_any_outputs(directory: Path) -> bool:
    return any(any(directory.glob(p)) for p in ["*_vs_nt.tsv", "*_vs_univec.tsv", "abricate.summary.tsv", "*.coverage.tsv"])


def collect_sample_dirs(root: Path) -> List[Path]:
    dirs = []
    if has_any_outputs(root):
        dirs.append(root)
    for d in sorted(root.rglob("*")):
        if d.is_dir() and has_any_outputs(d):
            dirs.append(d)
    seen = set()
    out = []
    for d in dirs:
        key = str(d.resolve())
        if key not in seen:
            seen.add(key)
            out.append(d)
    return out


def pattern_categories(texts: Iterable[str], pattern_map: Dict[str, List[str]]) -> Tuple[Dict[str, int], int]:
    compiled = {k: [re.compile(p, re.IGNORECASE) for p in pats] for k, pats in pattern_map.items()}
    counts = {k: 0 for k in pattern_map}
    total = 0
    for text in texts:
        value = str(text or "")
        for cat, pats in compiled.items():
            if any(p.search(value) for p in pats):
                counts[cat] += 1
                total += 1
    counts = {k: v for k, v in counts.items() if v > 0}
    return counts, total


def evidence_texts_from_df(df: pd.DataFrame, columns: Iterable[str]) -> List[str]:
    if df is None or df.empty:
        return []
    out = []
    for col in columns:
        if col in df.columns:
            out.extend(df[col].dropna().astype(str).tolist())
    return out


def filtered_blast_texts(df: pd.DataFrame, text_col: str, max_rows: int, min_qcovs: float, min_pident: float) -> List[str]:
    if df is None or df.empty or text_col not in df.columns:
        return []
    work = coerce_numeric(df.copy(), ["qcovs", "pident", "bitscore"])
    if "bitscore" in work.columns:
        work = work.sort_values("bitscore", ascending=False, na_position="last")
    if {"qcovs", "pident"}.issubset(work.columns):
        work = work[(work["qcovs"].fillna(0) >= min_qcovs) & (work["pident"].fillna(0) >= min_pident)]
    return work.head(max_rows)[text_col].dropna().astype(str).tolist()


def classify_sample(ov: Dict) -> Tuple[str, List[str]]:
    reasons = []
    red = 0
    yellow = 0

    if ov["strong_univec_hits"] > 0:
        red += 2
        reasons.append(f"{ov['strong_univec_hits']} high-confidence UniVec hit(s)")

    if ov["engineering_evidence_count"] > 0:
        red += 2
        reasons.append(f"specific engineering/vector evidence: {ov['engineering_evidence_categories']}")

    if ov["univec_hits"] > 0 and ov["strong_univec_hits"] == 0:
        yellow += 1
        reasons.append("weak/low-coverage UniVec hit(s); review manually")

    if ov["plasmidfinder_hits"] > 0:
        yellow += 1
        reasons.append("PlasmidFinder hit(s) detected")

    mobile_cat_count = len([x for x in ov["mobile_evidence_categories"].split(", ") if x])
    if mobile_cat_count >= 2:
        yellow += 1
        reasons.append(f"specific mobile-element evidence across {mobile_cat_count} categories: {ov['mobile_evidence_categories']}")

    if ov["has_coverage"] and ov["relative_depth_max"] >= 5:
        yellow += 1
        reasons.append(f"high relative contig depth; max/median_nonzero_depth={ov['relative_depth_max']:.1f}")

    if not ov["has_coverage"]:
        reasons.append("coverage/depth files missing or empty; depth-based interpretation unavailable")

    if red >= 2:
        return "Red", reasons
    if yellow >= 1:
        return "Yellow", reasons
    return "Green", reasons or ["no configured vector/mobile-element flags detected"]


def summarize_sample(sample_dir: Path, limits: Dict[str, int | None]) -> Dict:
    sample = sample_dir.name
    nt_path = find_first(sample_dir, ["*_vs_nt.tsv"])
    univec_path = find_first(sample_dir, ["*_vs_univec.tsv"])
    cov_path = find_first(sample_dir, ["*.coverage.tsv"])
    depth_path = find_first(sample_dir, ["*.depth.tsv"])
    amrfinder_path = find_first(sample_dir, ["*.amrfinderplus.tsv", "*amrfinder*.tsv"])

    nt = safe_read_tsv(nt_path, NT_COLS, headerless_default=True)
    nt = add_subject_coverage(nt)
    nt = coerce_numeric(nt, ["pident", "length", "qcovs", "scovs", "slen", "sstart", "send", "evalue", "bitscore"])
    univec = safe_read_tsv(univec_path, UNIVEC_COLS, headerless_default=True)
    univec = add_subject_coverage(univec)
    univec = coerce_numeric(univec, ["pident", "length", "qcovs", "scovs", "slen", "sstart", "send", "evalue", "bitscore"])
    cov = coerce_numeric(normalize_output_columns(safe_read_tsv(cov_path, COVERAGE_COLS)), ["numreads", "covbases", "coverage", "meandepth", "meanbaseq", "meanmapq"])

    abr_summary = coerce_numeric(safe_read_tsv(sample_dir / "abricate.summary.tsv", ABRICATE_SUMMARY_COLS), ["NUM_FOUND"])

    abr_hits = {}
    for db in ["plasmidfinder", "card", "vfdb"]:
        abr_hits[db] = coerce_numeric(
            safe_read_tsv(sample_dir / f"abricate.{db}.tsv", ABRICATE_HIT_COLS),
            ["%COVERAGE", "%IDENTITY"],
        )

    amrfinder = safe_read_tsv(amrfinder_path, AMRFINDER_COLS, headerless_default=False)
    amrfinder = coerce_numeric(
        amrfinder,
        ["Start", "Stop", "Target length", "Reference sequence length",
         "% Coverage of reference sequence", "% Identity to reference sequence",
         "Alignment length"],
    )

    nt_top = ""
    if not nt.empty and "stitle" in nt.columns:
        nt_top = str(nt.sort_values(["bitscore", "qcovs", "pident"], ascending=[False, False, False], na_position="last").iloc[0].get("stitle", ""))

    engineering_texts = []
    engineering_texts += filtered_blast_texts(univec, "stitle", max_rows=200, min_qcovs=20, min_pident=80)
    engineering_texts += filtered_blast_texts(nt, "stitle", max_rows=50, min_qcovs=70, min_pident=95)
    eng_cats, eng_count = pattern_categories(engineering_texts, ENGINEERING_PATTERNS)

    mobile_texts = []
    for db in ["plasmidfinder", "card", "vfdb"]:
        mobile_texts += evidence_texts_from_df(abr_hits[db], ["GENE", "PRODUCT", "RESISTANCE"])
    mobile_texts += evidence_texts_from_df(amrfinder, ["Gene symbol", "Sequence name", "Element subtype", "Class", "Subclass", "HMM description"])
    mobile_texts += filtered_blast_texts(nt, "stitle", max_rows=75, min_qcovs=70, min_pident=90)
    mob_cats, mob_count = pattern_categories(mobile_texts, MOBILE_PATTERNS)

    has_coverage = not cov.empty and "rname" in cov.columns and "meandepth" in cov.columns
    median_nonzero_depth = 0.0
    relative_depth_max = 0.0
    if has_coverage:
        depths = pd.to_numeric(cov["meandepth"], errors="coerce").fillna(0)
        nonzero = depths[depths > 0]
        if len(nonzero):
            median_nonzero_depth = float(nonzero.median())
            relative_depth_max = float(depths.max() / max(median_nonzero_depth, 1e-9))

    strong_univec_hits = 0
    if not univec.empty and {"qcovs", "pident"}.issubset(univec.columns):
        strong_univec_hits = int(((univec["qcovs"].fillna(0) >= 60) & (univec["pident"].fillna(0) >= 95)).sum())

    overview = {
        "sample": sample,
        "sample_dir": str(sample_dir),
        "sample_anchor": slugify(sample),
        "nt_hits": int(len(nt)),
        "nt_top_hit": nt_top,
        "univec_hits": int(len(univec)),
        "strong_univec_hits": strong_univec_hits,
        "plasmidfinder_hits": int(len(abr_hits["plasmidfinder"])),
        "card_hits": int(len(abr_hits["card"])),
        "vfdb_hits": int(len(abr_hits["vfdb"])),
        "amrfinderplus_hits": int(len(amrfinder)),
        "amrfinderplus_amr_hits": int((amrfinder["Element type"].fillna("").astype(str).str.upper() == "AMR").sum()) if "Element type" in amrfinder.columns and not amrfinder.empty else 0,
        "amrfinderplus_stress_hits": int((amrfinder["Element type"].fillna("").astype(str).str.lower().str.contains("stress")).sum()) if "Element type" in amrfinder.columns and not amrfinder.empty else 0,
        "amrfinderplus_virulence_hits": int((amrfinder["Element type"].fillna("").astype(str).str.lower().str.contains("virulence")).sum()) if "Element type" in amrfinder.columns and not amrfinder.empty else 0,
        "contigs_in_coverage": int(len(cov)),
        "mean_depth_avg": float(cov["meandepth"].mean()) if has_coverage else 0.0,
        "max_depth": float(cov["meandepth"].max()) if has_coverage else 0.0,
        "median_nonzero_depth": median_nonzero_depth,
        "relative_depth_max": relative_depth_max,
        "mean_coverage_pct": float(cov["coverage"].mean()) if has_coverage and "coverage" in cov.columns else 0.0,
        "max_univec_qcovs": float(univec["qcovs"].max()) if not univec.empty and "qcovs" in univec.columns else 0.0,
        "max_univec_scovs": float(univec["scovs"].max()) if not univec.empty and "scovs" in univec.columns and univec["scovs"].notna().any() else 0.0,
        "max_univec_pident": float(univec["pident"].max()) if not univec.empty and "pident" in univec.columns else 0.0,
        "engineering_evidence_count": eng_count,
        "engineering_evidence_categories": ", ".join(sorted(eng_cats.keys())),
        "mobile_evidence_count": mob_count,
        "mobile_evidence_categories": ", ".join(sorted(mob_cats.keys())),
        "engineering_keyword_hits": eng_count,
        "mobile_keyword_hits": mob_count,
        "has_nt": bool(nt_path and nt_path.exists()),
        "has_univec": bool(univec_path and univec_path.exists()),
        "has_coverage": bool(has_coverage),
        "has_depth": bool(depth_path and depth_path.exists()),
        "has_abricate_summary": bool((sample_dir / "abricate.summary.tsv").exists()),
        "has_plasmidfinder": bool((sample_dir / "abricate.plasmidfinder.tsv").exists()),
        "has_card": bool((sample_dir / "abricate.card.tsv").exists()),
        "has_vfdb": bool((sample_dir / "abricate.vfdb.tsv").exists()),
        "has_amrfinderplus": bool(amrfinder_path and amrfinder_path.exists()),
        "embedded_nt_rows": int(min(len(nt), limits["nt"])) if limits["nt"] is not None else int(len(nt)),
        "embedded_univec_rows": int(min(len(univec), limits["univec"])) if limits["univec"] is not None else int(len(univec)),
        "embedded_coverage_rows": int(min(len(cov), limits["coverage"])) if limits["coverage"] is not None else int(len(cov)),
        "nt_file": str(nt_path) if nt_path else "",
        "univec_file": str(univec_path) if univec_path else "",
        "coverage_file": str(cov_path) if cov_path else "",
        "depth_file": str(depth_path) if depth_path else "",
        "amrfinderplus_file": str(amrfinder_path) if amrfinder_path else "",
    }
    triage, reasons = classify_sample(overview)
    overview["triage_flag"] = triage
    overview["triage_reasons"] = "; ".join(reasons)

    nt_sorted = nt.sort_values(["bitscore", "qcovs", "pident"], ascending=[False, False, False], na_position="last") if not nt.empty else pd.DataFrame(columns=NT_COLS)
    univec_sorted = univec.sort_values(["bitscore", "qcovs", "pident"], ascending=[False, False, False], na_position="last") if not univec.empty else pd.DataFrame(columns=UNIVEC_COLS)
    cov_sorted = cov.sort_values("meandepth", ascending=False, na_position="last") if has_coverage else pd.DataFrame(columns=COVERAGE_COLS)

    tables = {
        "nt": nt_sorted.head(limits["nt"]) if limits["nt"] is not None else nt_sorted,
        "univec": univec_sorted.head(limits["univec"]) if limits["univec"] is not None else univec_sorted,
        "coverage": cov_sorted.head(limits["coverage"]) if limits["coverage"] is not None else cov_sorted,
        "abricate_summary": abr_summary,
        "amrfinderplus": amrfinder.head(limits["amrfinder"]) if limits.get("amrfinder") is not None else amrfinder,
        "abricate_top": {},
    }
    for db, df in abr_hits.items():
        if not df.empty:
            sort_cols = [c for c in ["%IDENTITY", "%COVERAGE"] if c in df.columns]
            abr_sorted = df.sort_values(sort_cols, ascending=[False] * len(sort_cols), na_position="last") if sort_cols else df
            tables["abricate_top"][db] = abr_sorted.head(limits["abricate"]) if limits["abricate"] is not None else abr_sorted
        else:
            tables["abricate_top"][db] = pd.DataFrame(columns=ABRICATE_HIT_COLS)
    return {"overview": overview, "tables": tables}


TABLE_COUNTER = 0


def df_to_html_table(df: pd.DataFrame, max_rows: int | None = None, escape_cells: bool = False) -> str:
    """Render a dataframe as a sortable, searchable, paginated HTML table.

    max_rows=None means include all rows in the HTML and paginate client-side.
    """
    global TABLE_COUNTER

    if df is None or df.empty:
        return '<p class="empty">No data available.</p>'

    TABLE_COUNTER += 1
    table_id = f"tbl_{TABLE_COUNTER}"

    x = df.copy()
    if max_rows is not None:
        x = x.head(max_rows)

    for col in x.columns:
        if pd.api.types.is_numeric_dtype(x[col]):
            x[col] = x[col].map(lambda v: "" if pd.isna(v) else f"{v:.3f}" if isinstance(v, float) else str(v))

    table = x.to_html(
        index=False,
        classes="data-table sortable-table",
        table_id=table_id,
        border=0,
        escape=escape_cells,
    )

    return f"""
<div class="table-widget" data-table-id="{table_id}">
  <div class="table-controls">
    <label>Search: <input class="table-search" type="search" placeholder="Filter table..." oninput="filterDashboardTable('{table_id}')"></label>
    <label>Rows:
      <select class="page-size" onchange="setDashboardPageSize('{table_id}', this.value)">
        <option value="10">10</option>
        <option value="25" selected>25</option>
        <option value="50">50</option>
        <option value="100">100</option>
        <option value="all">All</option>
      </select>
    </label>
    <button type="button" onclick="prevDashboardPage('{table_id}')">Prev</button>
    <span class="page-status" id="{table_id}_status">Page 1</span>
    <button type="button" onclick="nextDashboardPage('{table_id}')">Next</button>
  </div>
  <div class="table-scroll">{table}</div>
</div>
"""


def fig_div(fig, min_width_px: int | None = None) -> str:
    div = plot(fig, output_type="div", include_plotlyjs=False, config={"responsive": True})
    if min_width_px:
        return f'<div class="plot-scroll"><div class="plot-wide" style="min-width:{min_width_px}px">{div}</div></div>'
    return f'<div class="plot-scroll">{div}</div>'


def warn_panel(title: str, body: str) -> str:
    return f'<div class="warn-panel"><strong>{html.escape(title)}</strong><br>{html.escape(body)}</div>'


def make_flag_panel(ov: Dict) -> str:
    cls = {"Red": "flag-red", "Yellow": "flag-yellow", "Green": "flag-green"}.get(ov["triage_flag"], "flag-green")
    reasons = [r for r in str(ov.get("triage_reasons", "")).split("; ") if r]
    items = "".join(f"<li>{html.escape(r)}</li>" for r in reasons)
    return f"""
<div class="flag-panel {cls}">
  <div class="flag-title">Triage flag: {html.escape(str(ov['triage_flag']))}</div>
  <ul>{items}</ul>
</div>
"""


def source_links(ov: Dict) -> str:
    entries = []
    for label, key in [("Sample directory", "sample_dir"), ("NT TSV", "nt_file"), ("UniVec TSV", "univec_file"), ("Coverage TSV", "coverage_file"), ("Depth TSV", "depth_file"), ("AMRFinderPlus TSV", "amrfinderplus_file")]:
        value = ov.get(key, "")
        if not value:
            continue
        uri = file_uri(value)
        if uri:
            entries.append(f'<li><strong>{html.escape(label)}:</strong> <a href="{html.escape(uri)}"><code>{html.escape(value)}</code></a></li>')
        else:
            entries.append(f'<li><strong>{html.escape(label)}:</strong> <code>{html.escape(value)}</code></li>')
    return "<ul>" + "".join(entries) + "</ul>" if entries else '<p class="empty">No source files captured.</p>'


def overview_table(overview_df: pd.DataFrame) -> str:
    cols = [
        "sample", "triage_flag", "triage_reasons", "nt_hits", "univec_hits", "strong_univec_hits",
        "plasmidfinder_hits", "card_hits", "vfdb_hits", "amrfinderplus_hits", "amrfinderplus_amr_hits", "amrfinderplus_stress_hits", "amrfinderplus_virulence_hits", "engineering_evidence_count",
        "engineering_evidence_categories", "mobile_evidence_count", "mobile_evidence_categories",
        "has_coverage", "contigs_in_coverage", "embedded_nt_rows", "embedded_univec_rows", "embedded_coverage_rows",
        "mean_depth_avg", "max_depth", "relative_depth_max", "max_univec_qcovs", "max_univec_scovs", "max_univec_pident", "nt_top_hit",
    ]
    df = overview_df[[c for c in cols if c in overview_df.columns]].copy()
    if "sample" in df.columns:
        df["sample"] = df.apply(lambda r: f'<a href="#sample-{html.escape(str(r.get("sample_anchor", slugify(r["sample"]))))}">{html.escape(str(r["sample"]))}</a>', axis=1)
    return df_to_html_table(df, max_rows=None, escape_cells=False)


def bar_chart(df: pd.DataFrame, x: str, y, title: str, height: int = 520, orientation: str = "v", min_width: int | None = None):
    fig = px.bar(df, x=x, y=y, title=title, barmode="group")
    fig.update_layout(height=height, xaxis_title=x, yaxis_title="Value" if isinstance(y, list) else y)
    return fig_div(fig, min_width_px=min_width)


def make_dashboard(results: List[Dict], root_dir: Path, title: str, embed_plotly: bool = True) -> Tuple[str, pd.DataFrame]:
    overview_df = pd.DataFrame([r["overview"] for r in results]).sort_values("sample")
    numeric_cols = [
        "nt_hits", "univec_hits", "strong_univec_hits", "plasmidfinder_hits", "card_hits", "vfdb_hits", "amrfinderplus_hits", "amrfinderplus_amr_hits", "amrfinderplus_stress_hits", "amrfinderplus_virulence_hits",
        "engineering_evidence_count", "mobile_evidence_count", "engineering_keyword_hits", "mobile_keyword_hits",
        "contigs_in_coverage", "mean_depth_avg", "max_depth", "median_nonzero_depth", "relative_depth_max",
        "mean_coverage_pct", "max_univec_qcovs", "max_univec_scovs", "max_univec_pident",
    ]
    overview_df = ensure_numeric_columns(overview_df, numeric_cols)
    overview_df["plot_univec_size"] = positive_size(overview_df["univec_hits"])

    sample_count = len(overview_df)
    chart_min_width = max(1200, sample_count * 42)

    flag_counts = overview_df["triage_flag"].value_counts().to_dict()
    cards = "".join(
        f'<div class="card"><div class="label">{label}</div><div class="value">{value}</div></div>'
        for label, value in [
            ("Samples", sample_count),
            ("Red", flag_counts.get("Red", 0)),
            ("Yellow", flag_counts.get("Yellow", 0)),
            ("Green", flag_counts.get("Green", 0)),
            ("With coverage", int(overview_df["has_coverage"].sum())),
            ("Without coverage", int((~overview_df["has_coverage"].astype(bool)).sum())),
        ]
    )

    charts = []

    # 1. Triage count chart.
    triage_df = overview_df["triage_flag"].value_counts().rename_axis("triage_flag").reset_index(name="count")
    fig_triage = px.bar(triage_df, x="triage_flag", y="count", title="Triage flag counts")
    fig_triage.update_layout(height=360, xaxis_title="Triage", yaxis_title="Samples")
    charts.append(("Triage flag counts", fig_div(fig_triage)))

    # 2. Top-N evidence chart avoids a 70-sample unreadable all-metric bar plot.
    evidence_cols = ["strong_univec_hits", "plasmidfinder_hits", "amrfinderplus_amr_hits", "amrfinderplus_stress_hits", "amrfinderplus_virulence_hits", "engineering_evidence_count", "mobile_evidence_count", "card_hits", "vfdb_hits"]
    top_evidence = overview_df.copy()
    top_evidence["total_evidence"] = top_evidence[evidence_cols].sum(axis=1)
    top_evidence = top_evidence.sort_values("total_evidence", ascending=False).head(30)
    evidence_long = top_evidence[["sample"] + evidence_cols].melt(id_vars="sample", var_name="metric", value_name="count")
    fig_evidence = px.bar(evidence_long, x="sample", y="count", color="metric", barmode="group", title="Top 30 samples by evidence count")
    fig_evidence.update_layout(height=650, xaxis_title="Sample", yaxis_title="Count")
    charts.append(("Top 30 samples by evidence count", fig_div(fig_evidence, min_width_px=max(1400, len(top_evidence) * 70))))

    # 3. All-sample evidence heatmap.
    heat = overview_df.set_index("sample")[evidence_cols].astype(float)
    fig_heat = px.imshow(heat.T, labels=dict(x="Sample", y="Metric", color="Count"), title="All-sample evidence heatmap", aspect="auto")
    fig_heat.update_layout(height=480)
    charts.append(("All-sample evidence heatmap", fig_div(fig_heat, min_width_px=chart_min_width)))

    # 4. Coverage availability heatmap.
    availability_cols = ["has_nt", "has_univec", "has_coverage", "has_depth", "has_abricate_summary", "has_plasmidfinder", "has_card", "has_vfdb", "has_amrfinderplus"]
    for c in availability_cols:
        if c not in overview_df.columns:
            overview_df[c] = False
    avail = overview_df.set_index("sample")[availability_cols].astype(int)
    fig_avail = px.imshow(avail.T, labels=dict(x="Sample", y="Output type", color="Present"), title="Output/data availability heatmap", aspect="auto")
    fig_avail.update_layout(height=500)
    charts.append(("Output/data availability heatmap", fig_div(fig_avail, min_width_px=chart_min_width)))

    # 5. Depth charts.
    if overview_df["has_coverage"].astype(bool).any():
        depth_cols = ["mean_depth_avg", "max_depth", "relative_depth_max"]
        depth_top = overview_df.sort_values("relative_depth_max", ascending=False).head(30)
        depth_long = depth_top[["sample"] + depth_cols].melt(id_vars="sample", var_name="metric", value_name="value")
        fig_depth_top = px.bar(depth_long, x="sample", y="value", color="metric", barmode="group", title="Top 30 samples by relative depth")
        fig_depth_top.update_layout(height=620, xaxis_title="Sample", yaxis_title="Depth / ratio")
        charts.append(("Top 30 samples by relative depth", fig_div(fig_depth_top, min_width_px=max(1400, len(depth_top) * 70))))

        fig_cov_scatter = px.scatter(
            overview_df,
            x="mean_coverage_pct",
            y="mean_depth_avg",
            size="plot_univec_size",
            color="triage_flag",
            hover_data=["sample", "relative_depth_max", "triage_reasons"],
            title="Coverage vs depth, colored by triage flag",
        )
        fig_cov_scatter.update_layout(height=540, xaxis_title="Average coverage %", yaxis_title="Average mean depth")
        charts.append(("Coverage/depth triage scatter", fig_div(fig_cov_scatter)))
    else:
        charts.append(("Coverage/depth plots", warn_panel("Coverage unavailable", "No non-empty *.coverage.tsv files were found.")))

    # 6. UniVec scatter with safe marker sizes.
    fig_uv = px.scatter(
        overview_df,
        x="max_univec_qcovs",
        y="max_univec_pident",
        size="plot_univec_size",
        color="triage_flag",
        hover_data=["sample", "univec_hits", "strong_univec_hits", "triage_reasons"],
        title="Maximum UniVec identity vs query coverage",
    )
    fig_uv.update_layout(height=520, xaxis_title="Max UniVec query coverage %", yaxis_title="Max UniVec identity %")
    charts.append(("UniVec identity vs query coverage", fig_div(fig_uv)))

    if "max_univec_scovs" in overview_df.columns and overview_df["max_univec_scovs"].fillna(0).max() > 0:
        fig_uvs = px.scatter(
            overview_df,
            x="max_univec_scovs",
            y="max_univec_pident",
            size="plot_univec_size",
            color="triage_flag",
            hover_data=["sample", "univec_hits", "strong_univec_hits", "max_univec_qcovs", "triage_reasons"],
            title="Maximum UniVec identity vs subject coverage",
        )
        fig_uvs.update_layout(height=520, xaxis_title="Max UniVec subject coverage %", yaxis_title="Max UniVec identity %")
        charts.append(("UniVec identity vs subject coverage", fig_div(fig_uvs)))

    sidebar_items = []
    sample_sections = []
    for result in results:
        ov = result["overview"]
        anchor = ov["sample_anchor"]
        sidebar_items.append(f'<li class="sample-nav-item" data-sample="{html.escape(ov["sample"].lower())}"><a href="#sample-{html.escape(anchor)}">{html.escape(ov["sample"])}</a></li>')

        sample_cards = "".join(
            f'<div class="card"><div class="label">{html.escape(label)}</div><div class="value">{html.escape(str(value))}</div></div>'
            for label, value in [
                ("NT hits", ov["nt_hits"]),
                ("UniVec hits", ov["univec_hits"]),
                ("Strong UniVec", ov["strong_univec_hits"]),
                ("PlasmidFinder", ov["plasmidfinder_hits"]),
                ("AMRFinderPlus", ov["amrfinderplus_hits"]),
                ("AMRFinder AMR", ov["amrfinderplus_amr_hits"]),
                ("Engineering evidence", ov["engineering_evidence_count"]),
                ("Mobile evidence", ov["mobile_evidence_count"]),
                ("Coverage contigs", ov["contigs_in_coverage"]),
                ("Max/median depth", f"{ov['relative_depth_max']:.2f}" if ov["has_coverage"] else "n/a"),
            ]
        )

        cov_plot = warn_panel("No coverage data", "No non-empty coverage TSV was found for this sample.")
        cov = result["tables"]["coverage"]
        if ov["has_coverage"] and cov is not None and not cov.empty and "rname" in cov.columns:
            top_cov = ensure_numeric_columns(cov.head(100).copy(), ["meandepth", "coverage", "numreads", "covbases", "meanmapq"])
            fig_cov = px.bar(
                top_cov,
                x="rname",
                y="meandepth",
                hover_data=[c for c in ["coverage", "numreads", "covbases", "meanmapq"] if c in top_cov.columns],
                title=f"{ov['sample']}: top contig mean depth",
            )
            fig_cov.update_layout(height=520, xaxis_title="Contig", yaxis_title="Mean depth")
            cov_plot = fig_div(fig_cov, min_width_px=max(1200, min(100, len(top_cov)) * 36))

            fig_cov2 = px.scatter(
                top_cov,
                x="coverage",
                y="meandepth",
                size=positive_size(top_cov["numreads"]),
                hover_data=[c for c in ["rname", "numreads", "covbases", "meanmapq"] if c in top_cov.columns],
                title=f"{ov['sample']}: contig coverage vs depth",
            )
            fig_cov2.update_layout(height=460, xaxis_title="Coverage %", yaxis_title="Mean depth")
            cov_plot += fig_div(fig_cov2)

        univec_plot = '<p class="empty">No UniVec hits available.</p>'
        univec = result["tables"]["univec"]
        if univec is not None and not univec.empty and {"qcovs", "pident"}.issubset(univec.columns):
            univec = ensure_numeric_columns(univec.copy(), ["qcovs", "pident", "bitscore"])
            univec["plot_bitscore"] = positive_size(univec["bitscore"])
            fig = px.scatter(
                univec,
                x="qcovs",
                y="pident",
                size="plot_bitscore",
                hover_data=[c for c in ["qseqid", "sacc", "stitle", "length", "slen", "qcovs", "scovs", "evalue"] if c in univec.columns],
                title=f"{ov['sample']}: UniVec identity vs query coverage",
            )
            fig.update_layout(height=460, xaxis_title="Query coverage %", yaxis_title="Percent identity")
            univec_plot = fig_div(fig)

        nt_plot = '<p class="empty">No NT hits available.</p>'
        nt = result["tables"]["nt"]
        if nt is not None and not nt.empty and {"qcovs", "pident"}.issubset(nt.columns):
            nt = ensure_numeric_columns(nt.copy(), ["qcovs", "pident", "bitscore"])
            nt["plot_bitscore"] = positive_size(nt["bitscore"])
            fig = px.scatter(
                nt,
                x="qcovs",
                y="pident",
                size="plot_bitscore",
                hover_data=[c for c in ["qseqid", "sacc", "stitle", "length", "slen", "qcovs", "scovs", "evalue"] if c in nt.columns],
                title=f"{ov['sample']}: top NT hit identity vs query coverage",
            )
            fig.update_layout(height=460, xaxis_title="Query coverage %", yaxis_title="Percent identity")
            nt_plot = fig_div(fig)

        sample_sections.append(f"""
<section id="sample-{html.escape(anchor)}" class="sample-section">
  <div class="sample-header">
    <h2>{html.escape(ov['sample'])}</h2>
    <a class="back-link" href="#overview">Back to overview</a>
  </div>
  <p><strong>Directory:</strong> <code>{html.escape(ov['sample_dir'])}</code></p>
  {make_flag_panel(ov)}
  <p><strong>Top NT hit:</strong> {html.escape(ov['nt_top_hit']) if ov['nt_top_hit'] else '<span class="empty">No NT hits</span>'}</p>
  <div class="cards">{sample_cards}</div>
  <details open><summary>Coverage/depth plots</summary>{cov_plot}</details>
  <details><summary>UniVec identity vs coverage plot</summary>{univec_plot}</details>
  <details><summary>NT identity vs coverage plot</summary>{nt_plot}</details>
  <details><summary>Embedded table limits / full source files</summary>
    <p>This compact dashboard embeds capped table slices to prevent multi-GB HTML files. Use source TSV links below for complete raw output tables.</p>
    {source_links(ov)}
  </details>
  <details><summary>Top NT hits</summary>{df_to_html_table(result['tables']['nt'], None)}</details>
  <details><summary>Top UniVec hits</summary>{df_to_html_table(result['tables']['univec'], None)}</details>
  <details><summary>Coverage summary</summary>{df_to_html_table(result['tables']['coverage'], None)}</details>
  <details><summary>AMRFinderPlus hits</summary>{df_to_html_table(result['tables']['amrfinderplus'], None)}</details>
  <details><summary>Abricate summary</summary>{df_to_html_table(result['tables']['abricate_summary'], None)}</details>
  <details><summary>Abricate top hits: plasmidfinder</summary>{df_to_html_table(result['tables']['abricate_top']['plasmidfinder'], None)}</details>
  <details><summary>Abricate top hits: card</summary>{df_to_html_table(result['tables']['abricate_top']['card'], None)}</details>
  <details><summary>Abricate top hits: vfdb</summary>{df_to_html_table(result['tables']['abricate_top']['vfdb'], None)}</details>
</section>
""")

    plotly_tag = f"<script>{get_plotlyjs()}</script>" if embed_plotly else '<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>'
    charts_html = "\n".join(f"<details open><summary>{html.escape(name)}</summary>{body}</details>" for name, body in charts)
    generated_at = _dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    html_doc = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8" />
<meta name="viewport" content="width=device-width, initial-scale=1" />
<title>{html.escape(title)}</title>
{plotly_tag}
<style>
  :root {{
    --sidebar-width: 330px;
    --card-bg: #ffffff;
    --page-bg: #f5f7fb;
    --text: #1f2937;
  }}
  * {{ box-sizing: border-box; }}
  html {{ scroll-behavior: smooth; }}
  body {{ font-family: Arial, Helvetica, sans-serif; margin: 0; background: var(--page-bg); color: var(--text); overflow-x: auto; }}
  .layout {{ display: grid; grid-template-columns: var(--sidebar-width) minmax(0, 1fr); min-height: 100vh; }}
  aside {{ background: #111827; color: #fff; padding: 20px; position: sticky; top: 0; height: 100vh; overflow-y: auto; overflow-x: hidden; }}
  aside h1 {{ font-size: 1.08rem; margin-top: 0; line-height: 1.25; }}
  aside a {{ color: #bfdbfe; text-decoration: none; }}
  aside ul {{ padding-left: 18px; }}
  .sample-nav-item {{ overflow-wrap: anywhere; margin-bottom: 4px; }}
  #sampleFilter {{ width: 100%; padding: 8px; border-radius: 8px; border: 1px solid #374151; margin: 8px 0 12px; }}
  main {{ padding: 24px; min-width: 0; max-width: none; }}
  section {{ background: var(--card-bg); border-radius: 12px; padding: 18px; margin-bottom: 20px; box-shadow: 0 1px 3px rgba(0,0,0,0.08); min-width: 0; }}
  .sample-header {{ display: flex; justify-content: space-between; align-items: baseline; gap: 10px; }}
  .sample-header h2 {{ overflow-wrap: anywhere; }}
  .back-link {{ font-size: 0.9rem; white-space: nowrap; }}
  .cards {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(165px, 1fr)); gap: 10px; margin: 14px 0 18px 0; }}
  .card {{ background: #eff6ff; border: 1px solid #dbeafe; border-radius: 10px; padding: 12px; min-width: 0; }}
  .label {{ font-size: 0.85rem; color: #374151; }}
  .value {{ font-size: 1.25rem; font-weight: 700; margin-top: 4px; word-break: break-word; }}
  .flag-panel, .warn-panel {{ border-radius: 10px; padding: 12px 14px; margin: 12px 0; overflow-wrap: anywhere; }}
  .flag-panel ul {{ margin-bottom: 0; }}
  .flag-title {{ font-weight: 800; font-size: 1.05rem; }}
  .flag-red {{ background: #fee2e2; border: 1px solid #fecaca; }}
  .flag-yellow {{ background: #fef3c7; border: 1px solid #fde68a; }}
  .flag-green {{ background: #dcfce7; border: 1px solid #bbf7d0; }}
  .warn-panel {{ background: #fff7ed; border: 1px solid #fed7aa; }}
  .plot-scroll {{ width: 100%; max-width: 100%; overflow-x: auto; overflow-y: hidden; border: 1px solid #e5e7eb; border-radius: 10px; padding: 8px; background: #fff; margin-bottom: 12px; }}
  .plot-wide {{ width: max-content; }}
  .table-widget {{ margin: 8px 0 14px; }}
  .table-controls {{ display: flex; gap: 10px; flex-wrap: wrap; align-items: center; margin-bottom: 8px; font-size: 0.9rem; }}
  .table-controls input, .table-controls select {{ padding: 5px 7px; border: 1px solid #d1d5db; border-radius: 6px; background: #fff; }}
  .table-controls button {{ padding: 5px 9px; border: 1px solid #d1d5db; border-radius: 6px; background: #f9fafb; cursor: pointer; }}
  .table-controls button:hover {{ background: #eef2ff; }}
  .page-status {{ color: #4b5563; min-width: 120px; display: inline-block; }}
  .table-scroll {{ width: 100%; max-width: 100%; overflow-x: auto; border: 1px solid #e5e7eb; border-radius: 10px; margin: 8px 0 12px; }}
  .data-table {{ width: max-content; min-width: 100%; border-collapse: collapse; font-size: 0.88rem; }}
  .data-table th, .data-table td {{ border: 1px solid #e5e7eb; padding: 6px 8px; text-align: left; white-space: nowrap; }}
  .data-table th {{ background: #f3f4f6; position: sticky; top: 0; z-index: 1; cursor: pointer; user-select: none; }}
  .data-table th.sort-asc::after {{ content: " ▲"; font-size: 0.75rem; color: #2563eb; }}
  .data-table th.sort-desc::after {{ content: " ▼"; font-size: 0.75rem; color: #2563eb; }}
  details {{ margin: 12px 0; }}
  summary {{ cursor: pointer; font-weight: 700; margin-bottom: 8px; }}
  code {{ background: #f3f4f6; padding: 1px 4px; border-radius: 4px; overflow-wrap: anywhere; }}
  .muted {{ color: #d1d5db; }}
  .empty {{ color: #6b7280; font-style: italic; }}
  @media (max-width: 900px) {{
    .layout {{ grid-template-columns: 1fr; }}
    aside {{ position: relative; height: auto; }}
  }}
</style>
<script>
const dashboardTables = {{}};

function getDashboardTableState(tableId) {{
  if (!dashboardTables[tableId]) {{
    dashboardTables[tableId] = {{ page: 1, pageSize: 25, sortCol: null, sortDir: 1, filter: "" }};
  }}
  return dashboardTables[tableId];
}}

function parseDashboardCellValue(value) {{
  const text = (value || "").trim();
  const numeric = Number(text.replace(/,/g, ""));
  if (text !== "" && Number.isFinite(numeric)) return numeric;
  return text.toLowerCase();
}}

function dashboardTableRows(table) {{
  const tbody = table.tBodies[0];
  return tbody ? Array.from(tbody.rows) : [];
}}

function applyDashboardTable(tableId) {{
  const table = document.getElementById(tableId);
  if (!table || !table.tBodies.length) return;
  const state = getDashboardTableState(tableId);
  const tbody = table.tBodies[0];
  let rows = dashboardTableRows(table);

  if (state.sortCol !== null) {{
    rows.sort((a, b) => {{
      const av = parseDashboardCellValue(a.cells[state.sortCol]?.textContent || "");
      const bv = parseDashboardCellValue(b.cells[state.sortCol]?.textContent || "");
      if (av < bv) return -1 * state.sortDir;
      if (av > bv) return 1 * state.sortDir;
      return 0;
    }});
    rows.forEach(row => tbody.appendChild(row));
  }}

  const filtered = rows.filter(row => row.textContent.toLowerCase().includes(state.filter));
  const total = filtered.length;
  const pageSize = state.pageSize === "all" ? (total || 1) : Number(state.pageSize);
  const totalPages = Math.max(1, Math.ceil(total / pageSize));
  state.page = Math.min(Math.max(1, state.page), totalPages);

  const start = state.pageSize === "all" ? 0 : (state.page - 1) * pageSize;
  const end = state.pageSize === "all" ? total : start + pageSize;
  const visibleSet = new Set(filtered.slice(start, end));

  rows.forEach(row => {{ row.style.display = visibleSet.has(row) ? "" : "none"; }});

  const status = document.getElementById(`${{tableId}}_status`);
  if (status) status.textContent = `Page ${{state.page}} / ${{totalPages}} (${{total}} rows)`;
}}

function sortDashboardTable(tableId, colIdx) {{
  const table = document.getElementById(tableId);
  if (!table) return;
  const state = getDashboardTableState(tableId);

  if (state.sortCol === colIdx) state.sortDir *= -1;
  else {{ state.sortCol = colIdx; state.sortDir = 1; }}
  state.page = 1;

  Array.from(table.tHead?.rows[0]?.cells || []).forEach((th, idx) => {{
    th.classList.remove("sort-asc", "sort-desc");
    if (idx === state.sortCol) th.classList.add(state.sortDir === 1 ? "sort-asc" : "sort-desc");
  }});

  applyDashboardTable(tableId);
}}

function filterDashboardTable(tableId) {{
  const widget = document.querySelector(`[data-table-id="${{tableId}}"]`);
  const input = widget ? widget.querySelector(".table-search") : null;
  const state = getDashboardTableState(tableId);
  state.filter = input ? input.value.toLowerCase() : "";
  state.page = 1;
  applyDashboardTable(tableId);
}}

function setDashboardPageSize(tableId, value) {{
  const state = getDashboardTableState(tableId);
  state.pageSize = value === "all" ? "all" : Number(value);
  state.page = 1;
  applyDashboardTable(tableId);
}}

function prevDashboardPage(tableId) {{
  const state = getDashboardTableState(tableId);
  state.page -= 1;
  applyDashboardTable(tableId);
}}

function nextDashboardPage(tableId) {{
  const state = getDashboardTableState(tableId);
  state.page += 1;
  applyDashboardTable(tableId);
}}

function initDashboardTables() {{
  document.querySelectorAll("table.sortable-table").forEach(table => {{
    const tableId = table.id;
    if (!tableId) return;
    const headerRow = table.tHead?.rows[0];
    if (headerRow) {{
      Array.from(headerRow.cells).forEach((th, idx) => {{
        th.title = "Click to sort";
        th.addEventListener("click", () => sortDashboardTable(tableId, idx));
      }});
    }}
    applyDashboardTable(tableId);
  }});
}}

function filterSamples() {{
  const q = document.getElementById('sampleFilter').value.toLowerCase();
  for (const item of document.querySelectorAll('.sample-nav-item')) {{
    item.style.display = item.dataset.sample.includes(q) ? '' : 'none';
  }}
}}
function relayoutVisiblePlots() {{
  if (!window.Plotly) return;
  document.querySelectorAll('.plotly-graph-div').forEach(function(el) {{
    if (el.offsetParent !== null) {{
      try {{ Plotly.Plots.resize(el); }} catch(e) {{}}
    }}
  }});
}}
document.addEventListener('toggle', function(e) {{
  if (e.target.tagName === 'DETAILS' && e.target.open) {{
    setTimeout(relayoutVisiblePlots, 100);
  }}
}}, true);
window.addEventListener('resize', relayoutVisiblePlots);
window.addEventListener('load', function() {{ initDashboardTables(); setTimeout(relayoutVisiblePlots, 300); }});
</script>
</head>
<body>
<div class="layout">
  <aside>
    <h1>{html.escape(title)}</h1>
    <p class="muted">Root: <code>{html.escape(str(root_dir))}</code></p>
    <p class="muted">Generated: {html.escape(generated_at)}</p>
    <input id="sampleFilter" type="text" oninput="filterSamples()" placeholder="Filter samples..." />
    <ul>
      <li><a href="#overview">Overview</a></li>
      <li><a href="#charts">Summary charts</a></li>
      {''.join(sidebar_items)}
    </ul>
  </aside>
  <main>
    <section id="overview">
      <h2>Overview</h2>
      <p>This compact dashboard summarizes ge_screen outputs. Large source tables are capped in the HTML by default; use the source TSV links for complete data. Red/yellow/green flags are triage heuristics and require manual review.</p>
      <div class="cards">{cards}</div>
      {overview_table(overview_df)}
    </section>
    <section id="charts">
      <h2>Summary charts</h2>
      {charts_html}
    </section>
    {''.join(sample_sections)}
  </main>
</div>
</body>
</html>
"""
    return html_doc, overview_df


def main() -> None:
    parser = argparse.ArgumentParser(description="Build an enhanced HTML dashboard summarizing ge_screen outputs.")
    parser.add_argument("--input-dir", default=".", help="Root directory containing ge_screen outputs. Default: current directory.")
    parser.add_argument("--output-html", default=None, help="Output HTML dashboard. Default: <input-dir>/ge_screen_dashboard.html")
    parser.add_argument("--summary-tsv", default=None, help="Output summary TSV. Default: <input-dir>/ge_screen_dashboard_summary.tsv")
    parser.add_argument("--no-summary-tsv", action="store_true", help="Do not write summary TSV.")
    parser.add_argument("--title", default="Genetic Engineering Screen Dashboard", help="Dashboard title.")
    parser.add_argument("--embed-plotly", action="store_true", help="Embed Plotly JS for offline use. Adds several MB to the HTML.")
    parser.add_argument("--max-nt-rows", type=int, default=DEFAULT_MAX_NT_ROWS, help=f"Max NT rows embedded per sample. Default: {DEFAULT_MAX_NT_ROWS}. Use -1 for all.")
    parser.add_argument("--max-univec-rows", type=int, default=DEFAULT_MAX_UNIVEC_ROWS, help=f"Max UniVec rows embedded per sample. Default: {DEFAULT_MAX_UNIVEC_ROWS}. Use -1 for all.")
    parser.add_argument("--max-coverage-rows", type=int, default=DEFAULT_MAX_COVERAGE_ROWS, help=f"Max coverage rows embedded per sample. Default: {DEFAULT_MAX_COVERAGE_ROWS}. Use -1 for all.")
    parser.add_argument("--max-abricate-rows", type=int, default=DEFAULT_MAX_ABRICATE_ROWS, help=f"Max Abricate hit rows embedded per database/sample. Default: {DEFAULT_MAX_ABRICATE_ROWS}. Use -1 for all.")
    parser.add_argument("--max-amrfinder-rows", type=int, default=DEFAULT_MAX_AMRFINDER_ROWS, help=f"Max AMRFinderPlus rows embedded per sample. Default: {DEFAULT_MAX_AMRFINDER_ROWS}. Use -1 for all.")
    args = parser.parse_args()

    root = Path(args.input_dir).resolve()
    if not root.exists():
        raise SystemExit(f"Input directory does not exist: {root}")

    sample_dirs = collect_sample_dirs(root)
    if not sample_dirs:
        raise SystemExit(
            "No ge_screen output directories found. Expected files such as *_vs_nt.tsv, "
            "*_vs_univec.tsv, abricate.summary.tsv, or *.coverage.tsv under: "
            f"{root}"
        )

    def limit_value(value: int) -> int | None:
        return None if value is None or value < 0 else value

    limits = {
        "nt": limit_value(args.max_nt_rows),
        "univec": limit_value(args.max_univec_rows),
        "coverage": limit_value(args.max_coverage_rows),
        "abricate": limit_value(args.max_abricate_rows),
        "amrfinder": limit_value(args.max_amrfinder_rows),
    }

    if any(v is None for v in limits.values()):
        print("WARNING: one or more table limits are set to all rows; HTML may become very large.", file=sys.stderr)

    results = [summarize_sample(d, limits) for d in sample_dirs]
    html_doc, overview_df = make_dashboard(results, root, args.title, embed_plotly=args.embed_plotly)

    output_html = Path(args.output_html).resolve() if args.output_html else root / "ge_screen_dashboard.html"
    output_html.parent.mkdir(parents=True, exist_ok=True)
    output_html.write_text(html_doc, encoding="utf-8")
    print(f"Wrote dashboard: {output_html}")

    if not args.no_summary_tsv:
        summary_tsv = Path(args.summary_tsv).resolve() if args.summary_tsv else root / "ge_screen_dashboard_summary.tsv"
        summary_tsv.parent.mkdir(parents=True, exist_ok=True)
        overview_df.drop(columns=[c for c in ["plot_univec_size"] if c in overview_df.columns]).to_csv(summary_tsv, sep="\t", index=False)
        print(f"Wrote summary TSV: {summary_tsv}")


if __name__ == "__main__":
    main()
