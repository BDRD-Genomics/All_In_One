#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Combine per-sample ARC tables named like:
  <sample_id>_accurate_read_counts.tsv

New 9-col TSV header:
  lineage, family, below_family, Contig_Count, Total_SR_Counts,
  Total_LR_Counts, Total_Read_Count, normalized_family, normalized_rpm

Outputs (CSV):
  <project>_total_hits_merged.csv
  <project>_normalized_family_hits_merged.csv
  <project>_normalized_RPM_hits_merged.csv
"""

import pandas as pd
import os
import sys, getopt
import re

UNCLASSIFIED = "Unclassified"
FAM_RE = re.compile(r'([A-Za-z0-9][A-Za-z0-9\-_]*viridae)\b', re.IGNORECASE)

def family_from_family_col(s: str) -> str:
    """Extract the family token from the 'family' column."""
    if s is None:
        return ""
    m = FAM_RE.search(str(s))
    if not m:
        return ""
    fam = m.group(1)
    return fam[0].upper() + fam[1:]  # simple normalize

def family_from_full_lineage(full_lineage: str) -> str:
    """Fallback: scan Full_Lineage for the last …viridae token."""
    s = (full_lineage or "").replace("|", ";")
    tokens = [t.strip(" ;") for t in s.split(";") if t.strip()]
    fam = ""
    for tok in tokens:
        m = FAM_RE.search(tok)
        if m:
            fam = m.group(1)
    if fam:
        return fam[0].upper() + fam[1:]
    return ""

# CLI 
argv = sys.argv[1:]

opts, args = getopt.getopt(argv, "hr:p:o:", ["rc_path=", "project=", "outdir="])

project = "out"
rc_path = None
outdir = None

for opt,arg in opts:
    if opt == '-h':
        print("combine_ARC_hits.py -r <read_count path> -p <project name> [-o <outdir>]")
        sys.exit(0)
    elif opt in ("-r", "--rc_path"):
        rc_path = arg
    elif opt in ("-p", "--project"):
        project = arg
    elif opt in ("-o", "--outdir"):
        outdir = arg

if not rc_path:
    print("Error: --rc_path is required", file=sys.stderr)
    sys.exit(2)

rc_path = rc_path.rstrip("/")
outdir = (outdir or rc_path).rstrip("/")
os.makedirs(outdir, exist_ok=True)

pattern = "_accurate_read_counts.tsv"


arcs_tmp = os.popen("ls "+rc_path+"/*"+pattern+" 2>/dev/null").read().split('\n')
arc_files = [x for x in arcs_tmp if x]

# Base frames
main_frame = pd.DataFrame(columns=["Full_Lineage","Family", "Beyond_Family"])
main_frame_total = main_frame
main_frame_normalized = main_frame
main_frame_normalized_rpm = main_frame

empty_files=[]
all_samples = []

for file in arc_files:
    sample =  os.path.basename(file).replace("_accurate_read_counts.tsv","")
    all_samples.append(sample)

    if os.path.getsize(file) == 0:
        empty_files.append(sample)
        continue

    #  9-col header file 
    try:
        tmp_df = pd.read_csv(
            file, sep='\t', header=0, dtype=str, low_memory=False
        )
    except Exception as e:
        print(f"[WARN] Failed to read {file}: {e}", file=sys.stderr)
        empty_files.append(sample)
        continue

    # Standardize column names 
    tmp_df.columns = [c.strip() for c in tmp_df.columns]

    required = {"lineage","family","below_family","Total_Read_Count","normalized_family","normalized_rpm"}
    if not required.issubset(set(tmp_df.columns)):
        print(f"[WARN] Missing required columns in {file}: {required - set(tmp_df.columns)}", file=sys.stderr)
        empty_files.append(sample)
        continue

    tmp_df["lineage"] = tmp_df["lineage"].fillna("").astype(str).str.strip().str.rstrip(";")
    tmp_df["family"]  = tmp_df["family"].replace({"N/A": "","No family found": ""}).fillna("").astype(str).str.strip()
    tmp_df["below_family"] = tmp_df["below_family"].replace({"N/A": ""}).fillna("").astype(str).str.strip().str.rstrip(";")

    # taxonomy
    tmp_df["Full_Lineage"]  = tmp_df["lineage"].where(tmp_df["below_family"].eq(""),
                                                      tmp_df["lineage"] + ";" + tmp_df["below_family"])
    tmp_df["Beyond_Family"] = tmp_df["below_family"]

    #Use family column
    fam1 = tmp_df["family"].apply(family_from_family_col)
    fam2 = tmp_df["Full_Lineage"].apply(family_from_full_lineage)
    tmp_df["Family"] = fam1.where(fam1.ne(""), fam2)
    tmp_df.loc[tmp_df["Family"].eq(""), "Family"] = UNCLASSIFIED
    tmp_df = tmp_df[tmp_df["Full_Lineage"] != "TOTAL READ COUNT:"].copy()

    # convert to numeric metrics
    for c in ("Total_Read_Count", "normalized_family", "normalized_rpm"):
        tmp_df[c] = pd.to_numeric(tmp_df[c], errors="coerce")

    num_rows = len(tmp_df)
    if num_rows == 0:
        empty_files.append(sample)
        continue

    # Add sample-specific metric columns
    tmp_df[sample+"_total"] = tmp_df["Total_Read_Count"]
    tmp_df[sample+"_normalizedFamily"] = tmp_df["normalized_family"]
    tmp_df[sample+"_normalizedRPM"] = tmp_df["normalized_rpm"]

    # Subsets for merge
    subdf_normalizedFam = tmp_df[["Full_Lineage","Family","Beyond_Family", sample+"_normalizedFamily"]]
    subdf_normalizedRPM = tmp_df[["Full_Lineage","Family","Beyond_Family", sample+"_normalizedRPM"]]
    subdf_total         = tmp_df[["Full_Lineage","Family","Beyond_Family", sample+"_total"]]

    # Merge (outer) 
    main_frame_normalized     = pd.merge(main_frame_normalized, subdf_normalizedFam, how="outer",
                                         on=["Full_Lineage","Family","Beyond_Family"])
    main_frame_normalized_rpm = pd.merge(main_frame_normalized_rpm, subdf_normalizedRPM, how="outer",
                                         on=["Full_Lineage","Family","Beyond_Family"])
    main_frame_total          = pd.merge(main_frame_total, subdf_total, how="outer",
                                         on=["Full_Lineage","Family","Beyond_Family"])


add_cols = list(set([x for x in empty_files if empty_files.count(x)==3]))
for sam in add_cols:
    main_frame_total[sam] = None
    main_frame_normalized[sam] = None

print("Total number of samples:", len(all_samples))

# Write outputs 
main_frame_total.to_csv(os.path.join(outdir, f"{project}_total_hits_merged.csv"), index=False)
main_frame_normalized.to_csv(os.path.join(outdir, f"{project}_normalized_family_hits_merged.csv"), index=False)
main_frame_normalized_rpm.to_csv(os.path.join(outdir, f"{project}_normalized_RPM_hits_merged.csv"), index=False)