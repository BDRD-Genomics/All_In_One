#!/usr/bin/env python
"""
extract_review.py — pull the calls an analyst should review from an annotated
CSV: the model-only VF calls with NO homology support (review_priority
'3-review-no-homology'). These are the candidate-novel-VFs OR false positives —
the ones that actually need a human look, separated from the confirmed/known.

Writes a focused TSV/CSV sorted by model score (highest-confidence review items
first), plus a short console tally.

By default it targets '3-review-no-homology'. Use --priority to grab a different
review_priority bucket, or --all-review to grab everything that isn't a
confirmed known VF.

USAGE
  python extract_review.py --annotated sample_annotated.csv --out sample_review.tsv
  python extract_review.py --annotated sample_annotated.csv --out sample_review.tsv --all-review
"""
import argparse, csv, sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--annotated", required=True, help="annotated CSV from interpret_calls.py")
    ap.add_argument("--out", required=True, help="output file for the flagged rows")
    ap.add_argument("--priority", default="3-review-no-homology",
                    help="review_priority value to extract (default: 3-review-no-homology)")
    ap.add_argument("--all-review", action="store_true",
                    help="grab ALL non-confirmed calls needing review "
                         "(everything except priority starting with '1-' or '2-')")
    args = ap.parse_args()

    # annotated CSV has leading '#' comment lines — skip them
    with open(args.annotated) as f:
        lines = [ln for ln in f if not ln.startswith("#")]
    rdr = csv.DictReader(lines)
    fieldnames = rdr.fieldnames
    if "review_priority" not in fieldnames:
        sys.exit("ERROR: no 'review_priority' column — is this an annotated CSV "
                 "from interpret_calls.py?")

    rows = list(rdr)
    if args.all_review:
        keep = [r for r in rows
                if r["call"].upper() == "VF"
                and not r["review_priority"].startswith(("1-", "2-"))]
        label = "all non-confirmed review items"
    else:
        keep = [r for r in rows if r["review_priority"] == args.priority]
        label = f"review_priority == '{args.priority}'"

    # sort by model score descending — highest-confidence flags first
    def score(r):
        try: return float(r["model_score"])
        except (KeyError, ValueError): return 0.0
    keep.sort(key=score, reverse=True)

    # write (tab-separated for easy reading; keeps the useful columns up front)
    want_cols = ["protein_id", "model_score", "confidence_band", "model_category",
                 "homology_status", "review_priority", "interpretation"]
    out_cols = [c for c in want_cols if c in fieldnames] + \
               [c for c in fieldnames if c not in want_cols]
    delim = "\t" if args.out.endswith((".tsv", ".txt")) else ","
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=out_cols, delimiter=delim, extrasaction="ignore")
        w.writeheader()
        w.writerows(keep)

    # console tally
    total_vf = sum(1 for r in rows if r["call"].upper() == "VF")
    from collections import Counter
    cat_tally = Counter(r["model_category"] for r in keep if r.get("model_category"))
    print(f"[extract_review] {args.annotated}")
    print(f"  filter: {label}")
    print(f"  flagged for review: {len(keep)}  (of {total_vf} total VF calls)")
    if cat_tally:
        print(f"  by predicted category:")
        for cat, c in cat_tally.most_common():
            print(f"    {cat:<32} {c}")
    print(f"  wrote {args.out}  (sorted by model score, highest first)")
    print(f"  NOTE: these are model calls with NO VFDB homology — candidate novel")
    print(f"        VFs OR false positives. They REQUIRE analyst validation; they")
    print(f"        are not confirmed virulence factors.")


if __name__ == "__main__":
    main()
