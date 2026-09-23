#!/usr/bin/env python3
##################################################################################################
#                                                                                                #
#                                The purpose of this script is to                                #
#                   take the blast contigs containing at least 1 viral result and                #
#                           sort the output into ambiguous, phage and viral.                     #
#                                                                                                #
##################################################################################################

import sys
import os
import argparse
import time
import subprocess as sp

import pandas as pd
from ete3 import NCBITaxa


def parse_args():
    parser = argparse.ArgumentParser(
        prog='Filter Blast Viruses',
        usage=('filter_blast_viruses_no_biopython.py -t <viral_blast_tab> -c <contigs_fasta> '
               '-p <PE_covstats> -s <SE_covstats> -l <LR_covstats> '
               '-S <PE_QC_reads> -L <LR_QC_reads> '
               '-n <ncbidb> -v <virusFamily_genomeSizes.txt> -b <blacklist>'),
        description=('Filter putatively viral blast outputs to high confidence viral calls. '
                     'Outputs viral blast table, phage blast table and ambiguous blast tables. '
                     'Also outputs viral contigs and semiquantitative read counts.'),
        epilog='Parse the blast table from sequence comparative algorithm (MMSeqs, blastn, blastx)')
    # options
    parser.add_argument("--sample-id", required=True,
                        help="Sample identifier; will prefix all output files")
    parser.add_argument("-s", "--se_cov",
                        help="Single-end coverage stats (optional)")
    parser.add_argument("-p", "--pe_cov",
                        help="Paired-end coverage stats (optional)")
    parser.add_argument("-l", "--lr_cov",
                        help="Long-read coverage stats (optional)")
    parser.add_argument("-t", "--viral_blasttab",
                        help="Filtered viral BLAST/MMseqs table", required=True)
    parser.add_argument("-c", "--contigs",
                        help="Contigs FASTA file", required=True)
    parser.add_argument("-o", "--outdir",
                        help="Output directory (default: .)")
    parser.add_argument("-L", "--qc_count_LR",
                        help="File with a single number: total QC long reads")
    parser.add_argument("-S", "--qc_count_SR",
                        help="File with a single number: total QC short reads")
    parser.add_argument("-n", "--ncbidb",
                        help="Path to ete3 NCBITaxa sqlite DB", required=True)
    parser.add_argument("-v", "--viral_genome_sizes",
                        help="virusFamily_genomeSizes.txt", required=True)
    parser.add_argument("-b", "--blacklist",
                        help="Optional blacklist of problematic accessions (one per line)")
    return parser.parse_args()


def read_fasta_dict(fasta_path):
    """Read a FASTA into a dict {seq_id: sequence} without Biopython."""
    seqs = {}
    current_id = None
    chunks = []

    with open(fasta_path, 'r') as handle:
        for line in handle:
            line = line.rstrip('\n')
            if not line:
                continue
            if line.startswith('>'):
                # flush previous
                if current_id is not None:
                    seqs[current_id] = ''.join(chunks)
                # take the first token as the ID (up to first whitespace)
                current_id = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)

    # flush last record
    if current_id is not None:
        seqs[current_id] = ''.join(chunks)

    return seqs


def get_family(ncbi, full_taxid):
    parsed = full_taxid.split(';')
    # remove trailing empty field if present
    if parsed and parsed[-1] == '':
        parsed.pop()

    if not parsed:
        return "N/A", "N/A"

    try:
        taxid_dict = ncbi.get_name_translator([parsed[-1]])
        lineage = ncbi.get_lineage(taxid_dict.get(parsed[-1])[0])
        rank = ncbi.get_rank(lineage)
    except Exception:
        return "N/A", "N/A"

    family_key = list(rank.keys())
    family_value = list(rank.values())
    try:
        family_idx = family_value.index("family")
        family_dict = ncbi.get_taxid_translator([family_key[family_idx]])
        family = family_dict.get(family_key[family_idx])
        if family is None:
            return "N/A", "N/A"
        try:
            taxid_family_index = parsed.index(family)
            sub_family = ';'.join(parsed[taxid_family_index + 1:])
            return family, sub_family
        except ValueError:
            return family, "N/A"
    except ValueError:
        # 'family' not in ranks
        return "N/A", "N/A"
    except Exception:
        return "N/A", "N/A"


def normalization(counts, qc_reads, genome_length):
    if qc_reads == 0 or genome_length in (0, None):
        return "N/A"
    return ((counts / qc_reads) / genome_length) * 1000.0


def rpm_norm(counts, qc_reads):
    if qc_reads == 0:
        return 0.0
    return (counts / qc_reads) * 1_000_000.0


def main():
    args = parse_args()
    print(time.ctime())
    print(args.se_cov, args.pe_cov, args.lr_cov, args.viral_blasttab,
          args.contigs, args.outdir, args.qc_count_SR, args.qc_count_LR)

    filtered_blast_input = args.viral_blasttab
    contigs_file = args.contigs
    outdir = args.outdir if args.outdir else "."
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    SE_cov_file = args.se_cov if args.se_cov else None
    PE_cov_file = args.pe_cov if args.pe_cov else None
    LR_cov_file = args.lr_cov if args.lr_cov else None

    blacklist_file = args.blacklist if args.blacklist else None
    bad_acc_list = []
    if blacklist_file:
        print("Filtering by blacklist")
        with open(blacklist_file, 'r') as bad_acc_file:
            bad_acc_list = bad_acc_file.read().splitlines()

    # PE QC counts
    if args.qc_count_SR is None:
        print("PE_QC_counts")
        PE_QC_counts = 0
        print(PE_QC_counts)
    else:
        PE_QC_counts = pd.read_csv(args.qc_count_SR, header=None)
        PE_QC_counts = PE_QC_counts.iloc[0, 0]
        print("PE_QC_counts")
        print(PE_QC_counts)

    # LR QC counts
    if args.qc_count_LR is None:
        print("LR_QC_counts")
        LR_QC_counts = 0
        print(LR_QC_counts)
    else:
        LR_QC_counts = pd.read_csv(args.qc_count_LR, header=None)
        LR_QC_counts = LR_QC_counts.iloc[0, 0]
        print("LR_QC_counts")
        print(LR_QC_counts)

    print("PE_QC_counts")
    print(PE_QC_counts)
    print("LR_QC_counts")
    print(LR_QC_counts)
    QC_counts = PE_QC_counts + LR_QC_counts
    print("QC Counts:" + str(QC_counts))

    # coverage tables
    if PE_cov_file:
        PE_cov_df = pd.read_csv(
            PE_cov_file, delimiter="\t", low_memory=True, engine="pyarrow",
            names=["query_ID", "qlen", "SR_counts", "SR_unpaired_counts"]
        )
    else:
        PE_cov_df = None

    if LR_cov_file:
        LR_cov_df = pd.read_csv(
            LR_cov_file, delimiter="\t", low_memory=True, engine="pyarrow",
            names=["query_ID", "qlen", "LR_counts", "LR_unpaired_counts"]
        )
    else:
        LR_cov_df = None

    ###########################################################################
    # Step 1: filter putative viral results into: Viral, Phage, and Ambiguous #
    ###########################################################################
    blast_df = pd.read_csv(filtered_blast_input, delimiter="\t",
                           low_memory=True, engine="pyarrow")

    # blacklist flag column (always present so downstream code is simpler)
    if blacklist_file and bad_acc_list:
        blast_df["poor_accessions"] = blast_df["ref_description"].apply(
            lambda x: any(acc in str(x) for acc in bad_acc_list)
        )
    else:
        blast_df["poor_accessions"] = False

    phage_keywords = [
        "phage", "Caudoviricetes", "Microviridae",
        "unclassified bacterial viruses", "Leviviricetes", "Tubulavirales"
    ]

    sub_blast_df = blast_df.groupby("query_ID").agg(
        lineages=("lineage", lambda x: list(x)),
        descriptions=("ref_description", lambda x: list(x)),
        poor_accession_only=("poor_accessions", lambda x: all(list(x)))
    ).reset_index()

    sub_blast_df["cellular_organisms"] = sub_blast_df["lineages"].apply(
        lambda x: any("cellular organisms" in str(s) for s in x)
    )
    sub_blast_df["undefined_taxon"] = sub_blast_df["lineages"].apply(
        lambda x: any("undefined taxon" in str(s) for s in x)
    )
    sub_blast_df["phages"] = sub_blast_df["lineages"].apply(
        lambda x: any(k in str(s) for s in x for k in phage_keywords)
    )
    sub_blast_df["viral_only"] = sub_blast_df.eval(
        "cellular_organisms == False and phages == False and poor_accession_only == False"
    )

    viral_contig_ids = sub_blast_df.loc[
        (sub_blast_df["viral_only"] == True) &
        (sub_blast_df["undefined_taxon"] == False),
        "query_ID"
    ].tolist()
    undefined_contig_ids = sub_blast_df.loc[
        (sub_blast_df["viral_only"] == True) &
        (sub_blast_df["undefined_taxon"] == True),
        "query_ID"
    ].tolist()
    phage_contig_ids = sub_blast_df.loc[
        (sub_blast_df["phages"] == True) &
        (sub_blast_df["cellular_organisms"] == False),
        "query_ID"
    ].tolist()
    ambiguous_contig_ids = sub_blast_df.loc[
        sub_blast_df["cellular_organisms"] == True,
        "query_ID"
    ].tolist()
    poor_acc_contig_ids = sub_blast_df.loc[
        sub_blast_df["poor_accession_only"] == True,
        "query_ID"
    ].tolist()

    print("Undefined contig ID's: " + str(len(undefined_contig_ids)))

    # Reassign some undefined taxon contigs based on presence of 'virus' in description
    blast_undefined_temp = blast_df[
        blast_df["query_ID"].isin(undefined_contig_ids) &
        blast_df["lineage"].str.contains("undefined taxon")
    ]
    blast_undefined_temp["contains_virus"] = blast_undefined_temp["ref_description"].str.contains(
        "virus", case=False, na=False
    )
    blast_undefined_temp_sub = blast_undefined_temp.groupby("query_ID").agg(
        contains_virus=("contains_virus", lambda x: list(x))
    ).reset_index()
    blast_undefined_temp_sub["check_virus"] = blast_undefined_temp_sub["contains_virus"].apply(
        lambda x: all(x)
    )

    print("PRE: Viral contig ID's: " + str(len(viral_contig_ids)))
    print("PRE: Phage contig ID's: " + str(len(phage_contig_ids)))
    print("PRE: Ambiguous contig ID's: " + str(len(ambiguous_contig_ids)))
    print("PRE: Undefinied contig ID's: " + str(len(undefined_contig_ids)))

    viral_contig_ids += blast_undefined_temp_sub.loc[
        blast_undefined_temp_sub["check_virus"] == True, "query_ID"
    ].tolist()
    ambiguous_contig_ids += blast_undefined_temp_sub.loc[
        blast_undefined_temp_sub["check_virus"] == False, "query_ID"
    ].tolist()

    blast_viral = blast_df[blast_df["query_ID"].isin(viral_contig_ids)]

    if blacklist_file and bad_acc_list:
        bad_acc_pattern = "|".join(bad_acc_list)
        blast_viral = blast_viral[
            ~blast_viral['ref_description'].str.contains(
                bad_acc_pattern, case=False, na=False
            )
        ]

    viral_mask = ~blast_viral.duplicated(subset=["query_ID"], keep='first')
    blast_viral_top_results = blast_viral[viral_mask].sort_values(by="lineage")
    blast_viral_top_results.to_csv(
        os.path.join(outdir, f"{args.sample_id}_viral_blast_contigs.csv"),
        sep=",", index=False
    )

    blast_phage = blast_df[blast_df["query_ID"].isin(phage_contig_ids)]
    phage_mask = ~blast_phage.duplicated(subset=["query_ID"], keep='first')
    blast_phage_top_results = blast_phage[phage_mask].sort_values(by="lineage")
    blast_phage_top_results.to_csv(
        os.path.join(outdir, f"{args.sample_id}_phage_blast_contigs.csv"),
        sep=",", index=False
    )

    blast_ambiguous = blast_df[blast_df["query_ID"].isin(ambiguous_contig_ids)]
    blast_ambiguous.to_csv(
        os.path.join(outdir, f"{args.sample_id}_full_ambiguous_viral_blast_contigs.csv"),
        sep=",", index=False
    )

    blast_bad_acc = blast_df[blast_df["query_ID"].isin(poor_acc_contig_ids)]
    blast_bad_acc.to_csv(
        os.path.join(outdir, f"{args.sample_id}_full_bad_accessions_blast_contigs.csv"),
        sep=",", index=False
    )

    print("Viral contig ID's: " + str(len(viral_contig_ids)))
    print("Phage contig ID's: " + str(len(phage_contig_ids)))
    print("Ambiguous contig ID's: " + str(len(ambiguous_contig_ids)))
    print("Bad accessions only contig ID's: " + str(len(poor_acc_contig_ids)))

    ###########################################################################
    # Step 2: Extract contigs for high confidence viral results               #
    ###########################################################################

    print("Reading contigs file: " + contigs_file)
    seq_dict = read_fasta_dict(contigs_file)
    viral_fasta_file = "viral_contigs.fasta"

    try:
        out_fa_path = os.path.join(outdir, viral_fasta_file)
        outFile = open(out_fa_path, "w")
    except OSError:
        print("Can't create output file : " + viral_fasta_file)
        sys.exit(1)

    print("Writing contigs to output file: " + out_fa_path)
    missing = 0
    for seq_name in viral_contig_ids:
        seq = seq_dict.get(seq_name)
        if seq is None:
            missing += 1
            continue
        outFile.write(f">{seq_name}\n{seq}\n")
    outFile.close()
    if missing:
        print(f"Warning: {missing} viral contig IDs had no sequence in {contigs_file}")

    ###########################################################################
    # Step 3: Calculate semi-quantitative abundance for read counts report    #
    ###########################################################################

    # If we have no coverage files at all, stop here
    if not any([SE_cov_file, PE_cov_file, LR_cov_file]):
        print("No coverage stats provided; skipping read count calculation.")
        print("VS parsing finished.")
        print(time.ctime())
        return

    vf_dict = pd.read_csv(
        args.viral_genome_sizes, sep="\t", header=0, names=["fam", "size"]
    ).set_index('fam')['size'].to_dict()

    ncbi = NCBITaxa(dbfile=args.ncbidb)

    print("Length of pandas vf dataframe: " + str(len(blast_viral_top_results)))

    def calculate_read_counts():
        if PE_cov_df is not None:
            merged = pd.merge(
                blast_viral_top_results,
                PE_cov_df[["query_ID", "SR_counts"]],
                on="query_ID",
                how="left"
            )
        else:
            merged = blast_viral_top_results.copy()
            merged["SR_counts"] = 0

        if LR_cov_df is not None:
            merged2 = pd.merge(
                merged,
                LR_cov_df[["query_ID", "LR_counts"]],
                on="query_ID",
                how="left"
            )
        else:
            merged2 = merged.copy()
            merged2["LR_counts"] = 0

        merged2["SR_counts"] = merged2["SR_counts"].fillna(0).astype(int)
        merged2["LR_counts"] = merged2["LR_counts"].fillna(0).astype(int)

        merged2.loc[
            merged2["lineage"] == "undefined taxon", "lineage"
        ] = merged2["ref_description"]

        count_df = merged2.groupby("lineage").agg(
            Contig_Count=("lineage", "count"),
            Total_SR_Counts=("SR_counts", "sum"),
            Total_LR_Counts=("LR_counts", "sum")
        ).reset_index()
        count_df["Total_Read_Count"] = (
            count_df["Total_SR_Counts"] + count_df["Total_LR_Counts"]
        )
        count_df = count_df[count_df["Total_Read_Count"] > 0]
        return count_df

    if len(blast_viral_top_results) == 0:
        print("No viral results, skipping read count calculation")
        print("VS parsing finished.")
        print(time.ctime())
        return

    print("Calculating Read Counts....")
    count_df = calculate_read_counts()

    def qkt_function(row):
        if str(row["lineage"]).startswith("Viruses;"):
            family, below_family = get_family(ncbi, row["lineage"])
        else:
            family, below_family = "N/A", "N/A"

        if family != "N/A":
            mean_genome_family = vf_dict.get(family)
            if mean_genome_family:
                norm_family = normalization(
                    row["Total_Read_Count"], QC_counts, mean_genome_family
                )
            else:
                norm_family = "N/A"
        else:
            norm_family = "N/A"

        rpm = rpm_norm(row["Total_Read_Count"], QC_counts)
        return family, below_family, norm_family, rpm

    print("Adding metadata to counts table...")
    if len(blast_viral_top_results) != 0:
        count_df[["family", "below_family",
                  "normalized_family", "normalized_rpm"]] = count_df.apply(
            qkt_function, axis=1, result_type="expand"
        )
        col_order = [
            "lineage", "family", "below_family", "Contig_Count",
            "Total_SR_Counts", "Total_LR_Counts", "Total_Read_Count",
            "normalized_family", "normalized_rpm"
        ]
        out_counts = os.path.join(
            outdir, f"{args.sample_id}_accurate_read_counts.tsv"
        )
        count_df[col_order].to_csv(out_counts, sep="\t", index=False)
    else:
        with open(os.path.join(outdir, "no_viruses.txt"), 'w') as file:
            file.write("No viruses were found.\n")

    print("VS parsing finished.")
    print(time.ctime())


if __name__ == "__main__":
    main()
