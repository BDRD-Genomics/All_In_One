##################################################################################################
#                                                                                                #
#                                The purpose of this script is to                                #
#                   take the blast contigs containing at least 1 viral result and                #
#                           sort the output into ambiguous, phage and viral.                     #
#                                                                                                #
##################################################################################################
import sys
import os
import re
import pandas as pd
from Bio import SeqIO
import argparse
from ete3 import NCBITaxa
import time

parser = argparse.ArgumentParser(
    prog='Filter Blast Viruses',
    usage='filter_blast_viruses.py -t <viral_blast_tab> -c <contigs_fasta> -p <PE_covstats> -S <PE_QC_reads> -l <LR_covstats> -L <LR_QC_reads> -n <ncbidb> -v <virusFamily_genomeSizes.txt>',
    description='Filter putatively viral blast outputs to high confidence viral calls. Outputs viral blast table, phage blast table and ambiguous blast tables. Also outputs viral contigs and semiquantitative read counts.',
    epilog='Parse the blast table from sequence comparative algorithm (MMSeqs, blastn, blastx)'
)

parser.add_argument("--sample-id", required=True,
                    help="Sample identifier; will prefix all output files")
parser.add_argument("-s", "--se_cov")
parser.add_argument("-p", "--pe_cov")
parser.add_argument("-l", "--lr_cov")
parser.add_argument("-t", "--viral_blasttab")
parser.add_argument("-c", "--contigs")
parser.add_argument("-o", "--outdir")
parser.add_argument("-L", "--qc_count_LR")
parser.add_argument("-S", "--qc_count_SR")
parser.add_argument("-n", "--ncbidb")
parser.add_argument("-v", "--viral_genome_sizes")
parser.add_argument("-b", "--blacklist")

args = parser.parse_args()
print(time.ctime())
print(args.se_cov, args.pe_cov, args.lr_cov, args.viral_blasttab, args.contigs, args.outdir, args.qc_count_SR, args.qc_count_LR)

if args.viral_blasttab is None or args.contigs is None:
    print(parser.usage)
    sys.exit(0)

virus_genomes_file = args.viral_genome_sizes
ncbi_file = args.ncbidb
ncbi = NCBITaxa(dbfile=ncbi_file)

filtered_blast_input = args.viral_blasttab
contigs_file = args.contigs
outdir = args.outdir if args.outdir else "."
os.makedirs(outdir, exist_ok=True)

SE_cov_file = args.se_cov if args.se_cov else None
PE_cov_file = args.pe_cov if args.pe_cov else None
LR_cov_file = args.lr_cov if args.lr_cov else None  
blacklist_file = args.blacklist if args.blacklist else None

bad_acc_list = []
if blacklist_file:
    print("filtering by blacklist")
    with open(blacklist_file, 'r') as bad_acc_file:
        bad_acc_list = [line.strip() for line in bad_acc_file if line.strip()]


def read_count_file(path):
    if not path:
        return 0
    try:
        df = pd.read_csv(path, header=None)
        if df.empty:
            return 0
        return int(df.iloc[0, 0])
    except Exception as exc:
        print(f"WARNING: could not read QC count file {path}: {exc}")
        return 0


PE_QC_counts = read_count_file(args.qc_count_SR)
LR_QC_counts = read_count_file(args.qc_count_LR)

print("PE_QC_counts")
print(PE_QC_counts)
print("LR_QC_counts")
print(LR_QC_counts)
QC_counts = PE_QC_counts + LR_QC_counts
print("QC Counts:" + str(QC_counts))

if PE_cov_file:
    PE_cov_df = pd.read_csv(PE_cov_file, delimiter="\t", names=["qseqid", "qlen", "SR_counts", "SR_unpaired_counts"])
    PE_cov_df["qseqid"] = PE_cov_df["qseqid"].fillna("").astype(str)
else:
    PE_cov_df = None

if LR_cov_file:
    LR_cov_df = pd.read_csv(LR_cov_file, delimiter="\t", names=["qseqid", "qlen", "LR_counts", "LR_unpaired_counts"])
    LR_cov_df["qseqid"] = LR_cov_df["qseqid"].fillna("").astype(str)
else:
    LR_cov_df = None

###########################################################################
# Step 1: filter putative viral results into: Viral, Phage, and Ambiguous #
###########################################################################

blast_df = pd.read_csv(filtered_blast_input, delimiter=",", low_memory=False)

required_cols = ["qseqid", "sseqid", "sallgi", "lineage"]
missing_cols = [col for col in required_cols if col not in blast_df.columns]
if missing_cols:
    raise ValueError(f"Input blast table is missing required columns: {missing_cols}")

# FIXED: protect string operations from NaN/float values.
blast_df["lineage"] = blast_df["lineage"].fillna("undefined taxon").astype(str)
for col in required_cols:
    blast_df[col] = blast_df[col].fillna("").astype(str)

phage_keywords = [
    "phage",
    "Caudoviricetes",
    "Microviridae",
    "unclassified bacterial viruses",
    "Leviviricetes",
    "Tubulavirales",
]

if bad_acc_list:
    blast_df["poor_accessions"] = blast_df["sseqid"].apply(
        lambda x: any(acc in x for acc in bad_acc_list)
    )
else:
    blast_df["poor_accessions"] = False

sub_blast_df = (
    blast_df.groupby("qseqid")
    .agg(
        lineages=("lineage", lambda x: list(x)),
        descriptions=("sallgi", lambda x: list(x)),
        poor_accession_only=("poor_accessions", lambda x: all(list(x))),
    )
    .reset_index()
)

def assign_sequence_source(seqid):
    if seqid.startswith(("M07012","M09150","M09151","SH00683","VH02191","A01716","FS10001335")): # miseq, miseq, miseq, miseqi100, nextseq, novaseq, iseq
        return "short_read"
    elif seqid.startswith(("NODE","contig")): # | str(seqid).isdigit(): i'm going to have to come back to this
        return "contig"
    elif seqid.count("-") == 4:
        return "long_read"
    
print(sub_blast_df)
sub_blast_df["lineages"] = sub_blast_df["lineages"].apply(lambda x: ["" if pd.isna(i) else str(i) for i in x])
sub_blast_df["cellular_organisms"] = sub_blast_df["lineages"].apply(lambda x: any("cellular organisms" in s for s in x))
sub_blast_df["undefined_taxon"] = sub_blast_df["lineages"].apply(lambda x: any("undefined taxon" in s for s in x))
sub_blast_df["phages"] = sub_blast_df["lineages"].apply(lambda x: any(k in s for s in x for k in phage_keywords))
sub_blast_df["viral_only"] = sub_blast_df.eval(
    "cellular_organisms == False and phages == False and poor_accession_only == False"
)

viral_contig_ids = sub_blast_df.loc[(sub_blast_df["viral_only"] == True) & (sub_blast_df["undefined_taxon"] == False), "qseqid"].tolist()
undefined_contig_ids = sub_blast_df.loc[(sub_blast_df["viral_only"] == True) & (sub_blast_df["undefined_taxon"] == True), "qseqid"].tolist()
phage_contig_ids = sub_blast_df.loc[(sub_blast_df["phages"] == True) & (sub_blast_df["cellular_organisms"] == False), "qseqid"].tolist()
ambiguous_contig_ids = sub_blast_df.loc[sub_blast_df["cellular_organisms"] == True, "qseqid"].tolist()
poor_acc_contig_ids = sub_blast_df.loc[sub_blast_df["poor_accession_only"] == True, "qseqid"].tolist()

print("Undefined contig ID's: " + str(len(undefined_contig_ids)))

blast_undefined_temp = blast_df[
    blast_df["qseqid"].isin(undefined_contig_ids)
    & blast_df["lineage"].str.contains("undefined taxon", na=False)
].copy()

if not blast_undefined_temp.empty:
    blast_undefined_temp["contains_virus"] = blast_undefined_temp["sallgi"].str.contains("virus", case=False, na=False)
    blast_undefined_temp_sub = (
        blast_undefined_temp.groupby("qseqid")
        .agg(contains_virus=("contains_virus", lambda x: list(x)))
        .reset_index()
    )
    blast_undefined_temp_sub["check_virus"] = blast_undefined_temp_sub["contains_virus"].apply(lambda x: all(x))
else:
    blast_undefined_temp_sub = pd.DataFrame(columns=["qseqid", "contains_virus", "check_virus"])

print("PRE: Viral contig ID's: " + str(len(viral_contig_ids)))
print("PRE: Phage contig ID's: " + str(len(phage_contig_ids)))
print("PRE: Ambiguous contig ID's: " + str(len(ambiguous_contig_ids)))
print("PRE: Undefinied contig ID's: " + str(len(undefined_contig_ids)))

viral_contig_ids += blast_undefined_temp_sub.loc[blast_undefined_temp_sub["check_virus"] == True, "qseqid"].tolist()
ambiguous_contig_ids += blast_undefined_temp_sub.loc[blast_undefined_temp_sub["check_virus"] == False, "qseqid"].tolist()

blast_viral = blast_df[blast_df["qseqid"].isin(viral_contig_ids)].copy()
# fix 
if bad_acc_list:
    bad_acc_pattern = "|".join(re.escape(acc) for acc in bad_acc_list)
    blast_viral = blast_viral[
        ~blast_viral["sallgi"].str.contains(bad_acc_pattern, case=False, na=False, regex=True)
    ]

# remove excessive undefined taxon
if not blast_viral.empty:
    dropif = blast_viral["lineage"] == "undefined taxon"
    all_drop_if = blast_viral.groupby("qseqid")["lineage"].transform(lambda x: (x == "undefined taxon").all())
    rows_to_keep_mask = ~dropif | (dropif & all_drop_if)
    blast_viral_test = blast_viral[rows_to_keep_mask]
else:
    blast_viral_test = blast_viral

viral_mask = ~blast_viral_test.duplicated(subset=["qseqid"], keep="first")
blast_viral_top_results = blast_viral_test[viral_mask].sort_values(by="lineage")
blast_viral_top_results.to_csv(f"{outdir}/{args.sample_id}_viral_blast_contigs.csv", sep=",", index=False)
blast_viral_top_results["seq_source"] = blast_viral_top_results["qseqid"].apply(lambda x: assign_sequence_source(x))


blast_phage = blast_df[blast_df["qseqid"].isin(phage_contig_ids)].copy()
phage_mask = ~blast_phage.duplicated(subset=["qseqid"], keep="first")
blast_phage_top_results = blast_phage[phage_mask].sort_values(by="lineage")
blast_phage_top_results.to_csv(f"{outdir}/{args.sample_id}_phage_blast_contigs.csv", sep=",", index=False)

blast_ambiguous = blast_df[blast_df["qseqid"].isin(ambiguous_contig_ids)].copy()
blast_ambiguous.to_csv(f"{outdir}/{args.sample_id}_full_ambiguous_viral_blast_contigs.csv", sep=",", index=False)

blast_bad_accessions = blast_df[blast_df["qseqid"].isin(poor_acc_contig_ids)].copy()
blast_bad_accessions.to_csv(f"{outdir}/{args.sample_id}_full_bad_accessions_blast_contigs.csv", sep=",", index=False)

print("Viral contig ID's: " + str(len(viral_contig_ids)))
print("Phage contig ID's: " + str(len(phage_contig_ids)))
print("Ambiguous contig ID's: " + str(len(ambiguous_contig_ids)))
print("Bad accessions only contig ID's: " + str(len(poor_acc_contig_ids)))

###########################################################################
# Step 2: Extract contigs for high confidence viral results               #
###########################################################################


def read_FASTA_data(fastaFile):
    return SeqIO.index(fastaFile, "fasta")


print("Reading contigs file: " + contigs_file)
seq_dict = read_FASTA_data(contigs_file)
viral_fasta_file = "viral_contigs.fasta"
viral_fasta_path = os.path.join(outdir, viral_fasta_file)

print("Writing contigs to output file: " + viral_fasta_file)
missing_fasta_ids = []
with open(viral_fasta_path, "w") as outFile:
    for seq_name in viral_contig_ids:
        if seq_name in seq_dict:
            outFile.write(">" + seq_name + "\n")
            outFile.write(str(seq_dict[seq_name].seq) + "\n")
        else:
            missing_fasta_ids.append(seq_name)

if missing_fasta_ids:
    print(f"WARNING: {len(missing_fasta_ids)} viral IDs were not found in the FASTA and were skipped")

###########################################################################
# Step 3: Calculate semi-quantitative abundance for read counts report    #
###########################################################################

if not any([SE_cov_file, PE_cov_file, LR_cov_file]):
    print("No covstats files were provided; skipping read count calculation")
    print("VS parsing finished.")
    print(time.ctime())
    sys.exit(0)

vf_dict = pd.read_csv(virus_genomes_file, sep="\t", header=0, names=["fam", "size"]).set_index("fam")["size"].to_dict()


def get_family(full_taxid):
    parsed = str(full_taxid).split(";")
    if parsed and parsed[-1] == "":
        parsed.pop()
    if not parsed:
        return "N/A", "N/A"

    try:
        taxid = ncbi.get_name_translator([parsed[-1]])
        lineage = ncbi.get_lineage(taxid.get(parsed[-1])[0])
        rank = ncbi.get_rank(lineage)

        family_key = list(rank.keys())
        family_value = list(rank.values())
        family_idx = family_value.index("family")
        family_dict = ncbi.get_taxid_translator([family_key[family_idx]])
        family = family_dict.get(family_key[family_idx])
        try:
            taxid_family_index = parsed.index(family)
            sub_family = ";".join(parsed[taxid_family_index + 1:])
            return family, sub_family
        except ValueError:
            return family, "N/A"
    except Exception:
        return "N/A", "N/A"


def normalization(counts, QC_reads, genome_length):
    if QC_reads == 0 or genome_length in [0, None, ""]:
        return 0
    return ((counts / QC_reads) / genome_length) * 1000


def rpm_norm(counts, QC_reads):
    if QC_reads == 0:
        return 0
    return (counts / QC_reads) * 1000000


print("Length of pandas vf dataframe: " + str(len(blast_viral_top_results)))
#print(blast_viral_top_results)

def calculate_read_counts():
    if PE_cov_df is not None:
        test_blast_viral_top_results = pd.merge(
            blast_viral_top_results,
            PE_cov_df[["qseqid", "SR_counts"]],
            on="qseqid",
            how="left",
        )
    else:
        test_blast_viral_top_results = blast_viral_top_results.copy()
        test_blast_viral_top_results["SR_counts"] = 0

    if "SR_counts" in test_blast_viral_top_results.columns:
        test_blast_viral_top_results["SR_counts"] = test_blast_viral_top_results["SR_counts"].fillna(0)

    if LR_cov_df is not None:
        test_blast_viral_top_results2 = pd.merge(
            test_blast_viral_top_results,
            LR_cov_df[["qseqid", "LR_counts"]],
            on="qseqid",
            how="left",
        )
    else:
        test_blast_viral_top_results2 = test_blast_viral_top_results.copy()
        test_blast_viral_top_results2["LR_counts"] = 0

    if "LR_counts" in test_blast_viral_top_results2.columns:
        test_blast_viral_top_results2["LR_counts"] = test_blast_viral_top_results2["LR_counts"].fillna(0)

    test_blast_viral_top_results2.loc[
        test_blast_viral_top_results2["lineage"] == "undefined taxon", "lineage"
    ] = test_blast_viral_top_results2["sallgi"]

    #print(test_blast_viral_top_results2.head())

    test_blast_viral_top_results2.loc[test_blast_viral_top_results2["seq_source"] == 'short_read', 'SR_counts'] = 1
    test_blast_viral_top_results2.loc[test_blast_viral_top_results2["seq_source"] == 'long_read', 'LR_counts'] = 1

    #print(test_blast_viral_top_results2.groupby("lineage"))
    #test_blast_viral_top_results2.to_csv(f"{outdir}/{args.sample_id}_test.csv", sep=",", index=False)

#result = df.groupby('Store').agg({
#    'Sales': 'sum',
#    'Quantity': ['mean', 'max']
#})

    count_df = (
        test_blast_viral_top_results2.groupby("lineage")
        .agg(
            Contig_Count=("qseqid", lambda x: x.str.startswith("NODE").sum()),
            Total_SR_Counts=("SR_counts", "sum"),
            Total_LR_Counts=("LR_counts", "sum"),
        )
        .reset_index()
    )
    print(count_df.head())
    count_df["Total_Read_Count"] = count_df["Total_SR_Counts"] + count_df["Total_LR_Counts"]
    count_df = count_df[count_df["Total_Read_Count"] > 0]
    return count_df


if len(blast_viral_top_results) == 0:
    print("No viral results, skipping read count calculation")
    count_df = pd.DataFrame()
else:
    print("Calculating Read Counts....")
    count_df = calculate_read_counts()


def qkt_function(row):
    if str(row["lineage"]).startswith("Viruses;"):
        family, below_family = get_family(row["lineage"])
    else:
        family, below_family = "N/A", "N/A"

    if family != "N/A":
        mean_genome_family = vf_dict.get(family)
        norm_family = normalization(row["Total_Read_Count"], QC_counts, mean_genome_family) if mean_genome_family else "N/A"
    else:
        norm_family = "N/A"

    rpm = rpm_norm(row["Total_Read_Count"], QC_counts)
    return family, below_family, norm_family, rpm


print("Adding metadata to counts table...")
if len(blast_viral_top_results) != 0 and not count_df.empty:
    count_df[["family", "below_family", "normalized_family", "normalized_rpm"]] = count_df.apply(
        qkt_function, axis=1, result_type="expand"
    )
    col_order = [
        "lineage",
        "family",
        "below_family",
        "Contig_Count",
        "Total_SR_Counts",
        "Total_LR_Counts",
        "Total_Read_Count",
        "normalized_family",
        "normalized_rpm",
    ]
    count_df[col_order].to_csv(f"{outdir}/{args.sample_id}_accurate_read_counts.tsv", sep="\t", index=False)
else:
    with open(os.path.join(outdir, "no_viruses.txt"), "w") as file:
        file.write("No viruses were found.\n")

print("VS parsing finished.")
print(time.ctime())
