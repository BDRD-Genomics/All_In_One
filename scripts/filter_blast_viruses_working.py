##################################################################################################
#                                                                                                #
#                                The purpose of this script is to                                #
#                   take the blast contigs containing at least 1 viral result and                #
#                           sort the output into ambiguous, phage and viral.                     #
#                                                                                                #
##################################################################################################
import sys
import os
import pandas as pd
from Bio import SeqIO
import argparse
import subprocess as sp
from ete3 import NCBITaxa
import time

parser = argparse.ArgumentParser(
                    prog='Filter Blast Viruses',
				usage='filter_blast_viruses_working.py -t <viral_blast_tab>  -c <{contigs_fasta> -p <PE_covstats> -S <PE_QC_reads> -n <ncbidb -v <virusFamily_genomeSizes.txt>',
				description='Filter putatively viral blast outputs to high confidence viral calls. Outputs viral blast table, phage blast table and ambiguous blast tables. Also outputs viral contigs and semiquantitative read counts.',
				epilog='Parse the blast table from sequence comparative algorithm (MMSeqs,blastn,blastx)')
# add options
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
print(args.se_cov,args.pe_cov,args.lr_cov, args.viral_blasttab, args.contigs, args.outdir, args.qc_count_SR, args.qc_count_LR)
if (args.viral_blasttab == None or args.contigs == None):
        print(parser.usage)
        exit(0)

#virus_genomes_file = "/export/virusseeker/scripts/qkt_scripts/viralFamily_genomeSize.txt"
#ncbi = NCBITaxa(dbfile = "/export/database/taxonomy/taxa.sqlite")
virus_genomes_file= args.viral_genome_sizes
ncbi_file = args.ncbidb
ncbi = NCBITaxa(dbfile = ncbi_file)

filtered_blast_input = args.viral_blasttab
contigs_file =  args.contigs
outdir = args.outdir if args.outdir else "."
if not os.path.exists(outdir):
    os.mkdir(outdir)

SE_cov_file = args.se_cov if args.se_cov else None
PE_cov_file = args.pe_cov if args.pe_cov else None
LR_cov_file = args.se_cov if args.se_cov else None
blacklist_file = args.blacklist if args.blacklist else None
if blacklist_file:
    print("filtering by blacklist")
    with open(args.blacklist, 'r') as bad_acc_file:
        bad_acc_list = bad_acc_file.read().splitlines()


#PE_QC_counts = int(args.qc_count_SR) if args.qc_count_SR else 0
if args.qc_count_SR is None:
    #print("PE_QC_counts")
    #print(PE_QC_counts)
    PE_QC_counts = 0
    #print(PE_QC_counts)
else:
    PE_QC_counts = pd.read_csv(args.qc_count_SR,header=None)
    PE_QC_counts = PE_QC_counts.iloc[0,0]
    print("PE_QC_counts:")
    print(PE_QC_counts)
if args.qc_count_LR is None:
    #print("LR_QC_counts")
    #print(LR_QC_counts)
    LR_QC_counts = 0
    #print(LR_QC_counts)
else:
    LR_QC_counts = pd.read_csv(args.qc_count_LR,header=None)
    LR_QC_counts = LR_QC_counts.iloc[0,0]
    #print("LR_QC_counts:")
    #print(LR_QC_counts)

print("PE_QC_counts")
print(PE_QC_counts)
print("LR_QC_counts")
print(LR_QC_counts)
QC_counts = PE_QC_counts + LR_QC_counts
print("QC Counts:"+str(QC_counts))

if PE_cov_file:
    PE_cov_df = pd.read_csv(PE_cov_file, delimiter="\t", names=["query_ID", "qlen", "SR_counts", "SR_unpaired_counts"])
else:
    PE_cov_df=None
if LR_cov_file:
    LR_cov_df = pd.read_csv(LR_cov_file, delimiter="\t", names=["query_ID", "qlen", "LR_counts", "LR_unpaired_counts"])
else:
    LR_cov_df=None


###########################################################################
# Step 1: filter putative viral results into: Viral, Phage, and Ambiguous #
###########################################################################

#blast_df = pd.read_csv(filtered_blast_input, delimiter="\t", low_memory=True, engine="pyarrow")
blast_df = pd.read_csv(filtered_blast_input, delimiter=",")
phage_keywords = ["phage", "Caudoviricetes", "Microviridae", "unclassified bacterial viruses", "Leviviricetes", "Tubulavirales"]
if blacklist_file:
    blast_df["poor_accessions"] = blast_df["ref_description"].apply(lambda x: True if any(acc in x for acc in bad_acc_list) else False)

#print(bad_acc_list)
#print(blast_df[["query_ID","ref_description", "poor_accessions"]].tail())
sub_blast_df = blast_df.groupby("query_ID").agg(lineages=("lineage", lambda x: list(x)), descriptions=("ref_description", lambda x: list(x)), poor_accession_only=("poor_accessions", lambda x: all(list(x)))).reset_index()
#print(sub_blast_df.head())
#print(len(sub_blast_df["query_ID"].tolist()))
sub_blast_df["cellular_organisms"] = sub_blast_df["lineages"].apply(lambda x: any("cellular organisms" in s for s in x ))
sub_blast_df["undefined_taxon"] = sub_blast_df["lineages"].apply(lambda x: any("undefined taxon" in s for s in x ))
sub_blast_df["phages"] = sub_blast_df["lineages"].apply(lambda x: any(k in s for s in x for k in phage_keywords))
#if blacklist_file:
#    print(bad_acc_list)
#    print(sub_blast_df["descriptions"].apply(lambda desc_list: any(k in s for s in desc_list for k in bad_acc_list)))
#    sub_blast_df["poor_accessions"] = sub_blast_df["descriptions"].apply(lambda desc_list: all(k in s for s in desc_list for k in bad_acc_list))
sub_blast_df["viral_only"]=sub_blast_df.eval("cellular_organisms == False and phages == False and poor_accession_only == False")

#print(sub_blast_df.head())

viral_contig_ids = sub_blast_df.loc[(sub_blast_df["viral_only"] == True) & (sub_blast_df["undefined_taxon"] == False), "query_ID"].tolist()
undefined_contig_ids = sub_blast_df.loc[(sub_blast_df["viral_only"] == True) & (sub_blast_df["undefined_taxon"] == True), "query_ID"].tolist()
phage_contig_ids = sub_blast_df.loc[(sub_blast_df["phages"] == True) & (sub_blast_df["cellular_organisms"] == False), "query_ID"].tolist()
ambiguous_contig_ids = sub_blast_df.loc[sub_blast_df["cellular_organisms"] == True, "query_ID"].tolist()
poor_acc_contig_ids = sub_blast_df.loc[sub_blast_df["poor_accession_only"] == True, "query_ID"].tolist()

print("Undefined contig ID's: "+str(len(undefined_contig_ids)))

# Even with our list of high quality virus calls, there may be some assignments that align exclusively to accessions on a blacklist
# let's review virus calls and filter the results from bad accessions out

##now that we have a list of the undefined taxons let's see if we can suss out which ones go to virus and which go to ambiguous
blast_undefined_temp = blast_df[blast_df["query_ID"].isin(undefined_contig_ids) & blast_df["lineage"].str.contains("undefined taxon")]
blast_undefined_temp["contains_virus"] = blast_undefined_temp["ref_description"].str.contains("virus", case=False)
blast_undefined_temp_sub = blast_undefined_temp.groupby("query_ID").agg(contains_virus=("contains_virus", lambda x: list(x))).reset_index()
blast_undefined_temp_sub["check_virus"] = blast_undefined_temp_sub["contains_virus"].apply(lambda x: all( x )) 
#print(blast_undefined_temp_sub.head())

print("PRE: Viral contig ID's: "+str(len(viral_contig_ids)))
print("PRE: Phage contig ID's: "+str(len(phage_contig_ids)))
print("PRE: Ambiguous contig ID's: "+str(len(ambiguous_contig_ids)))
print("PRE: Undefinied contig ID's: "+str(len(undefined_contig_ids)))


viral_contig_ids += blast_undefined_temp_sub.loc[blast_undefined_temp_sub["check_virus"] == True, "query_ID"].tolist()
ambiguous_contig_ids += blast_undefined_temp_sub.loc[blast_undefined_temp_sub["check_virus"] == False, "query_ID"].tolist()

bad_acc_pattern='|'.join(bad_acc_list)
blast_viral = blast_df[blast_df["query_ID"].isin(viral_contig_ids)]
#print(len(blast_viral))
blast_viral = blast_viral[~blast_viral['ref_description'].str.contains(bad_acc_pattern, case=False, na=False)]
#blast_viral = blast_viral[blast_viral['lineages'].startswith("Viruses")
#print(blast_viral.head())
#print(len(blast_viral))
#print(len(blast_viral))
# remove excessive undefined taxon
dropif = (blast_viral['lineage'] == "undefined taxon")
all_drop_if = blast_viral.groupby("query_ID")['lineage'].transform(lambda x: (x == "undefined taxon").all())
rows_to_keep_mask = ~dropif | ( dropif & all_drop_if)
blast_viral_test = blast_viral[rows_to_keep_mask]
#print(blast_viral_test.head())
#grab first result
viral_mask = ~blast_viral_test.duplicated(subset=["query_ID"],keep='first')
blast_viral_top_results = blast_viral_test[viral_mask].sort_values(by="lineage")
#print(blast_viral_top_results.head(50))
#print(len(blast_viral_top_results))
blast_viral_top_results.to_csv(f"{outdir}/{args.sample_id}_viral_blast_contigs.csv", sep=",", index=False)


blast_phage = blast_df[blast_df["query_ID"].isin(phage_contig_ids)]
phage_mask = ~blast_phage.duplicated(subset=["query_ID"],keep='first')
blast_phage_top_results = blast_phage[phage_mask].sort_values(by="lineage")
blast_phage_top_results.to_csv(f"{outdir}/{args.sample_id}_phage_blast_contigs.csv", sep=",", index=False)

blast_ambiguous = blast_df[blast_df["query_ID"].isin(ambiguous_contig_ids)]
blast_ambiguous.to_csv(f"{outdir}/{args.sample_id}_full_ambiguous_viral_blast_contigs.csv",
                       sep=",", index=False)

blast_ambiguous = blast_df[blast_df["query_ID"].isin(poor_acc_contig_ids)]
blast_ambiguous.to_csv(f"{outdir}/{args.sample_id}_full_bad_accessions_blast_contigs.csv",
                       sep=",", index=False)

print("Viral contig ID's: "+str(len(viral_contig_ids)))
#print("Undefined contig ID's: "+str(len(undefined_contig_ids)))
print("Phage contig ID's: "+str(len(phage_contig_ids)))
print("Ambiguous contig ID's: "+str(len(ambiguous_contig_ids)))
print("Bad accessions only contig ID's: "+str(len(poor_acc_contig_ids)))


###########################################################################
# Step 2: Extract contigs for high confidence viral results               #
###########################################################################

def read_FASTA_data(fastaFile):
        fa_dict = SeqIO.index(fastaFile, "fasta")
        return fa_dict

print("Reading contigs file: " + contigs_file)
seq_dict = read_FASTA_data(contigs_file)
viral_fasta_file = "viral_contigs.fasta"

try:
       outFile = open(outdir+"/"+viral_fasta_file, "w")
except:
       print("Can't create output file : "+ viral_fasta_file)
       sys.exit(1)
print("Writing contigs to output file: "+ viral_fasta_file)
for seq_name in viral_contig_ids:
        outFile.write(">"+seq_name+"\n")
        seq=seq_dict[seq_name].seq
        outFile.write(str(seq)+"\n")

###########################################################################
# Step 3: Calculate semi-quanititative abundance for read counts report   #
###########################################################################

# Converting "Generate_assignment_report.pl" and "VS_Read_Counter.py" into a single script

#First we want to add a column to "blast_viral_top_results" with how many reads mapped to that contig

#check to make sure we have at least one covstats file

if not any([SE_cov_file, PE_cov_file, LR_cov_file]):
    sys.exit()

vf_dict = pd.read_csv(virus_genomes_file, sep="\t", header=0, names=["fam","size"]).set_index('fam')['size'].to_dict()

"""
def get_SR_count(q_id):
    cmd = "grep " + q_id + " " + PE_cov_file
    try:
        PE_line = sp.check_output(cmd, shell=True, encoding="utf-8")
        PE_total_reads = int(PE_line.split("\t")[2])
    except sp.CalledProcessError: # grep non-zero exit status indicates no match
        PE_total_reads=0
        pass
    return PE_total_reads

def get_LR_count(q_id):
    cmd = "grep " + q_id + " " + LR_cov_file
    try:
        LR_line = sp.check_output(cmd, shell=True, encoding="utf-8")
        LR_total_reads = int(SE_line.split("\t")[2])
    except sp.CalledProcessError: # grep non-zero exit status indicates no match
        LR_total_reads=0
        pass
    return LR_total_reads
"""

def get_family(full_taxid):
    parsed = full_taxid.split(';')
    parsed.pop()
    
    taxid=ncbi.get_name_translator([parsed[-1]])
    lineage=ncbi.get_lineage(taxid.get(parsed[-1])[0])
    rank = ncbi.get_rank(lineage)

    family_key=list(rank.keys())
    family_value=list(rank.values())
    try:
        family_idx=family_value.index("family")
        family_dict=ncbi.get_taxid_translator([family_key[family_idx]])
        family=family_dict.get(family_key[family_idx])
        #print(family)
        try:
            taxid_family_index = parsed.index(family)
            sub_family = ';'.join(parsed[taxid_family_index+1:])
            return family, sub_family
        except ValueError:
            #sub_family_arr.append('N/A')
            return family, "N/A"
    except:
         return "N/A", "N/A"

def normalization(counts, QC_reads, genome_length):
    return ((counts/QC_reads)/genome_length)*1000
def rpm_norm(counts, QC_reads):
    return (counts/QC_reads)*1000000

virus_dict={}

print("Length of pandas vf dataframe: " + str(len(blast_viral_top_results)))

def calculate_read_counts():
    if PE_cov_df is not None:
        test_blast_viral_top_results = pd.merge(blast_viral_top_results, PE_cov_df[["query_ID", "SR_counts"]])
    else:
        test_blast_viral_top_results = blast_viral_top_results
        test_blast_viral_top_results["SR_counts"] = 0

    if LR_cov_df is not None:
        test_blast_viral_top_results2 = pd.merge(test_blast_viral_top_results, LR_cov_df[["query_ID", "LR_counts"]])
    else:
        test_blast_viral_top_results2 = test_blast_viral_top_results
        test_blast_viral_top_results2["LR_counts"] = 0

    test_blast_viral_top_results2.loc[test_blast_viral_top_results2["lineage"] == "undefined taxon", "lineage"] = test_blast_viral_top_results2["ref_description"]
    #print(test_blast_viral_top_results2.head())
    count_df = test_blast_viral_top_results2.groupby("lineage").agg(
        Contig_Count=("lineage", "count"),
        Total_SR_Counts=("SR_counts", "sum"),
        Total_LR_Counts=("LR_counts", "sum")
        ).reset_index()
    count_df["Total_Read_Count"] = count_df["Total_SR_Counts"] + count_df["Total_LR_Counts"]
    count_df = count_df[count_df["Total_Read_Count"] > 0]
    #print(count_df)
    return count_df

if len(blast_viral_top_results) == 0:
    print("No viral results, skipping read count calculation")
else:
    print("Calculating Read Counts....")
    count_df = calculate_read_counts()

#print(count_df.head())
#count_df["Total_Read_Count"] = count_df["Total_SR_Counts"] + count_df["Total_LR_Counts"]

# now we must get family/sub-famly and normalized read counts

def qkt_function(row):
    #print(row["lineage"])
    if row["lineage"].startswith("Viruses;"):
        family, below_family = get_family(row["lineage"])
    else:
        family, below_family = "N/A", "N/A"
    if family != "N/A":
        try:
            mean_genome_family = vf_dict[family]
        except KeyError:
            mean_genome_family=None
            pass
        if mean_genome_family:
            norm_family = normalization(row["Total_Read_Count"], QC_counts, mean_genome_family)
        else:
            norm_family = "N/A"
    else:
        norm_family = "N/A"
    rpm = rpm_norm(row["Total_Read_Count"], QC_counts)
    #print(family, below_family, norm_family, rpm)
    return family, below_family, norm_family, rpm

print("Adding metadata to counts table...")
if len(blast_viral_top_results) != 0:
    count_df[["family", "below_family", "normalized_family", "normalized_rpm"]] = count_df.apply(qkt_function, axis=1, result_type="expand")
    col_order = ["lineage", "family", "below_family", "Contig_Count", "Total_SR_Counts", "Total_LR_Counts", "Total_Read_Count", "normalized_family", "normalized_rpm"]
    count_df[col_order].to_csv(f"{outdir}/{args.sample_id}_accurate_read_counts.tsv", sep="\t", index=False)
else:
    with open('no_viruses.txt', 'w') as file:
        file.write("No viruses were found.\n")
#col_order = ["lineage", "family", "below_family", "Contig_Count", "Total_SR_Counts", "Total_LR_Counts", "Total_Read_Count", "normalized_family", "normalized_rpm"]
#count_df[col_order].to_csv(outdir+"/accurate_read_counts.tsv", sep="\t", index=False)

print("VS parsing finished.")
print(time.ctime())

#print(count_df.head())
