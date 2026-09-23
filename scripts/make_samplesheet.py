#!/bin/python

import os
import glob
import argparse
import re
import pandas as pd
import numpy as np

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="VS_MD_parser",
        description="Parse blast table",
        usage="python3  -i <blast_input> -n <ncbi_database> -v <vhunter_database>",
        epilog="Parse the blast table from sequence comparative algorithm (MMSeqs, blastn, blastx).",
    )
    parser.add_argument("-r", "--reads", required=False, help="Path to directory with trimmed reads")
    parser.add_argument("-c", "--contigs", required=False, help="Path to directory with assembly files")
    parser.add_argument("-d", "--data_dir", required=False, help="Path to directory containing both fastq and assembly files")
    parser.add_argument("-o", "--outdir", required=False, default=".")
    parser.add_argument("-a", "--assembler", required=True)
    return parser.parse_args()

def find_reads(dir):
    #replace local path to absolute path
    #dir = os.path.abspath(dir) if dir.startswith(".") else dir
    dir = os.path.abspath(dir)
    #grab all fastq files
    read_list = glob.glob(dir+"/**/*.fastq.gz", recursive=True)
    #print(read_list)
    return read_list

def find_contigs(dir):
    #replace local path to absolute path
    dir = os.path.abspath(dir)
    #grab all assembly files
    contig_list = glob.glob(dir+"/**/*.fasta", recursive=True)+glob.glob(dir+"/**/*.fa", recursive=True)+glob.glob(dir+"/**/*.fna", recursive=True)
    return contig_list

def get_sample_name(full_file_list):
    #regex_file_patterns = [r".fasta", r".fna", r".fa"]
    #_S7_fastp_R1.fastq.gz
    regex_file_patterns = [
        r"_metaspades_contigs.fasta",
        r"_S\d*_R1(.*)\.fastq\.gz", r"_S\d*_R2(.*)\.fastq\.gz", 
        r"_R1(.*)\.fastq\.gz", r"_R2(.*)\.fastq\.gz",
        #r"_S\d*_fastp_R1\.fastq.gz", r"_S\d*_fastp_R2\.fastq.gz",
        r"\.LR\.trimmed\.fastq\.gz",
        r"_fastp_R1\.fastq\.gz", r"_fastp_R2\.fastq\.gz",
        r"_host_removed_sr_R1\.fastq\.gz", r"_host_removed_sr_R2\.fastq\.gz",
        r"_CLC\.fasta",
        r"_host_removed_LR\.fastq\.gz", 
        r"_S\d*_L001_R1(.*)\.fastq\.gz", r"_S\d*_L001_R2(.*)\.fastq\.gz",
        r"\.fastq\.gz",
        r"_contigs.fasta",
        r"\.assembly\.fasta",
        r"\.job\d+\.polished_assembly\.fasta",
        r".fasta", r".fna", r".fa"
        r"_raven",r"_dragonflye",r"_metspades", r"_spades"
        ]
    combined_patterns= "|".join(regex_file_patterns)
    #file_list = [fname.rsplit("/",1)[1] for fname in full_file_list]
    file_dict = {re.sub(combined_patterns, "", f.rsplit("/",1)[1]): f for f in full_file_list}
    return file_dict

def main():
    args = parse_args()
    print(args)

    if (args.reads is None or args.contigs is None) and args.data_dir is None:
        print("Missing required args.")
        print("  Required: -r/--reads and -c/--contigs OR -d/--data-dir")
        return 1

    reads_dir = args.reads if args.reads else None
    contigs_dir = args.contigs if args.contigs else None
    data_dir = args.data_dir if args.data_dir else None
    outdir = args.outdir if args.outdir else "./"

    reads_files = find_reads(reads_dir) if reads_dir else find_reads(data_dir)
    contig_files = find_contigs(contigs_dir) if contigs_dir else find_contigs(data_dir)

    fastq_R1_files = [fastq for fastq in reads_files if re.search("R1", fastq)]
    fastq_R2_files = [fastq for fastq in reads_files if re.search("R2", fastq)]
    lr_files = [fastq for fastq in reads_files if fastq not in fastq_R1_files and fastq not in fastq_R2_files]

    lr_dict = get_sample_name(lr_files)
    fastq_R1_dict = get_sample_name(fastq_R1_files)
    fastq_R2_dict = get_sample_name(fastq_R2_files)
    contig_dict = get_sample_name(contig_files)
    df = pd.DataFrame([fastq_R1_dict,fastq_R2_dict,lr_dict, contig_dict]).T
    df = df.reset_index()
    df.columns = ["sample_id","fastq_1", "fastq_2", "long_read", "contigs"]
    mode_conditions = [df["fastq_1"].notna() & df["fastq_2"].notna() & df["long_read"].notna(),
                       df["fastq_1"].notna() & df["fastq_2"].notna(),
                       df["long_read"].notna()]
    mode_values = ["hybrid", "short", "long"]
    df["mode"] = np.select(mode_conditions, mode_values, default="ERROR")
    df["assembler"] = args.assembler
    if "ERROR" in df["mode"].values:
        print("You have unpaired samples. Removing from samplesheet...")
        df = df[df["mode"] != "ERROR"]
    df = df.fillna("No_Read")
    print("Saving samplesheet to "+outdir+"/samplesheet.csv")
    df.to_csv(outdir+"/samplesheet.csv", index=False)

if __name__ == "__main__":
    raise SystemExit(main())
