#!/usr/bin/env python3

# Exercise full report generation
# Chanel Mosby-Tourtellot Updated:2026-07-22

import pandas as pd
from docx.shared import Inches
from docxtpl import DocxTemplate, InlineImage
import sys, getopt, os.path, pathlib, subprocess, json, glob

argv=sys.argv[1:]
opts,args = getopt.getopt(argv,"s:i:o:p:t:a:")
for opt,arg in opts:
    if opt == '-h':
        print("exercise_report_Full.py -o <output_directory> -s <sample_id> -t <template_file> -i <input directory> -p <project_id> -a <assembler>")
        sys.exit()
    elif opt in ("-o", "--outdir"):
        outdir = arg
    elif opt in ("-s", "--sample_id"):
        sample_id = arg
    elif opt in ("-i", "--input_dir"):
        input_dir = arg
    elif opt in ("-p", "--project_id"):
        project_id = arg
    elif opt in ("-t", "--template_file"):
        template_file = arg
    elif opt in ("-a", "--assembler"):
        assembler = arg



# Prepare qc results
#qc_path = glob.glob(os.path.join(input_dir,f'{sample_id}*/qc_stats/qc_stats_final.tsv'))
#qc_stats_final = pd.read_table(qc_path[0])
try:
    qc_stats_final = pd.read_table(os.path.join(input_dir, f'{sample_id}', 'qc_stats', 'qc_stats_final.tsv'))
except:
    #qc_stats_final = pd.DataFrame(names=['sample_id', 'raw_read_count' 'raw_read_total_length'])
    qc_stats_final = pd.DataFrame()
#print(qc_stats_final)
qc_subset = qc_stats_final.iloc[:, :9]
qc_subset = qc_subset.rename(columns={'Sample': 'sample_id', 'Raw_Num_Reads': 'raw_read_count', 'Sum_Length': 'raw_read_total_length', 'Avg_Length': 'raw_read_avg_length', 'GC': 'raw_read_gc', 'Trimmed_Num_Reads': 'trimmed_read_count', 'Sum_Length.1': 'trimmed_read_total_length', 'Avg_Length.1': 'trimmed_read_avg_length', 'GC.1': 'trimmed_read_gc'})
qc_subset = qc_subset.sort_values(by='sample_id')
#print(qc_subset)

# Adjusting to pull out only the sample for cases where multiple sample_ids are in the same all in one run
qc_subset = qc_subset.loc[qc_subset['sample_id'].str.contains(sample_id),]

#print(sample_id)
raw_read_count = qc_subset.loc[qc_subset['sample_id'] == sample_id, 'raw_read_count'].item()
trimmed_read_count = qc_subset.loc[qc_subset['sample_id'] == sample_id, 'trimmed_read_count'].item()
#raw_read_count = qc_subset.loc['raw_read_count'].item()
#trimmed_read_count = qc_subset.loc['trimmed_read_count'].item()


# Prepare mash results
mash_path = glob.glob(os.path.join(input_dir,'reads_taxonomic_classifier/mash', f'{sample_id}*.mash.txt'))
#full_mash = pd.concat((pd.read_table(f, sep ='\t', names=['identity', 'hash', 'multiplicity', 'pvalue', 'query_id', 'query_name']) for f in mash_path), ignore_index = True)
full_mash = pd.read_table(mash_path[0], names=['identity', 'hash', 'multiplicity', 'pvalue', 'query_id', 'query_name'])
mash = full_mash.iloc[:20]
mash['sample_id'] = sample_id
mash = mash.sort_values(by='sample_id')

# Prepare bracken results
bracken_path = glob.glob(os.path.join(input_dir,'reads_taxonomic_classifier/bracken', f'{sample_id}*.bracken.S.tsv'))
#full_mash = pd.concat((pd.read_table(f, sep ='\t', names=['identity', 'hash', 'multiplicity', 'pvalue', 'query_id', 'query_name']) for f in mash_path), ignore_index = True)
bracken = pd.read_table(bracken_path[0], names=['name', 'taxonomy_id', 'taxonomy_lvl', 'kraken_assigned_reads', 'added_reads', 'new_est_reads', 'fraction_total_reads'])
bracken = bracken.sort_values(by='fraction_total_reads', ascending=False)
bracken = bracken.iloc[1:21]



# Prepare BUSCO results 
#busco_path = os.path.join(input_dir, f'{sample_id}', 'busco', f'{assembler}', f'{sample_id}_{assembler}', f'short_summary.specific.*{sample_id}_{assembler}.json')
#busco_path = glob.glob(busco_path)


#with open(busco_path[0], 'r') as file: 
#    busco = json.load(file)
#busco_results = busco.get('results', {})

# Prepare CheckM results
try:
    checkm = pd.read_table(os.path.join(input_dir, f'{sample_id}', 'assembly_verification', f'{assembler}','checkm/checkm_full_stats.tsv'))
except:
    checkm =pd.DataFrame(columns=['Bin Id', 'Completeness', 'Contamination', 'Genome size (bp)', '# contigs', 'N50 (contigs)', 'Longest contig (bp)'])
checkm_filtered = checkm[['Bin Id', 'Completeness', 'Contamination', 'Genome size (bp)', '# contigs', 'N50 (contigs)', 'Longest contig (bp)']]
checkm_filtered = checkm_filtered.rename(columns={'Bin Id': 'Bin_ID', 'Completeness': 'CheckM_Completeness', 'Contamination': 'CheckM_Contamination', 'Genome size (bp)': 'Total_Size', '# contigs': 'Num_Contigs', 'N50 (contigs)': 'N50', 'Longest contig (bp)': 'Longest_Contig'})
checkm_filtered['sample_id'] = sample_id
checkm_filtered = checkm_filtered.sort_values(by='sample_id')

# Prepare CheckM2 results
try:
    checkm2 = pd.read_table(os.path.join(input_dir, f'{sample_id}', 'assembly_verification', f'{assembler}', 'checkm2/quality_report.tsv'))
except:
    checkm2 = pd.DataFrame(columns=['Name', 'Completeness_Specific', 'Contamination'])
checkm2_filtered = checkm2[['Name', 'Completeness_Specific', 'Contamination']]
checkm2_filtered = checkm2_filtered.rename(columns={'Name': 'Bin_ID', 'Completeness_Specific': 'CheckM2_Completeness', 'Contamination': 'CheckM2_Contamination'})
checkm2_filtered['sample_id'] = sample_id
checkm2_filtered = checkm2_filtered.sort_values(by='sample_id')

completeness_df = checkm_filtered.merge(checkm2_filtered, on='sample_id', how='outer')
#completeness_df = completeness_df[['sample_id', 'CheckM_Completeness',  'CheckM_Contamination', 'CheckM2_Completeness', 'CheckM2_Contamination']]
#completeness_df["BUSCO"] = f"{busco_results.get('Complete percentage')}"


# Prepre BLAST table results
blast = pd.read_table(os.path.join(input_dir, f'{sample_id}','blast', f'{assembler}', f'{sample_id}_{assembler}_core_nt_blastn.tsv'), header=None)
blast.columns = ['contig', 'qlen', 'qstart', 'qend', 'accession', 'hit', 'sstart', 'ssend', 'align_len', 'slen', 'ident', 'mismatch', 'orient', 'gap', 'evalue', 'bitscore', 'score' ]
blast = blast[['contig', 'qlen', 'accession', 'hit', 'ident']]
blast_filtered = blast.drop_duplicates(subset=['contig'], keep='first')
blast_filtered = blast_filtered.sort_values(by='qlen', ascending=False)
blast_filtered = blast_filtered.iloc[:5]
blast_filtered['sample_id'] = sample_id

# Prepare rgi results
try:
    rgi = pd.read_table(os.path.join(input_dir, f'{sample_id}', 'rgi', f'{assembler}', f'{sample_id}_{assembler}_rgi_nuc.txt'))
except:
    rgi = pd.DataFrame(columns=['Contig', 'Best_Hit_ARO', 'Best_Identities', 'Drug Class'])
rgi_filtered = rgi[['Contig', 'Best_Hit_ARO', 'Best_Identities', 'Drug Class']]
rgi_filtered = rgi_filtered.rename(columns={'Contig': 'Contig', 'Best_Hit_ARO': 'Best_Hit_ARO', 'Best_Identities': 'Best_Identities', 'Drug Class': 'Drug_Class'})
rgi_filtered['sample_id'] = sample_id
rgi_filtered = rgi_filtered.sort_values(by='sample_id')

# Prepare amrfinder results
amrfinder_path = glob.glob(os.path.join(input_dir, f'{sample_id}', 'amrfinder', f'{assembler}', f'{sample_id}_{assembler}_amrfinder.tsv'))
try:
    amrfinder = pd.read_table(amrfinder_path[0])
except:
    amrfinder = pd.DataFrame(columns=['Contig id', 'Element symbol', '% Identity to reference', 'Type'])
amr_filtered = amrfinder[amrfinder['Type'] == "VIRULENCE"]
amr_filtered = amr_filtered[['Contig id', 'Element symbol', '% Identity to reference']]
amr_filtered = amr_filtered.rename(columns={'Contig id': 'Contig_id', 'Element symbol': 'Element_symbol', '% Identity to reference': 'Ident'})
amr_filtered['sample_id'] = sample_id
amr_filtered = amr_filtered.sort_values(by='sample_id')


# Import template document
template = DocxTemplate(template_file)

# Create dictionaries for template 
qc_dict = qc_subset.to_dict(orient='records')
mash_dict = mash.to_dict(orient='records')
bracken_dict = bracken.to_dict(orient='records')
checkm_dict = checkm_filtered.to_dict(orient='records')
completeness_dict = completeness_df.to_dict(orient='records')
blast_dict = blast_filtered.to_dict(orient='records')
rgi_dict = rgi_filtered.to_dict(orient='records')
amr_dict = amr_filtered.to_dict(orient='records')
template_dict = {'project_id': project_id, 'sample_id': sample_id, 'assembler': assembler, 'qc': qc_dict, 'mash': mash_dict, 'bracken':bracken_dict, 'checkm': checkm_dict, 'completeness': completeness_dict, 'core_nt_blast': blast_dict, 'rgi': rgi_dict, 'amrfinder': amr_dict, 'raw_read_count': raw_read_count, 'trimmed_read_count': trimmed_read_count}

# Note time stamps for QC and taxonomy are going to be copied along with the updated text from the subreports and are not recaptured here

# # Grab average quality score for trimmed reads
# trimmed_path = glob.glob(os.path.join(input_dir, "*", 'q_stats', f'*_seqkit_trimmedReads.txt'))
# trimmed_qual = pd.concat((pd.read_csv(f, sep ='\t') for f in trimmed_path), ignore_index = True)
# trimmed_avg_quality =  round(trimmed_qual['AvgQual'].mean(),2)
# template_dict['trimmed_avg_quality'] = trimmed_avg_quality

# Set up images for dictionary and template
qc_plot_path = os.path.join(input_dir, f'{sample_id}', 'qc_plots', f'{sample_id}_long_raw_reads_vs_Q1_trimmed_reads.jpeg')
if os.path.exists(qc_plot_path):
    qc_plot = InlineImage(template, f'{qc_plot_path}', Inches(6))
    template_dict['qc_plot'] = qc_plot
else:
    print(f"File '{qc_plot_path}' does not exist.")





# multiqc_pre_path = os.path.join(input_dir, 'multiqc/pretrim/multiqc_plots/png/fastqc_per_base_sequence_quality_plot.png')
# multiqc_pre = InlineImage(template, f'{multiqc_pre_path}', Inches(6))
# template_dict['multiqc_pre'] = multiqc_pre

# multiqc_post_path = os.path.join(input_dir, 'multiqc/post_trim/multiqc_plots/png/fastqc_per_base_sequence_quality_plot.png')
# multiqc_post = InlineImage(template, f'{multiqc_post_path}', Inches(6))
# template_dict['multiqc_post'] = multiqc_post

krona_plot_path = os.path.join(input_dir, 'reads_taxonomic_classifier/kraken2/snapshots', f'{sample_id}_long.krona.snapshot.png')
if os.path.exists(krona_plot_path):
    krona_plot = InlineImage(template, f'{krona_plot_path}', Inches(6))
    template_dict['krona_plot'] = krona_plot
else:
    print(f"File '{krona_plot_path}' does not exist.")

#analysis_timestamp = subprocess.getoutput("date")
#template_dict['analysis_timestamp'] = analysis_timestamp

# Render the report
template.render(template_dict)
savepath = os.path.join(outdir, f'{project_id}_{sample_id}_{assembler}_Report_Full.docx')
print(savepath)
template.save(savepath)
