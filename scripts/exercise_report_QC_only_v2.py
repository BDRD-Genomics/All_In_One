#!/usr/bin/env python3

# Exercise QC subreport generation
# Chanel Mosby-Tourtellot Updated:2026-05-14

import pandas as pd
from docx.shared import Inches
from docxtpl import DocxTemplate, InlineImage
import sys, getopt, os.path, pathlib, subprocess, glob

argv=sys.argv[1:]
opts,args = getopt.getopt(argv,"s:i:o:p:t:")
for opt,arg in opts:
    if opt == '-h':
        print("exercise_report_QC_only.py -o <output_directory> -s <sample_id> -t <template_file> -i <input directory> -p <project name>")
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

env = os.environ.copy()
env["TZ"]="EST5EDT"

# Prepare qc results
qc_stats_final = pd.read_table(os.path.join(input_dir, 'qc_stats', 'qc_stats_final.tsv'))
qc_subset = qc_stats_final.iloc[:, :9]
qc_subset = qc_subset.rename(columns={'Sample': 'sample_id', 'Raw_Num_Reads': 'raw_read_count', 'Sum_Length': 'raw_read_total_length', 'Avg_Length': 'raw_read_avg_length', 'GC': 'raw_read_gc', 'Trimmed_Num_Reads': 'trimmed_read_count', 'Sum_Length.1': 'trimmed_read_total_length', 'Avg_Length.1': 'trimmed_read_avg_length', 'GC.1': 'trimmed_read_gc'})
qc_subset = qc_subset.sort_values(by='sample_id')

# Import template document
template = DocxTemplate(template_file)

# Create dictionaries for template 
qc_dict = qc_subset.to_dict(orient='records')
template_dict = {'project_id': project_id, 'sample_id': sample_id, 'qc': qc_dict}

# Date and time capture
qc_timestamp = subprocess.getoutput ('date +"%A %B %d, %Y, %T %Z"')
template_dict['qc_timestamp'] = qc_timestamp

# Grab average qual # Updated for multiple samples
trimmed_path = glob.glob(os.path.join(input_dir, "*", 'q_stats', f'*_seqkit_trimmedReads.txt'))
trimmed_qual = pd.concat((pd.read_csv(f, sep ='\t') for f in trimmed_path), ignore_index = True)
trimmed_avg_quality =  round(trimmed_qual['AvgQual'].mean(),2)
template_dict['trimmed_avg_quality'] = trimmed_avg_quality


# Set up images for dictionary and template
qc_plot_path = os.path.join(input_dir, 'qc_plots', f'*raw_reads_vs_trimmed_reads.jpeg')
qc_plot_path = glob.glob(qc_plot_path)
qc_plot = InlineImage(template, f'{qc_plot_path[0]}', Inches(6))
template_dict['qc_plot'] = qc_plot

multiqc_pre_path = os.path.join(input_dir, 'multiqc/pretrim/multiqc_plots/png/fastqc_per_base_sequence_quality_plot.png')
multiqc_pre = InlineImage(template, f'{multiqc_pre_path}', Inches(6))
template_dict['multiqc_pre'] = multiqc_pre

multiqc_post_path = os.path.join(input_dir, 'multiqc/post_trim/multiqc_plots/png/fastqc_per_base_sequence_quality_plot.png')
multiqc_post = InlineImage(template, f'{multiqc_post_path}', Inches(6))
template_dict['multiqc_post'] = multiqc_post

# Render the report
template.render(template_dict)
savepath = os.path.join(outdir, f'{project_id}_{sample_id}_Exercise_Report_QC.docx')
template.save(savepath)

