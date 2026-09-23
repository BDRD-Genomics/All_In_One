#!/usr/bin/env python3

# Exercise Taxonomy subreport generation
# Chanel Mosby-Tourtellot Updated:2026-05-14

import pandas as pd
from docx.shared import Inches
from docxtpl import DocxTemplate, InlineImage
import sys, getopt, os.path, pathlib, subprocess

argv=sys.argv[1:]
opts,args = getopt.getopt(argv,"s:i:o:n:t:")
for opt,arg in opts:
    if opt == '-h':
        print("exercise_report_Taxonomy_only.py -o <output_directory> -s <dataset_id> -t <template_file> -i <input directory> -n <sub name>")
        sys.exit()
    elif opt in ("-o", "--outdir"):
        outdir = arg
    elif opt in ("-s", "--dataset_id"):
        dataset_id = arg
    elif opt in ("-i", "--input_dir"):
        input_dir = arg
    elif opt in ("-n", "--sub_name"):
        sub_name = arg
    elif opt in ("-t", "--template_file"):
        template_file = arg



# Prepare mash results
full_mash = pd.read_table(os.path.join(input_dir, 'reads_taxonomic_classifier/mash', f'{dataset_id}_{sub_name}_long.mash.txt'), names=['identity', 'hash', 'multiplicity', 'pvalue', 'query_id', 'query_name'])
mash = full_mash.iloc[:5]
mash['sample_id'] = sub_name
mash = mash.sort_values(by='sample_id')

# Import template document
template = DocxTemplate(template_file)

# Create dictionaries for template 
mash_dict = mash.to_dict(orient='records')
template_dict = {'dataset_id': dataset_id, 'sub_name': sub_name, 'mash': mash_dict}

# Date and time capture
taxonomy_timestamp = subprocess.getoutput("date")
template_dict['taxonomy_timestamp'] = taxonomy_timestamp

# Set up images for dictionary and template
krona_plot_path = os.path.join(input_dir, 'reads_taxonomic_classifier/kraken2', f'{dataset_id}_{sub_name}_long_krona.png')
if os.path.exists(krona_plot_path):
    krona_plot = InlineImage(template, f'{krona_plot_path}', Inches(6))
    template_dict['krona_plot'] = krona_plot
else:
    print(f"File '{krona_plot_path}' does not exist.")

# Render the report
template.render(template_dict)
savepath = os.path.join(outdir, f'{dataset_id}_{sub_name}_Exercise_Report_Taxonomy.docx')
template.save(savepath)
