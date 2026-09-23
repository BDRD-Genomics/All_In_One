#!/usr/bin/env python3

# Exercise Taxonomy subreport generation
# Chanel Mosby-Tourtellot Updated:2026-05-14

import pandas as pd
from docx.shared import Inches
from docxtpl import DocxTemplate, InlineImage
import sys, getopt, os.path, pathlib, subprocess, glob

argv=sys.argv[1:]
opts,args = getopt.getopt(argv,"s:i:o:p:t:")
for opt,arg in opts:
    if opt == '-h':
        print("exercise_report_Taxonomy_only.py -o <output_directory> -s <sample_id> -t <template_file> -i <input directory> -p <project_id>")
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

# Prepare mash results
mash_path = glob.glob(os.path.join(input_dir,'reads_taxonomic_classifier/mash', f'{sample_id}*.mash.txt'))
#full_mash = pd.concat((pd.read_table(f, sep ='\t', names=['identity', 'hash', 'multiplicity', 'pvalue', 'query_id', 'query_name']) for f in mash_path), ignore_index = True)
full_mash = pd.read_table(mash_path[0], names=['identity', 'hash', 'multiplicity', 'pvalue', 'query_id', 'query_name'])
mash = full_mash.iloc[:5]
mash['sample_id'] = sample_id
mash = mash.sort_values(by='sample_id')

# Import template document
template = DocxTemplate(template_file)

# Create dictionaries for template 
mash_dict = mash.to_dict(orient='records')
template_dict = {'project_id': project_id, 'sample_id': sample_id, 'mash': mash_dict}

# Date and time capture
taxonomy_timestamp = subprocess.getoutput ('date +"%A %B %d, %Y, %T %Z"')
template_dict['taxonomy_timestamp'] = taxonomy_timestamp

# Set up images for dictionary and template
krona_plot_path = os.path.join(input_dir, 'reads_taxonomic_classifier/kraken2/snapshots', f'{sample_id}_long.krona.snapshot.png')
if os.path.exists(krona_plot_path):
    krona_plot = InlineImage(template, f'{krona_plot_path}', Inches(6))
    template_dict['krona_plot'] = krona_plot
else:
    print(f"File '{krona_plot_path}' does not exist.")

# Render the report
template.render(template_dict)
savepath = os.path.join(outdir, f'{project_id}_{sample_id}_Exercise_Report_Taxonomy.docx')
template.save(savepath)
