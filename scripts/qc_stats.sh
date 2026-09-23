nextflow run  /export/nextflow/bdrd/all_in_one_pipeline/qc_stats.nf --samplesheet_csv=/mnt/genomics/common_projects/NAMRU-6/P0031_24_NM/250225_OROV_CPE/qc_stats/samplesheet.csv \
                          --outdir=/mnt/genomics/common_projects/NAMRU-6/P0031_24_NM/250225_OROV_CPE/qc_stats/ \
                          --fastqc_threads=100 \
                          --project_id="250225_OROV_CPE" \
                          -profile cluster_normal



#nextflow run  /export/nextflow/bdrd/all_in_one_pipeline/qc_stats.nf --samplesheet_csv=/mnt/genomics/common_projects/NAVSEA/241203_Data_QC/sheet/samplesheet.csv \
#       			  --outdir=/mnt/genomics/common_projects/NAVSEA/241203_Data_QC/pipeline_output/ \
#			  --fastqc_threads=100 \
#			  --project_id="NAVSEA" \
#			  -profile cluster_normal 
