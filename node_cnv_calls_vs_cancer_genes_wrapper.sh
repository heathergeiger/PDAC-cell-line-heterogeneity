#!/bin/bash

CNV_data_dir="../download/NYGC_collab/CNV"
summary_file_dir="./scCNV_5Mb_intervals/summary_metrics_files" #Copied all sample-level files per_cell_summary_metrics.csv here as e.g. HPAF2_per_cell_summary_metrics.csv.

sample=$1
bed=${CNV_data_dir}/${sample}/node_cnv_calls.bed
summary=${summary_file_dir}/${sample}_per_cell_summary_metrics.csv

#Outputs a file ${sample}_non_noisy_cell_node_cnv_calls.bed to subdirectory scCNV_cancer_genes below the current working directory.
#This contains a subset of the lines in node_cnv_calls.bed, minus those with is_noisy == 1 or mean_ploidy < 1.5 in per_cell_summary_metrics.csv.
#Also some formatting to make it a true BED file (chr/start/end/name, with name containing all other relevant info like copy number and confidence score).

Rscript node_cnv_calls_subset_non_noisy_and_write_bed.R $sample $bed $summary

#Use bedtools intersectBed to intersect the resulting BED from the above with the cancer gene coordinates.

./node_cnv_calls_vs_cancer_genes.sh $sample

#This file will contain one or more lines for each combination of cell-interval info and gene.
#Convert to wide format cell x gene, with the copy number value per cell per gene.

Rscript node_cnv_calls_vs_cancer_genes_spread.R $sample
