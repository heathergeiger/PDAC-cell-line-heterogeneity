#!/bin/bash

export PATH=/nfs/sw/bedtools/bedtools-2.17.0/bin:$PATH #Set path to intersectBed

sample=$1

#Assuming already ran this from other code-subdirectory:
#Rscript scripts/node_cnv_calls_subset_non_noisy_and_write_bed.R $sample $bed $summary

#Creates file ${sample}_non_noisy_cell_node_cnv_calls.bed.
#File ${sample}_mappable_regions.bed should just be a copy of 10X scCNV software output mappable_regions.bed.

grep -v '^#' ${sample}_mappable_regions.bed | awk '($3 - $2) >= 2000000' > ${sample}_mappable_regions_2Mb.bed
intersectBed -wo -a ${sample}_non_noisy_cell_node_cnv_calls.bed -b ${sample}_mappable_regions_2Mb.bed > ${sample}_intersect_mappable_regions.txt

#First few lines of intersectBed result:
#chr1    820000  2800000 chr1:820000-2800000_cell460_1_15    chr1    2780000 12860000    20000
#chr1    820000  2800000 chr1:820000-2800000_cell754_4_15    chr1    2780000 12860000    20000
#chr1    820000  2820000 chr1:820000-2820000_cell1581_4_15   chr1    2780000 12860000    40000

#Get the total amount of overlap for each combination of cell + copy number + region.
#This is the first part here combined with the awk '{sum[$1]+=$2} END {for(i in sum)print i"\t"sum[i]}'. 
#Then, take only the first occurence of each combination of cell + region after reverse sort by amount of overlap (sort -n -r -k2,2).
#The "take only the first occurence" part here is the awk '!seen[$1]++'.
#Finally, re-separate cell, region, and copy number with the "tr" command.

awk '{ OFS="\t"}{split($4,a,"_");cellnum=a[2];copynum=a[3];region=$5":"$6"-"$7;overlap=$NF;print cellnum"_"region"_"copynum,overlap}' ${sample}_intersect_mappable_regions.txt | \
 awk '{sum[$1]+=$2} END {for(i in sum)print i"\t"sum[i]}' | \
 sort -n -r -k2,2 | \
 awk '{ OFS="\t"}{split($1,a,"_");print a[1]"_"a[2],a[3]}' | \
 awk '!seen[$1]++' \
 tr '_' '\t' > ${sample}_intersect_mappable_regions_one_copy_number_per_cell_region_combination.txt
 
#Get number of genes overlapping each mappable region.
#Will use this in the R script run immediately after this.

intersectBed -wo -a ${sample}_mappable_regions_2Mb.bed -b Ensembl_84_minus_non_canonical_and_chrM_unique_names_select_biotypes.bed | \
 awk '{print $1":"$2"-"$3}' | \
 awk '{seen[$1]++} END {for(i in seen)print i"\t"seen[i]}' > ${sample}_mappable_regions_2Mb_num_intersecting_genes.txt
 
#Use an R script to convert from long to wide format (cell x region copy number matrix), and then go from here to define clones all within that script.
#Will eventually make this as a markdown file so can show plots.
 
Rscript scripts/CNV_to_clones_for_clonealign_corresp_wrapper.R
 
#End result of this script looks like this (for MP2-B):
#chr7	7480000	22520000	4_3
#chr7	22540000	25040000	4_3
#chr7	25620000	29580000	4_3
#chr7	34820000	44000000	4_3
#chr7	44040000	46820000	4_3
#chr7	46840000	48840000	4_3
#chr7	49700000	54220000	4_3
#chr10	760000	4960000	2_3
#chr10	6380000	11720000	2_3
#...

#This is a BED file with the name field containing clone1-copy-num_clone2-copy-num.
#Intersect this BED file with gene coordinates. Here using Ensembl 84, excluding certain chromosomes and gene biotypes to match the GTF used for scRNA, and giving a unique gene symbol per gene ID.

#Afer the intersect, awk command to convert to CSV gene x clone output.
#The awk '!seen[$0]++' uniques for the few cases where the gene overlaps multiple segments.

gene_coord_bed=Ensembl_84_minus_non_canonical_and_chrM_unique_names_select_biotypes.bed

echo "Gene,clone1,clone2" > ${sample}_gene_by_clone.csv
intersectBed -wo -a ${sample}_segment_by_clone.bed -b $gene_coord_bed | \
 awk '{ OFS=","}{split($4,a,"_");print $8,a[1],a[2]}' | \
 awk '!seen[$0]++' >> ${sample}_gene_by_clone.csv
