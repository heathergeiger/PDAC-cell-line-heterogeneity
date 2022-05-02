#!/bin/bash

sample=$1

input=${sample}_non_noisy_cell_node_cnv_calls.bed
ref=PDAC_cancer_genes.bed
output=${sample}_non_noisy_cell_node_cnv_calls_intersect_cancer_genes.txt

intersectBed -wo -a $input -b $ref > $output

#First few lines of output look like this:
#chr1    820000  11440000    chr1:820000-11440000_cell990_0_29   chr1    11106535    11262507    MTOR    155972
#chr1    820000  11560000    chr1:820000-11560000_cell126_2_17   chr1    11106535    11262507    MTOR    155972
#chr1    820000  11580000    chr1:820000-11580000_cell1263_1_17  chr1    11106535    11262507    MTOR    155972

#So field 8 is the gene name.
#Get number of occurences per cell + gene, then number of cells with 1,2,3, and so on intervals for each gene.

awk '{split($4,a,"_");print a[2],$8}' $output | sort | uniq -c | awk '{print $1,$NF}' | sort | uniq -c

#For MP2-B at least, we get the following summary:
#21 genes have only one interval, for every single cell.
#32 genes have a small number of cells with a second interval (up to 15 cells, which is still not a lot out of 2455).
#MTOR has 3 cells with a second interval, and 1 cell with a third interval.

#Since any given gene has so few cells with multiple intervals, just take the interval with the most overlap.
#The part awk '!seen[$4$7]++' is equivalent to doing a sort | uniq specifically to take only the first occurence for each combination of cell + gene.

output_final=${sample}_non_noisy_cell_node_cnv_calls_intersect_cancer_genes_most_overlapping_intervals.txt

awk '{ OFS="\t"}{split($4,a,"_");split(a[2],b,"cell");print $1,$2,$3,b[2],a[3],a[4],$8,$NF}' $output | sort -n -r -k8,8 | awk '!seen[$4$7]++' > $output_final

#First few lines of output_final look like this:
#chr8    9940000 43940000    1302    3   17  TUSC3   349434
#chr8    9900000 16900000    886 4   16  TUSC3   349434
#chr8    9800000 16220000    1183    4   16  TUSC3   349434
#chr8    9760000 35260000    894 2   17  TUSC3   349434
#chr8    9680000 17080000    1122    3   16  TUSC3   349434

#Where first three columns are the coordinates of the interval from node_cnv_calls.bed.
#Column 4 is the cell ID number.
#Column 5 is the copy number value of the interval.
#Column 6 is the event confidence.
#Column 7 is the gene name.
#And column 8 is the number of bp of overlap.
