library(tidyr)

args <- commandArgs(trailingOnly=TRUE)
sample <- args[1]

info <- read.table(paste0(sample,"_non_noisy_cell_node_cnv_calls_intersect_cancer_genes_most_overlapping_intervals.txt"),header=FALSE)
colnames(info) <- c("chr","start","end","cell","copies","conf","gene","overlap")
info$cell <- paste0("cell",info$cell)
info <- info[,c("cell","copies","gene")]
info <- spread(info,key="gene",value="copies")
rownames(info) <- info$cell
info <- info[,2:ncol(info)]

write.csv(info,file=paste0(sample,"_copy_number_per_cell_cancer_genes_from_node_cnv_calls.csv"))
