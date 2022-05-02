library(pheatmap)

samples <- c("MP2B","BXPC3","HPAF2","HPDE","HPNE","MP2C","PANC1")

#Read in sample-level files just created in node_cnv_calls_vs_cancer_genes_wrapper.sh.

for(i in 1:length(samples))
{
    this_sample <- read.csv(paste0(samples[i],"_copy_number_per_cell_cancer_genes_from_node_cnv_calls.csv"),header=TRUE,row.names=1,check.names=FALSE)
    rownames(this_sample) <- paste0(samples[i],"_",rownames(this_sample))
    if(i == 1){combined <- this_sample}
    if(i > 1){combined <- rbind(combined,this_sample)}
}

#Read in per_cell_summary_metrics.csv files per sample, and match up cell names (created to be unique per sample).

for(i in 1:length(samples))
{
    sample <- samples[i]
    this_sample <- read.csv(paste0(sample,"_per_cell_summary_metrics.csv"),header=TRUE,row.names=2,check.names=FALSE)
    rownames(this_sample) <- paste0(samples[i],"_cell",rownames(this_sample))
    if(i == 1){summary <- this_sample}
    if(i > 1){summary <- rbind(summary,this_sample)}
}

summary <- summary[rownames(combined),]

summary <- data.frame(sample = sapply(strsplit(rownames(summary),"_"),"[[",1),
    ploidy = summary$mean_ploidy,
    row.names=rownames(summary))

#Subset a random 1000 cells per sample.

set.seed(1392)

subset_cells <- c()

for(sample in samples)
{
    myind <- which(summary$sample == sample)
    subset_cells <- c(subset_cells,sample(myind,1000,replace=FALSE))
}

#Set max copy number value to no more than 8 for visualization purposes.
#Similary set max mean ploidy value from summary to 4.

combined[combined > 8] <- 8
summary$ploidy[summary$ploidy > 4] <- 4

#Set colorblind-friendly color vector for samples.

mycol <- c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#999999")

sample_colors <- mycol[1:6]
names(sample_colors) <- setdiff(levels(factor(summary$sample)),"HPNE")
sample_colors <- list(sample = sample_colors)

#Output heatmap for HPNE.

pheatmap(combined[summary$sample == "HPNE",],
    show_rownames=FALSE,
    fontsize_col=6,
    legend=FALSE)

#Output heatmap for all other samples.

combined <- combined[subset_cells,]
summary <- summary[subset_cells,]

combined <- combined[summary$sample != "HPNE",]
summary <- summary[summary$sample != "HPNE",]

pheatmap(combined,
    show_rownames=FALSE,
    annotation_row=summary,
    annotation_colors=sample_colors,
    fontsize_col=6)
