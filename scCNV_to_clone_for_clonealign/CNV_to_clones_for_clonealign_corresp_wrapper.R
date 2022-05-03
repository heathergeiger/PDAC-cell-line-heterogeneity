library(tidyr);library(dplyr);library(pheatmap)

#The code below is based on looking at clones within MP2-B.
#Will eventually convert this to an Rmarkdown so can see the logic at each step once visualizations were used to make decisions.

sample="MP2B"

info <- read.table(paste0(sample,"_intersect_mappable_regions_one_copy_number_per_cell_region_combination.txt"),header=FALSE)
colnames(info) <- c("cell","region","copy_num")
info <- spread(info,key="region",value="copy_num")
rownames(info) <- info$cell
info <- info[,2:ncol(info)]
info <- info[,grep('chrX',colnames(info),invert=TRUE)]
info <- info[,grep('chrY',colnames(info),invert=TRUE)]
chr_per_region <- sapply(strsplit(colnames(info),":"),"[[",1)
region_metadata <- data.frame(chr = sapply(strsplit(colnames(info),":"),"[[",1),
    start_and_end = sapply(strsplit(colnames(info),":"),"[[",2))
region_metadata$start <- as.numeric(sapply(strsplit(region_metadata$start_and_end,"-"),"[[",1))
region_metadata$end <- as.numeric(sapply(strsplit(region_metadata$start_and_end,"-"),"[[",2))
region_metadata$length <- region_metadata$end - region_metadata$start
region_metadata$chr_numeric <- as.numeric(sapply(strsplit(region_metadata$chr,"chr"),"[[",2))
myorder <- order(region_metadata$chr_numeric,region_metadata$start,region_metadata$end)
region_metadata <- region_metadata[myorder,]
info <- info[,myorder]
mean_per_chr <- matrix(NA,nrow=nrow(info),ncol=22,dimnames=list(rownames(info),paste0("chr",1:22)))
for(i in 1:22)
{
    chr <- paste0("chr",i)
    chr_indices <- which(region_metadata$chr == chr)
    info_this_chr <- info[,chr_indices]
    region_metadata_this_chr <- region_metadata[chr_indices,]
    region_metadata_this_chr$prop_chr <- region_metadata_this_chr$length/sum(region_metadata_this_chr$length)
    info_this_chr <- sweep(info_this_chr,2,region_metadata_this_chr$prop_chr,FUN="*")
    mean_per_chr[,i] <- apply(info_this_chr,1,sum)
}

median_of_means_per_chr <- apply(mean_per_chr,1,median)

min_median_of_means_per_chr <- mean(median_of_means_per_chr) - sd(median_of_means_per_chr)
max_median_of_means_per_chr <- mean(median_of_means_per_chr) + sd(median_of_means_per_chr)

info <- info[median_of_means_per_chr >= min_median_of_means_per_chr & median_of_means_per_chr <= max_median_of_means_per_chr,]
mean_per_chr <- mean_per_chr[rownames(info),]

mean_per_chr[mean_per_chr > 6] <- 6
pheatmap(mean_per_chr,cluster_cols=FALSE,show_rownames=FALSE)

#The following chr seem to warrant further review:
#chr1,3,7,10,12,21

info_trunc <- info
info_trunc[info_trunc > 6] <- 6

rownames(region_metadata) <- colnames(info)

region_indices <- which(region_metadata$chr_numeric %in% c(1,3,7,10,12,21))

pheatmap(info_trunc[,region_indices],
    show_rownames=FALSE,cluster_cols=FALSE,show_colnames=FALSE,
    annotation=region_metadata[,c("chr","start","end","length")])

pheatmap(info_trunc[,region_metadata$chr_numeric == 1],
    show_rownames=FALSE,cluster_cols=FALSE,
    annotation=region_metadata[,c("start","end","length")])

pheatmap(info_trunc[,region_metadata$chr_numeric == 3],
    show_rownames=FALSE,cluster_cols=FALSE,
    annotation=region_metadata[,c("start","end","length")])

pheatmap(info_trunc[,region_metadata$chr_numeric == 7],
    show_rownames=FALSE,cluster_cols=FALSE,
    annotation=region_metadata[,c("start","end","length")])

pheatmap(info_trunc[,region_metadata$chr_numeric == 12,],
    show_rownames=FALSE,cluster_cols=FALSE,
    annotation=region_metadata[,c("start","end","length")])

#Remove chr1p, starting at coordinate 149,880,000  to get region where copy number 3 vs. 3.
#Remove chr3p, which is mostly consistent between cells.
#Remove chr7q which is mostly the same between cells, or has lower frequency amplifications that are not consistent with the division into 2-3 clones we otherwise seem to be seeing.
#Also remove chr12p and parts of chr12q that are overwhelmingly copy number=2 rather than a mix of 2 and 3.

region_indices <- which(region_metadata$chr_numeric %in% c(1,3,7,10,12,21))

region_indices <- setdiff(region_indices,which(region_metadata$chr_numeric == 1 & region_metadata$start < 149e6))
region_indices <- setdiff(region_indices,which(region_metadata$chr_numeric == 3 & region_metadata$start < 90e6))
region_indices <- setdiff(region_indices,which(region_metadata$chr_numeric == 7 & region_metadata$start > 60e6))
region_indices <- setdiff(region_indices,which(region_metadata$chr_numeric == 12 & region_metadata$start < 80e6))

pheatmap(info_trunc[,region_indices],
    show_rownames=FALSE,cluster_cols=FALSE,show_colnames=FALSE,
    annotation=region_metadata[,c("chr","start","end","length")])

#There is a clear subset of cells with chr7p copy number=3 vs. 4.
#The cells with chr7p copy number=4 also have chr10 copy number=2 and chr12q subset copy number=3.
#Meanwhile, at least an easily defined subset of cells with chr7p copy number=3 also have chr10 copy number=3 and chr12q subset copy number=2.

#While we also see some correspondence between chr7/10/12 copy number values and those in chr1,3, and 21, I later tried assorted variations of including these, and did not get as clean a result for clonealign.
#Take a final look at mappable regions vs. copy number for chr7, 10, and 12.

region_indices <- setdiff(region_indices,which(region_metadata$chr_numeric %in% c(1,3,21)))

pheatmap(info_trunc[,region_indices],
    show_rownames=FALSE,cluster_cols=FALSE,show_colnames=TRUE,
    annotation=region_metadata[,c("chr","start","end","length")])

#Looks like chr7 up to 4,580,000 is amplified less frequently than the rest of chr7p.
#For chr12, amplification seems to be less frequent in 80260000-83660000, 83680000-87740000, 87760000-90520000, and 127100000-131200000.

region_names <- rownames(region_metadata)[region_indices]
to_remove <- match(c("chr7:240000-4580000","chr12:80260000-83660000","chr12:83680000-87740000","chr12:87760000-90520000","chr12:127100000-131200000"),region_names)
region_indices <- region_indices[setdiff(1:length(region_indices),to_remove)]

#Let's define clone1 as the one with chr7p=4, chr10=2,chr12=3 and clone2 as the one with chr7p=3,chr10=3,chr12=2.
#Output in bed format so can intersect with gene coordinates.

bed_output <- data.frame(region_metadata[region_indices,c("chr","start","end")],
    clone1 = rep(NA,times=length(region_indices)),
    clone2 = rep(NA,times=length(region_indices)))

bed_output$clone1[bed_output$chr == "chr7"] <- 4
bed_output$clone2[bed_output$chr == "chr7"] <- 3
bed_output$clone1[bed_output$chr == "chr10"] <- 2
bed_output$clone2[bed_output$chr == "chr10"] <- 3
bed_output$clone1[bed_output$chr == "chr12"] <- 3
bed_output$clone2[bed_output$chr == "chr12"] <- 2

bed_output$name <- paste0(bed_output$clone1,"_",bed_output$clone2)

write.table(bed_output[,c("chr","start","end","name")],file=paste0(sample,"_segment_by_clone.bed"),row.names=FALSE,quote=FALSE,sep="\t",col.names=FALSE)
