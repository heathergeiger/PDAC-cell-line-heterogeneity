#' Load libraries.
#' 
## -----------------------------------------------------------------------------
library(tidyr)
library(dplyr)
library(pheatmap)

#' 
#' Set sample name variable.
#' 
#' Here for demonstration purposes, we will be looking at MP2B.
#' 
#' Then, read in the text file created by scCNV_to_clones_for_clonealign_wrapper.sh.
#' 
## -----------------------------------------------------------------------------
sample="MP2B"

info <- read.table(paste0(sample,"_intersect_mappable_regions_one_copy_number_per_cell_region_combination.txt.gz"),header=FALSE,stringsAsFactors=FALSE)

#' 
#' Do some basic formatting including set column names, convert from long to wide (cell x region), remove non-autosomes.
#' 
## -----------------------------------------------------------------------------
colnames(info) <- c("cell","region","copy_num")
info <- spread(info,key="region",value="copy_num")
rownames(info) <- info$cell
info <- info[,2:ncol(info)]
info <- info[,grep('chrX',colnames(info),invert=TRUE)]
info <- info[,grep('chrY',colnames(info),invert=TRUE)]

#' 
#' Also get the chr, start, and end info for each region in a separate table.
#' 
#' Region names are currently eg. "chr1:1-20000", want to separate the chr1, 1, and 20000.
#' 
#' Also add a field for length (end - start), which we will use to weight regions for the mean per chromosome.
#' 
## -----------------------------------------------------------------------------
chr_per_region <- sapply(strsplit(colnames(info),":"),"[[",1)
start_and_end_per_region <- sapply(strsplit(colnames(info),":"),"[[",2)
start_per_region <- as.numeric(sapply(strsplit(start_and_end_per_region,"-"),"[[",1))
end_per_region <- as.numeric(sapply(strsplit(start_and_end_per_region,"-"),"[[",2))
chr_per_region_numeric <- as.numeric(sapply(strsplit(chr_per_region,"chr"),"[[",2))
region_metadata <- data.frame(chr = chr_per_region_numeric,
	start = start_per_region,  
	end = end_per_region,
	row.names=colnames(info),
	stringsAsFactors=FALSE)
region_metadata$length <- region_metadata$end - region_metadata$start

#' 
#' Get number of genes per region and add this to table.
#' 
## -----------------------------------------------------------------------------
genes_per_region <- read.table(paste0(sample,"_mappable_regions_2Mb_num_intersecting_genes.txt"),
                               header=FALSE,stringsAsFactors=FALSE)
colnames(genes_per_region) <- c("region_name","num_genes")
region_metadata$region_name <- rownames(region_metadata)
old_regions <- region_metadata$region_name
region_metadata <- merge(region_metadata,genes_per_region,by="region_name")
setdiff(old_regions,region_metadata$region_name)

#' 
#' Looks like chr1:122500000-124760000 was removed because it has 0 genes matched to scRNA gene annotation. That's fine.
#' 
#' Let's proceed. Order regions by chr1,2,3, then start coordinate.
#' 
## -----------------------------------------------------------------------------
info <- info[,region_metadata$region_name]
myorder <- order(region_metadata$chr,region_metadata$start)
region_metadata <- region_metadata[myorder,]
info <- info[,myorder]

#' 
#' Get the mean per chromosome by taking a weighted average based on the length of regions.
#' 
## -----------------------------------------------------------------------------
mean_per_chr <- matrix(NA,nrow=nrow(info),ncol=22,dimnames=list(rownames(info),paste0("chr",1:22)))

for(chr in 1:22)
{
	chr_indices <- chr_indices <- which(region_metadata$chr == chr)
	info_this_chr <- info[,chr_indices]
	region_metadata_this_chr <- region_metadata[chr_indices,]
	region_metadata_this_chr$prop_chr <- region_metadata_this_chr$length/sum(region_metadata_this_chr$length)
	info_this_chr <- sweep(info_this_chr,2,region_metadata_this_chr$prop_chr,FUN="*")
	mean_per_chr[,chr] <- apply(info_this_chr,1,sum)
}

#' 
#' Get the median of these values across chromosomes, for each cell.
#' 
#' This is our custom "ploidy" value.
#' 
#' Then, remove cells with custom "ploidy" value more than one SD from the mean.
#' 
#' First include these cells on heatmap so we can see what they look like, then remove and replot heatmap.
#' 
## -----------------------------------------------------------------------------
median_of_means_per_chr <- apply(mean_per_chr,1,median)

min_median_of_means_per_chr <- mean(median_of_means_per_chr) - sd(median_of_means_per_chr)
max_median_of_means_per_chr <- mean(median_of_means_per_chr) + sd(median_of_means_per_chr)

custom_ploidy_category_per_cell <- rep("In_range",times=nrow(info))
custom_ploidy_category_per_cell[median_of_means_per_chr < min_median_of_means_per_chr] <- "Low"
custom_ploidy_category_per_cell[median_of_means_per_chr > max_median_of_means_per_chr] <- "High"

custom_ploidy_category_per_cell_dat <- data.frame(ploidy = custom_ploidy_category_per_cell,
	row.names=rownames(info))

#' 
## -----------------------------------------------------------------------------
mean_per_chr[mean_per_chr > 6] <- 6 #Truncate to show range in lower values better.

pheatmap(mean_per_chr,
	annotation_row=custom_ploidy_category_per_cell_dat,
	show_rownames=FALSE,
	cluster_cols=FALSE)

#' 
## -----------------------------------------------------------------------------
info <- info[median_of_means_per_chr >= min_median_of_means_per_chr & median_of_means_per_chr <= max_median_of_means_per_chr,]
mean_per_chr <- mean_per_chr[rownames(info),]

pheatmap(mean_per_chr,
	show_rownames=FALSE,
	cluster_cols=FALSE)

#' 
#' The following chromosomes seem to warrant further review: 1,3,7,10,12,21
#' 
#' Let's plot the actual regions for these.
#' 
## -----------------------------------------------------------------------------
info_trunc <- info
info_trunc[info_trunc > 6] <- 6

region_indices <- which(region_metadata$chr %in% c(1,3,7,10,12,21))

region_metadata$chr_char <- paste0("chr",region_metadata$chr)

rownames(region_metadata) <- region_metadata$region_name

#' 
## -----------------------------------------------------------------------------
pheatmap(info_trunc[,region_indices],
    show_rownames=FALSE,cluster_cols=FALSE,show_colnames=FALSE,
    annotation=region_metadata[,c("chr_char","end","length")])

#' 
#' Let's look at each chromosome individually for some of these, so we can read the actual region coordinates.
#' 
## -----------------------------------------------------------------------------
pheatmap(info_trunc[,region_metadata$chr == 1],
    show_rownames=FALSE,cluster_cols=FALSE)

#' 
## -----------------------------------------------------------------------------
pheatmap(info_trunc[,region_metadata$chr == 3],
    show_rownames=FALSE,cluster_cols=FALSE)

#' 
## -----------------------------------------------------------------------------
pheatmap(info_trunc[,region_metadata$chr == 7],
    show_rownames=FALSE,cluster_cols=FALSE)

#' 
## -----------------------------------------------------------------------------
pheatmap(info_trunc[,region_metadata$chr == 12],
    show_rownames=FALSE,cluster_cols=FALSE)

#' 
#' In chr1, looks like the region with copy number 3 vs. 4 starts at 149,880,000. This happens to also mean that we want to focus on the q arm.
#' 
#' In chr3, the earlier coordinates are mostly consistent between cells (copy number=2). Then, let's focus only on the q arm (starts at >90Mb).
#' 
#' In chr7, later coordinates either are mostly the same between cells (copy number=3), or have lower frequency amplifications that are not consistent with the division into 2-3 clones we otherwise seem to be seeing. Let's subset to just the p arm (ends at ~60Mb).
#' 
#' Finally in chr12, the p arm (up to ~35Mb) seems to be all the same (amplified) copy number. The start of the q arm also shows little variability (all copy number=2). So focus on only the end of the q arm.
#' 
## -----------------------------------------------------------------------------
region_indices <- setdiff(region_indices,which(region_metadata$chr == 1 & region_metadata$start < 149e6))
region_indices <- setdiff(region_indices,which(region_metadata$chr == 3 & region_metadata$start < 90e6))
region_indices <- setdiff(region_indices,which(region_metadata$chr == 7 & region_metadata$start > 60e6))
region_indices <- setdiff(region_indices,which(region_metadata$chr == 12 & region_metadata$start < 80e6))

#' 
#' Remake heatmap.
#' 
## -----------------------------------------------------------------------------
pheatmap(info_trunc[,region_indices],
    show_rownames=FALSE,cluster_cols=FALSE,show_colnames=FALSE,
    annotation=region_metadata[,c("chr_char","end","length")])

#' 
#' Let's also look at chr3 again.
#' 
## -----------------------------------------------------------------------------
pheatmap(info_trunc[,region_metadata$chr == 3 & 1:nrow(region_metadata) %in% region_indices],
    show_rownames=FALSE,cluster_cols=FALSE)

#' 
#' Looks like ~93.7-118.9Mb and ~148.58Mb onward show the most consistent pattern across adjacent segments. Let's remove the regions in between.
#' 
## -----------------------------------------------------------------------------
chr3_indices_to_remove <- which(region_metadata$chr == 3 & region_metadata$start >120e6 & region_metadata$start < 148e6)
region_indices <- setdiff(region_indices,chr3_indices_to_remove)

#' 
#' Also look at chr12 again.
#' 
## -----------------------------------------------------------------------------
pheatmap(info_trunc[,region_metadata$chr == 12 & 1:nrow(region_metadata) %in% region_indices],
         show_rownames=FALSE,cluster_cols=FALSE)

#' 
#' We see amplification in the most cells if we subset to 93.5-123.6Mb.
#' 
## -----------------------------------------------------------------------------
chr12_indices_to_remove <- which(region_metadata$chr == 12 & (region_metadata$start <93e6 | region_metadata$start > 126e6))
region_indices <- setdiff(region_indices,chr12_indices_to_remove)

#' 
#' Let's get the number of genes per chromosome after all this subsetting.
#' 
## -----------------------------------------------------------------------------
region_metadata[region_indices,] %>% group_by(chr) %>% summarize(genes_per_chr = sum(num_genes))

#' 
#' And redo the heatmap yet again.
#' 
## -----------------------------------------------------------------------------
pheatmap(info_trunc[,region_indices],
         show_rownames=FALSE,cluster_cols=FALSE,
         annotation=region_metadata[,c("chr_char","end","length")],
         show_colnames=FALSE)

#' 
#' What are the numbers of cells with each number of regions amplified (defined as the higher of the two main copy numbers present) for each chromosome?
#' 
## -----------------------------------------------------------------------------
chr_of_interest <- c(1,3,7,10,12,21)
higher_copy_num <- c(4,3,4,3,3,4)

for(i in 1:length(chr_of_interest))
{
  indices_this_chr <- intersect(region_indices,which(region_metadata$chr == chr_of_interest[i]))
  regions_amp_per_cell <- apply(info[,indices_this_chr],1,FUN=function(x)sum(x >= higher_copy_num[i]))
  print(chr_of_interest[i])
  print(table(regions_amp_per_cell))
}

#' 
#' Let's define the following number of regions to say that a chromosome/chromosome arm has high copy number.
#' 
#' chr1 - 16/20
#' 
#' chr3 - 12/14
#' 
#' chr7 - 6/8
#' 
#' chr10 - 14/18
#' 
#' chr12 - 4/5
#' 
#' chr21 - 3/4
#' 
#' Then, see what percent of cells have higher/lower copy number for each.
#' 
## -----------------------------------------------------------------------------
amp_binary_matrix <- matrix(NA,nrow=nrow(info),ncol=6)
rownames(amp_binary_matrix) <- rownames(info)
colnames(amp_binary_matrix) <- paste0("chr",c(1,3,7,10,12,21))

chr_of_interest <- c(1,3,7,10,12,21)
higher_copy_num <- c(4,3,4,3,3,4)
num_regions_required <- c(16,12,6,14,4,3)

for(i in 1:length(chr_of_interest))
{
  indices_this_chr <- intersect(region_indices,which(region_metadata$chr == chr_of_interest[i]))
  regions_amp_per_cell <- apply(info[,indices_this_chr],1,FUN=function(x)sum(x >= higher_copy_num[i]))
  regions_amp_per_cell <- as.numeric(as.vector(regions_amp_per_cell))
  amp_binary_matrix[,i] <- ifelse(regions_amp_per_cell >= num_regions_required[i],1,0)
}

colSums(amp_binary_matrix)

#' 
#' We know that the following properties are needed for a segment x clone matrix that will result in the most confident clonealign calls:
#' 
#' 1. Largest possible number of genes to define each clone, since a lot of genes will not be expressed at all in sparse scRNA data.
#' 2. Amplifications (higher copy number) easier to detect than deletions, also due to sparsity. So ideally have one large set of genes amplified in clone A relative to clone B, and another amplified in clone B relative to clone A.
#' 3. Breaking up into a larger number of clones may mean a lower number of amplified genes to define each clone. In which case, a 2+ clone model may be worse than a 3+ clone model.
#' 
#' It looks from the heatmap like one set of cells has the following set of chromosomes amplified: 1,3,10. Versus the other set of cells has the following: 7,12,21.
#' 
#' So, we should include at least one chromosome from each of these groups in our final clone definitions.
#' 
#' We also find that chr1 and chr10 have relatively fewer cells with amplification. So perhaps the clone with these cells amplified will be our secondary clone.
#' 
#' Let's look at some confusion matrices for number of cells with amplification in one vs. both, for different pairs of chromosomes.
#' 
#' Also get a correlation coefficient.
#' 
## -----------------------------------------------------------------------------
table(amp_binary_matrix[,"chr1"],
      amp_binary_matrix[,"chr3"])
cor(amp_binary_matrix[,"chr1"],amp_binary_matrix[,"chr3"])

table(amp_binary_matrix[,"chr1"],
      amp_binary_matrix[,"chr10"])
cor(amp_binary_matrix[,"chr1"],amp_binary_matrix[,"chr10"])

table(amp_binary_matrix[,"chr10"],
      amp_binary_matrix[,"chr3"])
cor(amp_binary_matrix[,"chr3"],amp_binary_matrix[,"chr10"])

#' 
## -----------------------------------------------------------------------------
table(amp_binary_matrix[,"chr7"],
      amp_binary_matrix[,"chr12"])
cor(amp_binary_matrix[,"chr7"],amp_binary_matrix[,"chr12"])

table(amp_binary_matrix[,"chr7"],
      amp_binary_matrix[,"chr21"])
cor(amp_binary_matrix[,"chr7"],amp_binary_matrix[,"chr21"])

table(amp_binary_matrix[,"chr12"],
      amp_binary_matrix[,"chr21"])
cor(amp_binary_matrix[,"chr12"],amp_binary_matrix[,"chr21"])

#' 
## -----------------------------------------------------------------------------
table(amp_binary_matrix[,"chr1"],
      amp_binary_matrix[,"chr7"])
cor(amp_binary_matrix[,"chr1"],amp_binary_matrix[,"chr7"])

table(amp_binary_matrix[,"chr3"],
      amp_binary_matrix[,"chr7"])
cor(amp_binary_matrix[,"chr3"],amp_binary_matrix[,"chr7"])

table(amp_binary_matrix[,"chr10"],
      amp_binary_matrix[,"chr7"])
cor(amp_binary_matrix[,"chr10"],amp_binary_matrix[,"chr7"])

#' 
#' Some notes:
#' 
#' 1. Looking at chromosomes with cells from the smaller clone amplified (1,3,10), we find many cells amplified in chr3 but not chr1/chr10.
#' 2. If we want to use two chromosomes to define the larger clone (which we probably should since none of chr7,12,21 have a ton of genes on their own), chr 7 and 12 have the best correlation.
#' 3. Chromosomes 1 and 10 have a lot more genes (1122 and 960 respectively) compared to chr3 (502). If we want to choose one of each of chr1/3/10 and chr7/12/21 as a first pass for defining clones, we should thus probably start with chr1 or chr10 over chr3.
#' 4. Finally, chr10 has a very slightly better correlation rate with chr7 than chr1 does, so that more cells may be included in one of the two clones (1025+521=1546 compare to 1032+451=1483). So let's start with using chr7 and chr10 to define clones.
#' 
#' Define clones, then remake heatmap with this info.
#' 
## -----------------------------------------------------------------------------
clones <- rep("Undefined",times=nrow(info))
clones[amp_binary_matrix[,"chr7"] == 1 & amp_binary_matrix[,"chr10"] == 0] <- "clone1"
clones[amp_binary_matrix[,"chr7"] == 0 & amp_binary_matrix[,"chr10"] == 1] <- "clone2"

clones_dat <- data.frame(clone = clones,row.names=rownames(info))

pheatmap(info_trunc[,region_indices],
         show_rownames=FALSE,cluster_cols=FALSE,show_colnames=FALSE,
         annotation_row=clones_dat,
         annotation=region_metadata[,c("chr_char","end","length")])

#' 
#' And get the confusion matrix with remaining chromosomes.
#' 
## -----------------------------------------------------------------------------
table(clones,amp_binary_matrix[,"chr1"])
table(clones,amp_binary_matrix[,"chr3"])
table(clones,amp_binary_matrix[,"chr12"])
table(clones,amp_binary_matrix[,"chr21"])

#' 
#' At first glance, none of these accuracy rates are that bad.
#' 
#' However, consider the lowest of the two rates if we separate clone1 and clone2.
#' 
#' chr1 - 335/(335 + 186)=64%
#' 
#' chr3 - 862/(862 + 163)=84%
#' 
#' chr12 - 880/(880 + 145)=86%
#' 
#' chr21 - 759/(759 + 266)=74%
#' 
#' A fairly large percent of clone2 cells (as we have defined them) do not have chr1 amplification.
#' 
#' While a large percent of clone1 cells do not have chr21 amplification.
#' 
#' In theory, it seems like we could have chr3 included as another chromosome amplified in clone2, though for whatever reason I did not for the publication. For this markdown, I will also print the coordinates for that chromosome, and will run clonealign including it. I expect to see that including chr3 increases the number of cells with sufficient expression for clonealign, but should change very few calls for the cells included in the paper (figures 2e-f).
#' 
#' Make another heatmap including only chr 7,10, and 12.
#' 
## -----------------------------------------------------------------------------
region_indices_flt <- region_indices[region_metadata$chr[region_indices] %in% c(7,10,12)]

pheatmap(info_trunc[,region_indices_flt],
         show_rownames=FALSE,cluster_cols=FALSE,fontsize_col=8,
         annotation_row=clones_dat,
         annotation=region_metadata[,c("chr_char","end","length")])

#' 
#' Fewer cells have amplification for first region of chr7. Just remove that and then will have the set of coordinates used for clonealign.
#' 
## -----------------------------------------------------------------------------
region_indices_flt <- region_indices_flt[2:length(region_indices_flt)]

#' 
#' Let's print a version also with chr3 included, just to see how that looks in clonealign.
#' 
## -----------------------------------------------------------------------------
region_indices_flt <- c(region_indices_flt,region_indices[region_metadata$chr[region_indices] == 3])

to_print <- region_metadata[region_indices_flt,]
to_print <- data.frame(chr = paste0("chr",to_print$chr),
                       to_print[,c("start","end")],
                       stringsAsFactors=FALSE)
clone1_copy_num <- rep(NA,times=nrow(to_print))
clone2_copy_num <- rep(NA,times=nrow(to_print))
clone1_copy_num[to_print$chr == "chr3"] <- 2
clone2_copy_num[to_print$chr == "chr3"] <- 3 #Really 2 vs. 4 for part of the chromosome in some cells, but others have only 3, so we have to set it to the lowest value since this is a two-clone model.
clone1_copy_num[to_print$chr == "chr10"] <- 2
clone2_copy_num[to_print$chr == "chr10"] <- 3
clone1_copy_num[to_print$chr == "chr7"] <- 4
clone2_copy_num[to_print$chr == "chr7"] <- 3
clone1_copy_num[to_print$chr == "chr12"] <- 3
clone2_copy_num[to_print$chr == "chr12"] <- 2

to_print$name <- paste0(clone1_copy_num,"_",clone2_copy_num)

write.table(to_print,file="MP2B_segment_by_clone_incl_chr3.bed",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")

to_print <- to_print[to_print$chr != "chr3",]

write.table(to_print,file="MP2B_segment_by_clone.bed",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")

