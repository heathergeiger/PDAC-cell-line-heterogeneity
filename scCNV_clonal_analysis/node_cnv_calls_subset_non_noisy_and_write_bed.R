args <- commandArgs(trailingOnly=TRUE)
sample <- args[1]
bed <- args[2]
summary <- args[3]

bed <- read.table(bed,header=TRUE,skip=2,comment.char="")
colnames(bed)[1] <- "chrom"
bed$id <- paste0("cell",bed$id)

summary <- read.csv(summary,row.names=2,check.names=FALSE)
rownames(summary) <- paste0("cell",rownames(summary))

non_noisy_cells <- rownames(summary)[summary$is_noisy == 0 & summary$mean_ploidy >= 1.5]

bed <- bed[bed$id %in% non_noisy_cells,]

bed <- data.frame(bed[,c("chrom","start","end")],
    name=paste0(bed$chrom,":",bed$start,"-",bed$end,"_",bed$id,"_",bed$copy_number,"_",bed$event_confidence))

write.table(bed,file=paste0(sample,"_non_noisy_cell_node_cnv_calls.bed"),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
