library(pheatmap)
.libPaths("/nfs/sw/clonealign/clonealign-2.0/R/lib64/R/library")
suppressMessages(library(Matrix))
suppressMessages(library(tensorflow))
suppressMessages(library(clonealign))
sample="BOSTON" #This was the path for MP2B on NYGC system.
RNA_path=paste0("../download/NYGC_collab/3prime_matrices/",sample,"_filtered_feature_bc_matrix/")
RNA_counts <- readMM(paste0(RNA_path,"matrix.mtx.gz"))
RNA_features <- read.table(paste0(RNA_path,"features.tsv.gz"),header=FALSE,stringsAsFactors=FALSE)
rownames(RNA_counts) <- make.unique(RNA_features[,2],sep="-")
colnames(RNA_counts) <- readLines(paste0(RNA_path,"barcodes.tsv.gz"))
#clone_copy_num <- read.csv("MP2B_gene_by_clone.csv",header=TRUE,row.names=1,check.names=FALSE)
clone_copy_num <- read.csv("MP2B_gene_by_clone_incl_chr3.csv",header=TRUE,row.names=1,check.names=FALSE)
RNA_counts <- RNA_counts[intersect(rownames(RNA_counts),rownames(clone_copy_num)),]
clone_copy_num <- clone_copy_num[rownames(RNA_counts),]
RNA_counts <- Matrix(as.matrix(RNA_counts),sparse=TRUE)
ca_data <- preprocess_for_clonealign(gene_expression_data = t(RNA_counts),
    copy_number_data = clone_copy_num)
RNA_counts_flt <- t(ca_data$gene_expression_data)
copy_num_flt <- ca_data$copy_number_data
set.seed(1392)
model_from_clonealign <- run_clonealign(ca_data$gene_expression_data,
    ca_data$copy_number_data,
    clone_call_probability=0.90)
to_print <- data.frame(cell = colnames(RNA_counts_flt),
    clone = model_from_clonealign$clone)
#write.csv(to_print,file="Results_clonealign/MP2B_clonealign_calls_05.05.22.csv",row.names=FALSE,quote=FALSE)
write.csv(to_print,file="Results_clonealign/MP2B_clonealign_calls_incl_chr3_05.05.22.csv",row.names=FALSE,quote=FALSE)
