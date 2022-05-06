#' Using MP2B as an example in this script.
#' 
#' First, process the scRNA data the usual way.
#' 
#' For now, include all cells, regardless of whether or not they have sufficient expression to be called by clonealign (with or without chr3).
#' 
## ----eval=FALSE---------------------------------------------------------------
## library(Seurat)

#' 
## ----eval=TRUE----------------------------------------------------------------
library(ggplot2)

#' 
## ----eval=FALSE---------------------------------------------------------------
## RNA_path="../download/NYGC_collab/3prime_matrices/BOSTON_filtered_feature_bc_matrix/"
## raw_counts <- Read10X(RNA_path)
## seurat.obj <- CreateSeuratObject(raw_counts,min.cells=10,min.features=200)
## seurat.obj[["percent.mito"]] <- PercentageFeatureSet(seurat.obj,pattern="^MT-")
## seurat.obj <- subset(seurat.obj,percent.mito < 10) #Very few cells with higher mito than this, looking at distribution.
## seurat.obj <- SCTransform(seurat.obj,conserve.memory=TRUE)
## seurat.obj <- RunPCA(seurat.obj)
## seurat.obj <- FindNeighbors(seurat.obj,reduction = "pca",dims = 1:30)
## seurat.obj <- RunUMAP(seurat.obj,reduction = "pca",dims = 1:30)
## 
## umap_coords <- Embeddings(seurat.obj,reduction="umap")
## write.csv(umap_coords,file="MP2B_only_umap_coords_very_lightly_filtered_RNA.csv",row.names=TRUE,quote=FALSE)

#' 
#' Read in clonealign output, both with and without chr3.
#' 
## -----------------------------------------------------------------------------
excl_chr3 <- read.csv("MP2B_clonealign_calls_05.05.22.csv",header=TRUE,stringsAsFactors=FALSE)
incl_chr3 <- read.csv("MP2B_clonealign_calls_incl_chr3_05.05.22.csv",header=TRUE,,stringsAsFactors=FALSE)

incl_chr3_only_cells <- setdiff(incl_chr3$cell,excl_chr3$cell)
to_rbind <- data.frame(cell = incl_chr3_only_cells,clone = "no_call")
excl_chr3 <- rbind(excl_chr3,to_rbind)

rownames(excl_chr3) <- excl_chr3$cell
excl_chr3 <- excl_chr3[incl_chr3$cell,]

clone_info <- data.frame(excl_chr3 = excl_chr3$clone,
                         incl_chr3 = incl_chr3$clone,
                         row.names=excl_chr3$cell,
                         check.names=FALSE)

table(clone_info$excl_chr3,clone_info$incl_chr3)

#' 
#' We find that very few cells change their clone assignment whether we include or exclude chr3, but we get a lot more cells able to be assigned if we include it.
#' 
#' Let's therefore use the calls incl. chr3 from here on out, and plot vs. the scRNA UMAP including all cells after light filtering (not just those from clonealign output).
#' 
## -----------------------------------------------------------------------------
umap_coords <- read.csv("MP2B_only_umap_coords_very_lightly_filtered_RNA.csv",header=TRUE,row.names=1,check.names=FALSE)

clonealign_missing_cells <- setdiff(rownames(umap_coords),rownames(clone_info))

to_rbind <- data.frame(excl_chr3 = rep("no_call",times=length(clonealign_missing_cells)),
                       incl_chr3 = "no_call",
                       row.names=clonealign_missing_cells)

clone_info <- rbind(clone_info,to_rbind)
clone_info <- clone_info[rownames(umap_coords),]

umap_coords$clone <- clone_info$incl_chr3

ggplot(umap_coords,aes(UMAP_1,UMAP_2,colour=clone)) + geom_point() + scale_colour_manual(values=c("#009292","#FF6DB6","#000000","#999999"))

#' 
#' We find that clone1 and clone2 show strong separation even when we include all cells in calculating UMAP (not just those output by clonealign).
#' 
#' It is likely that the gray cells in the green cluster show similar expression to Clone1 cells, while the gray cells in the pink cluster show similar expression to Clone2 cells.
