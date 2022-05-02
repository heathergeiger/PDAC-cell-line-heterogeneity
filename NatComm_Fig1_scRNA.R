library(dplyr)
library(reshape2)
library(stringr)
library(umap)
library(Seurat)
## read in objects
readInData <- function(data, projname, add.ids = "cols"){
  if(add.ids == "rows"){
    temp <- melt(data, id.vars = "CellId") %>% dcast(variable ~ CellId)
    temp$variable <-  sapply(temp$variable, function(x) str_replace(x, "[.]", "-"))
    rownames(temp) <- temp$variable
    temp$variable <- NULL
    colnames(temp) <- paste(projname, colnames(temp), sep="_")
    data <- temp
  } else {
    colnames(data) <- paste(projname, colnames(data), sep="_") 
  }
  return(data %>% CreateSeuratObject(project=projname, min.cells = 3, min.features = 200)
         %>% NormalizeData(verbose=FALSE)
         %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
  )
  data <- as.data.frame(data)
  data$Gene <- rownames(data)
  
  return(data)
}
Panc1.3D <- Read10X("~/Box/Single_Cell_10x_Data/Cell_Lines/3prime/panc2d3d_rna_122019/panc3d_redone/filtered_feature_bc_matrix/") %>% readInData(projname = "Panc1.3D")
Panc1.2D <- Read10X("~/Box/Single_Cell_10x_Data/Cell_Lines/3prime/panc2d3d_rna_122019/panc2d_redone/filtered_feature_bc_matrix") %>% readInData(projname = "Panc1.2D")
MP2.A <- Read10X("~/Box/Single_Cell_10x_Data/Cell_Lines/3prime/3prime_matrices/ATCC_filtered_feature_bc_matrix/") %>% readInData(projname = "MP2.A")
BXPC3 <- Read10X("~/Box/Single_Cell_10x_Data/Cell_Lines/3prime/3prime_matrices/BXPC3_filtered_feature_bc_matrix/") %>% readInData(projname = "BXPC3")
MP2.C <- Read10X("~/Box/Single_Cell_10x_Data/Cell_Lines/3prime/3prime_matrices/LAB_filtered_feature_bc_matrix/") %>% readInData(projname = "MP2.C")
MP2.B <- Read10X("~/Box/Single_Cell_10x_Data/Cell_Lines/3prime/3prime_matrices/BOSTON_filtered_feature_bc_matrix/") %>% readInData(projname = "MP2.B")
HPDE <- Read10X("~/Box/Single_Cell_10x_Data/Cell_Lines/3prime/3prime_matrices/HPDE_filtered_feature_bc_matrix/") %>% readInData(projname = "HPDE")
HPNE <- Read10X("~/Box/Single_Cell_10x_Data/Cell_Lines/3prime/3prime_matrices/HPNE_filtered_feature_bc_matrix/") %>% readInData(projname = "HPNE")
HPAF2 <- Read10X("~/Box/Single_Cell_10x_Data/Cell_Lines/3prime/3prime_matrices/HPAF2_filtered_feature_bc_matrix/") %>% readInData(projname = "HPAF2")

##begin analysis here
all_combined <- merge(Panc1.2D, y = c(Panc1.3D, MP2.A, MP2.B, MP2.C, BXPC3, HPAF2, HPDE, HPNE), project = "all_merged")
options(future.globals.maxSize = 8000 * 1024^2)
all_combined[["percent.mt"]] <- PercentageFeatureSet(all_combined, pattern = "^MT-")
all_combined <- subset(all_combined, subset = percent.mt < 20)
FeaturePlot(all_combined, features = "percent.mt", pt.size = 1, split.by = "orig.ident", reduction = "umap")
all_combined <- NormalizeData(object = all_combined, normalization.method = "LogNormalize", scale.factor = 10000)

all_combined <- FindVariableFeatures(object = all_combined, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(x = all_combined)

all_combined <- ScaleData(object = all_combined, verbose = FALSE)
all_combined <- RunPCA(object = all_combined, npcs = 30, verbose = FALSE)
#newpanc <- write.csv(Embeddings(object = all_combined, reduction = "pca"), file = "PCA_embeddings_NO_oldPanc2D.csv")


# UMAP and Clustering
#all_combined <- RunTSNE(object = all_combined, reduction ="pca", dims=1:10)
all_combined <- RunUMAP(object = all_combined, reduction = "pca", dims = 1:10)
all_combined <- FindNeighbors(object = all_combined, reduction = "pca", dims = 1:10)
all_combined <- FindClusters(all_combined, resolution = 0.5)

lines_p1 <- DimPlot(object = all_combined, reduction = "umap", group.by = "orig.ident")
lines_sampsplit <- DimPlot(object = all_combined, reduction = "umap", split.by = "orig.ident")
lines_p2 <- DimPlot(object = all_combined, reduction = "umap", label = TRUE)

lines_fp <- FeaturePlot(all_combined, c("EPCAM", "VIM", "KRT8","CD44", "NES", "NOTCH2" ), min.cutoff = "q9", reduction="umap", ncol = 3)
deer_fp <- FeaturePlot(all_combined, c("COL1A1", "COL4A1", "FN1", "LAMC1", "LAMB2", "VEGFA", "KRAS", "SMAD4", "PTGS2", "TP53", "CDKN2A", "CXCL8"), min.cutoff = "q9", reduction = "umap")
moffat_fp <- FeaturePlot(all_combined, c("KRAS", "CDKN2A", "CDK6", "MYST3", "KAT6A", "SSX4", "PRKAA1", "TRAF6", "TLR4", "SKAP1", "EIF3H"), min.cutoff = "q9", reduction = "umap")

vahid_vln <- VlnPlot(all_combined, c("CD151", "CLDN4", "EPCAM", "CD9", "CD63", "CD81"), group.by = "orig.ident")


lines_fp_violin <- VlnPlot(all_combined, c("EPCAM", "VIM", "KRT8", "MUC1", "CD44"), group.by = 'orig.ident')
lines_epcam_qc <- AverageExpression(all_combined, assays = "RNA", add.ident = "orig.ident", features = "EPCAM")
write.csv(lines_epcam_qc, "epcam_gexp_per_sample.csv")

lines_vim_qc <- AverageExpression(all_combined, assays = "RNA", add.ident = "orig.ident", features = "VIM")
write.csv(lines_vim_qc, "vimentin_gexp_per_sample.csv")


lines_krt8_qc <- AverageExpression(all_combined, assays = "RNA", add.ident = "orig.ident", features = "KRT8")
write.csv(lines_krt8_qc, "krt8_gexp_per_sample.csv")


lines_cd44_qc <- AverageExpression(all_combined, assays = "RNA", add.ident = "orig.ident", features = "CD44")
write.csv(lines_cd44_qc, "cd44_gexp_per_sample.csv")


######
quant_test1 <- table(all_combined@meta.data$orig.ident)

test_plot <- ggplot(all_combined@meta.data, aes(x = orig.ident, fill = seurat_clusters)) + 
  geom_bar(position = "fill", width = .7, color = "black", lwd = .2) + 
  theme(panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.line=element_line(color = "black"),
        axis.ticks=element_line(color = "black"),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.key = element_blank(),
        aspect.ratio = 1.5) + 
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(drop = F) +
  xlab("Sample") +
  ylab("Proportion of sample")


###########################Remaking BRCA1 and AURKA vln plots for HPDE HPNE FIG1
normals <- merge(HPNE, y = HPDE, project = "normals")
normals[["percent.mt"]] <- PercentageFeatureSet(normals, pattern = "^MT-")
normals <- subset(normals, subset = percent.mt < 20)
FeaturePlot(normals, features = "percent.mt", pt.size = 1, split.by = "orig.ident", reduction = "umap")
normals <- NormalizeData(object = normals, normalization.method = "LogNormalize", scale.factor = 10000)

normals <- FindVariableFeatures(object = normals, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(x = normals)

normals <- ScaleData(object = normals, verbose = FALSE)
normals <- RunPCA(object = normals, npcs = 30, verbose = FALSE)

VlnPlot(normals, features = c("BRCA1", "AURKA"))

#######generating the QC plots for reviewer 5
# Visualize QC metrics as a violin plot
VlnPlot(all_combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(all_combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all_combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
