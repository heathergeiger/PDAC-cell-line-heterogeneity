#new panc2d/3d scRNA data, as per experimental redo 12/2019


library(Seurat)
library(dplyr)
library(reshape2)
library(stringr)
## read in object
# FNAs <- readRDS("~/Downloads/PDAC_tutorial.rds")
# 
# DimPlot(object=FNAs, reduction = "tsne")
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
#panc3d <- Read10X("/Users/maitralab/Desktop/panc1_redone/panc2d3d_rna_122019/panc3d_redone/filtered_feature_bc_matrix") %>% readInData(projname = "Panc1.3D")
#panc2d <- Read10X("/Users/maitralab/Desktop/panc1_redone/panc2d3d_rna_122019/panc2d_redone/filtered_feature_bc_matrix") %>% readInData(projname = "Panc1.2D")

panc2d <- Read10X("~/Box/Single_Cell_10x_Data/Cell_Lines/3prime/panc2d3d_rna_122019/panc2d_redone/filtered_feature_bc_matrix/") %>% readInData(projname = "Panc1.2D")
panc3d <- Read10X("~/Box/Single_Cell_10x_Data/Cell_Lines/3prime/panc2d3d_rna_122019/panc3d_redone/filtered_feature_bc_matrix/") %>% readInData(projname = "Panc1.3D")


all_mp2 <- merge(panc2d, y = panc3d, project = "panc1")

all_mp2[["percent.mt"]] <- PercentageFeatureSet(all_mp2, pattern = "^MT-")
all_mp2 <- subset(all_mp2, subset = percent.mt < 20)
#FeaturePlot(all_mp2, features = "percent.mt", pt.size = 1, split.by = "orig.ident", reduction = "umap")

mito.genes <- grep(pattern = "^MT-", x = rownames(x = all_mp2@meta.data), value = TRUE)

# AddMetaData adds columns to object@meta.data, and is a great place to
# stash QC stats

all_mp2 <- NormalizeData(object = all_mp2, normalization.method = "LogNormalize", scale.factor = 10000)
all_mp2 <- FindVariableFeatures(object = all_mp2, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(x = all_mp2)
all.processed <- ScaleData(object = all_mp2, verbose = FALSE)
all.processed <- RunPCA(object = all.processed, npcs = 30, verbose = FALSE)

#View PCA outputs
print(all.processed[["pca"]], dims = 1:5, nfeatures = 10)
VizDimLoadings(all.processed, dims = 1:2, reduction = "pca")

# t-SNE and Clustering
all.processed <- RunTSNE(object = all.processed, reduction ="pca", dims=1:10)
all.processed <- RunUMAP(object = all.processed, reduction = "pca", dims = 1:10)
all.processed <- FindNeighbors(object = all.processed, reduction = "pca", dims = 1:10)
all.processed <- FindClusters(all.processed, resolution = 0.4)

p1 <- DimPlot(object = all.processed, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(object = all.processed, reduction = "umap", label = TRUE)
sampsplit <- DimPlot(object = all.processed, reduction = "umap", split.by = "orig.ident", label = TRUE)
plot_emt <- FeaturePlot(all.processed, c("EPCAM", "VIM", "KRT8","GATA6"), min.cutoff = "q9", reduction="umap")


Idents(all.processed) <- all.processed$orig.ident
Idents(all.processed) <- all.processed$seurat_clusters

panc2.markers <- FindMarkers(object = all.processed, ident.1 = c(3, 5, 6), min.pct = 0.25)
write.csv(panc2.markers,  file = "panc2_markers.csv")
panc3.markers <- FindMarkers(object = all.processed,  ident.1 = c(0, 1, 2, 4), min.pct = 0.25)
write.csv(panc3.markers,  file = "panc3_markers.csv")

#looking for relevant expression in genes that correlate to snATACseq motif networks. 
#DotPlot(pbmc, features = features) + RotatedAxis()

sp_motif <- DotPlot(all.processed, features = c('SOD1', 'PHGDH', 'ID1', 'CYP1B1', 'NDUFV1', 'NDUFV2'))
jun_motif <- DotPlot(all.processed, features = c('TGFB1', 'DDIT3', 'SMAD3', 'SMAD4', 'RCAN1') + RotatedAxis())


#looking for panc2d-specific markers here. Did not run on 3/18/2021. 
#panc2_c3.markers <- FindMarkers(object = all.processed, ident.1 = 3, ident.2 = c(0,1,2,4,5,6))
#panc2_c5.markers <- FindMarkers(object = all.processed, ident.1 = 5, ident.2 = c(0,1,2,3,4,6))
#panc2_c6.markers <- FindMarkers(object = all.processed, ident.1 = 6, ident.2 = c(0,1,2,3,4,5))

#looking for panc3d-specific markers here (comparisons to 2d clusters 2,6,7). Also did not run on 3/18/2021. 
#panc3_c0.markers <- FindMarkers(object = all.processed, ident.1 = 0, ident.2 = c(1,2,3,4,5,6))
#panc3_c2.markers <- FindMarkers(object = all.processed, ident.1 = 2, ident.2 = c(0,1,3,4,5,6))
#panc3_c1.markers <- FindMarkers(object = all.processed, ident.1 = 1, ident.2 = c(0,2,3,4,5,6))
#panc3_c12.markers <- FindMarkers(object = all.processed, ident.1 = c(1,2), ident.2 = c(0,3,4,5,6))
#panc3_c4.markers <- FindMarkers(object = all.processed, ident.1 = 4, ident.2 = c(0,1,2,3,5,6))

#in clue, clusters 2 and 1 have very few markers, i will combine them here. Did not run on 3/18/2021. 
#panc3_cx.markers <- FindMarkers(object = all.processed, ident.1 = c(1, 2), ident.2 = c(0,3,4,5,6))

#a 3d collection vs a bulky 2d collection
#panc1_subset.markers <- FindMarkers(object = all.processed, ident.1 = c(1,3,0), ident.2 = c(2,6,7))

#panc1_subset.markers <- merge(panc3.markers, panc2.markers, by.y = "orig.ident")



#fat csv writing going here 
#write.csv(panc2_c3.markers, file = "panc2_c3.csv")
#write.csv(panc2_c5.markers, file = "panc2_c5.csv")
#write.csv(panc2_c6.markers, file = "panc2_c6.csv")
#write.csv(panc3_c0.markers, file = "panc3_c0.csv")
#write.csv(panc3_c1.markers, file = "panc3_c1.csv")
#write.csv(panc3_c2.markers, file = "panc3_c2.csv")
#write.csv(panc3_c4.markers, file = "panc3_c4.csv")

#read in all marker files as lists to be used in GSEA. Example with hallmark below: 
hallmark <- gmtPathways("~/Downloads/h.all.v7.0.symbols.gmt")
go <- gmtPathways("~/Downloads/c5.all.v7.0.symbols.gmt")
kegg <- gmtPathways("~/Downloads/c2.cp.kegg.v7.0.symbols.gmt")
c2 <- gmtPathways("~/Downloads/c2.cp.v7.0.symbols.gmt")


#to use fgsea on both sets
library(fgsea)
library(ggplot2)

#for cluster panc3d
#panc3.markers <- write.csv(panc3.markers, file = "panc3d_markers.csv")

test4 = panc3.markers
test4 <- as.data.frame(test4)
test4$y=(-log10(test4$p_val_adj))
test4$k=is.infinite(test4$y)

d4 = test4
d4$y[d4$k]=runif(sum(d4$k), max(d4$y[!d4$k]), max(d4$y[!d4$k])*1.25)

#next is changing the y value to negative if the -log value is, in fact, negative
d4$y = ifelse(d4$avg_logFC<0, -d4$y, d4$y)

#order the d4 by the y
d4 = d4[order(d4$y), ]

#convert the d4 to a matrix
mat4 = as.matrix(d4)

stats4 = mat4[,6]    #must read "Named number" when done

c4_fgseaRes <- fgsea(pathways = hallmark, 
                     stats = stats4,
                     minSize=5,
                     maxSize=Inf,
                     nperm=10000)


c4_fgseaRes_kegg <- fgsea(pathways = kegg, 
                         stats = stats4, minSize=5, maxSize=Inf, nperm=10000)

c4_fgseaRes_go <- fgsea(pathways = go, 
                        stats = stats4,
                        minSize=5,
                        maxSize=Inf,
                        nperm=10000)

c4_fgseaRes_all <- fgsea(pathways = c2, 
                        stats = stats4,
                        minSize=5,
                        maxSize=Inf,
                        nperm=10000)

x <- as.data.frame(c4_fgseaRes[1:44, 1:7])
write.csv(x, file = "panc3d_hallmark.csv")

#for panc2d markers
testc5 = panc2.markers
testc5 <- as.data.frame(testc5)
testc5$y=(-log10(testc5$p_val_adj))
testc5$k=is.infinite(testc5$y)

d5 = testc5
d5$y[d5$k]=runif(sum(d5$k), max(d5$y[!d5$k]), max(d5$y[!d5$k])*1.25)

#next is changing the y value to negative if the -log value is, in fact, negative
d5$y = ifelse(d5$avg_logFC<0, -d5$y, d5$y)

#order the d5 by the y
d5 = d5[order(d5$y), ]

#convert the d5 to a matrix
mat5 = as.matrix(d5)

stats5 = mat5[,6]    #must read "Named number" when done

c5_fgseaRes <- fgsea(pathways = hallmark, 
                     stats = stats5,
                     minSize=5,
                     maxSize=Inf,
                     nperm=10000)


c5_fgseaRes_kegg <- fgsea(pathways = kegg, stats = stats5, minSize=5, maxSize=Inf, nperm=10000)

c5_fgseaRes_go <- fgsea(pathways = go, 
                        stats = stats5,
                        minSize=5,
                        maxSize=Inf,
                        nperm=10000)


c5_fgseaRes_all <- fgsea(pathways = go, 
                        stats = stats5,
                        minSize=5,
                        maxSize=Inf,
                        nperm=10000)

y <- as.data.frame(c5_fgseaRes[1:44, 1:7])
write.csv(y, file = "panc2d_hallmark.csv")

############### IMPLEMENTING DOTPLOT SCRIPT HERE################

#made a merged file of all the samples labeled by orig.ident and pathways. Thus was done by converting the fgsea results into csv files and manipulating them in Excel. 

ggplot(panc1_gsea_combined, aes(orig.ident, pathway)) + # change replacement to pathway prn
  geom_point(aes(size = -log10(padj), col = NES)) +
  scale_color_gradientn(colours = c("blue2","blue","white","red","red2"), na.value = "grey80", limits = c(-1,1)*max(abs(panc1_gsea_combined$NES))) +
  scale_size(name = "FDR", 
             breaks = c(-log10(0.005), -log10(0.05), -log10(0.5)), 
             labels = c(0.005, 0.05, 0.5)) +
  theme(panel.border=element_blank(),
        panel.background=element_rect(fill="white", color="black", size=0.5, linetype="solid"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        # plot.title=element_text(size=rel(1.25), hjust=0.5),
        axis.line=element_blank(),
        axis.text=element_text(color="black", size=rel(.75)),
        axis.text.x=element_text(angle=90, vjust=.5),
        axis.title.x=element_text(color="black", size=rel(1), angle=0),
        axis.title.y=element_text(color="black", size=rel(1), angle=90),
        axis.ticks.length=unit(0.1, "cm"),
        axis.ticks=element_line(color="black"),
        legend.key = element_blank(),
        aspect.ratio = 4) + # ratio y / x
  labs(x="Origin", y="Pathway")



# dotplot for gsea DID NOT RUN. 
# upload df named "pathways" that contains pathways to be plotted
```{r echo = F}

library(ggplot2)
df <- rbind(c4_fgseaRes, c5_fgseaRes) # merge all desired fgsea results


df <- df[df$pathway %in% pathways$pathways, ]

df$pathway <- factor(df$pathway, levels = c("Cell cycle", "Cell cycle checkpoints", "Regulation of cell cycle", "DNA repair", "DNA replication", "DNA synthesis", "Response to DNA damage", "Response to stress", "Proteasome", "Protein catabolism", "Amine metabolism", "Amino acid metabolism", "Carbohydrate metabolism", "Cellular respiration", "Lipid metabolism", "Monosaccharide metabolism", "OXPHOS", "Microtubule", "Microtubule-based process", "Microtubule organization", "mRNA processing", "RNA splicing", "Transcription"))

# alternatively
pathway_replacement = read.csv("C:/Users/jjlee3/Box/Data/NGS/scPDAC/Organoids/Analysis_new/pathway_replacement.csv")

toPlot <- fgseaResObject[fgseaResObject$pathway %in% pathway_replacement$pathway, ]
df <- merge(pathway_replacement, toPlot, by = "pathway")
df[is.na(df)] <- 0
df$replacement <- factor(df$replacement, levels = c("Cell cycle", "G2M checkpoint", "Cell proliferation", "DNA replication", "DNA repair", "Response to DNA damage", "Response to stress", "Mitotic spindle", "Microtubule-based process", "Apoptosis", "Cell junction", "Cell-substrate junction", "Glycolysis", "TCA cycle and respiration", "Electron transport chain", "OXPHOS", "ROS", "Hypoxia", "mTOR signaling", "Ser-Thr kinase activity", "Protein kinase activity", "Transcription", "Translation", "Ubiquitination", "Proteasome", "TNFa", "EMT", "Interferon signaling", "KRAS signaling", "MYC", "p53 pathway"))
#######################


ggplot(df, aes(Cluster, replacement)) + # change replacement to pathway prn
  geom_point(aes(size = -log10(padj), col = NES)) +
  scale_color_gradientn(colours = c("blue2","blue","white","red","red2"), na.value = "grey80", limits = c(-1,1)*max(abs(df$NES))) +
  scale_size(name = "FDR", 
             breaks = c(-log10(0.005), -log10(0.05), -log10(0.5)), 
             labels = c(0.005, 0.05, 0.5)) +
  theme(panel.border=element_blank(),
        panel.background=element_rect(fill="white", color="black", size=0.5, linetype="solid"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        # plot.title=element_text(size=rel(1.25), hjust=0.5),
        axis.line=element_blank(),
        axis.text=element_text(color="black", size=rel(.75)),
        axis.text.x=element_text(angle=90, vjust=.5),
        axis.title.x=element_text(color="black", size=rel(1), angle=0),
        axis.title.y=element_text(color="black", size=rel(1), angle=90),
        axis.ticks.length=unit(0.1, "cm"),
        axis.ticks=element_line(color="black"),
        legend.key = element_blank(),
        aspect.ratio = 4) + # ratio y / x
  labs(x="Cluster", y="Pathway")

ggplot(gsea_mp2_combined, aes(orig.ident, pathway)) + # change replacement to pathway prn
  geom_point(aes(size = -log10(padj), col = NES)) +
  scale_color_gradientn(colours = c("blue2","blue","white","red","red2"), na.value = "grey80", limits = c(-1,1)*max(abs(gsea_mp2_combined$NES))) +
  scale_size(name = "FDR", 
             breaks = c(-log10(0.005), -log10(0.05), -log10(0.5)), 
             labels = c(0.005, 0.05, 0.5)) +
  theme(panel.border=element_blank(),
        panel.background=element_rect(fill="white", color="black", size=0.5, linetype="solid"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        # plot.title=element_text(size=rel(1.25), hjust=0.5),
        axis.line=element_blank(),
        axis.text=element_text(color="black", size=rel(.75)),
        axis.text.x=element_text(angle=90, vjust=.5),
        axis.title.x=element_text(color="black", size=rel(1), angle=0),
        axis.title.y=element_text(color="black", size=rel(1), angle=90),
        axis.ticks.length=unit(0.1, "cm"),
        axis.ticks=element_line(color="black"),
        legend.key = element_blank(),
        aspect.ratio = 4) + # ratio y / x
  labs(x="Origin", y="Pathway")











###############

#Integrating the datasets here to see if differential pathways appear 

panc_anchors <- FindIntegrationAnchors(object.list = c(panc2d, panc3d), dims = 1:30)

panc_integrated <- IntegrateData(anchorset = panc_anchors, dims = 1:30)

DefaultAssay(panc_integrated) <- "integrated"


panc_integrated <- ScaleData(panc_integrated, verbose = FALSE)
panc_integrated <- RunPCA(panc_integrated, npcs = 30, verbose = FALSE)

panc_integrated <- RunUMAP(panc_integrated, reduction = "pca", dims = 1:20)
panc_integrated <- FindNeighbors(panc_integrated, reduction = "pca", dims = 1:20)
panc_integrated <- FindClusters(panc_integrated, resolution = 0.5)

# Visualization
panc_int1 <- DimPlot(panc_integrated, reduction = "umap", group.by = "orig.ident")
panc_int2 <- DimPlot(panc_integrated, reduction = "umap", label = TRUE)
panc_int3 <- DimPlot(panc_integrated, reduction = "umap", label = TRUE, group.by = c("seurat_clusters", "orig.ident"))



#########
#here we will try inverCNV to see if any copy number changes occurred between 3d and 2d 
library(infercnv)
library(Seurat)

#data = Read10X(data.dir = "10x_data_dir/") 
#seurat_obj = CreateSeuratObject(raw.data=data, min.cells=3, min.genes=200)
#counts_matrix = as.matrix(seurat_obj@raw.data[,seurat_obj@cell.names])

#use the pre-processed all_mp2 object here

panc_raw <- as.matrix(GetAssayData(all_mp2[["RNA"]], slot = "counts"))

write.table(panc_raw, "panc_combined_raw_matrix.txt", sep = "\t", quote = FALSE)

#this will parse cells by 2d or 3d groups
panc_annotation <- FetchData(all_mp2, vars = "orig.ident")

#confirm colnames are by the same name as fetchassay here 
head(panc_annotation)
#also gene annotation file is needed 
gene_order_file=("/Users/maitralab/Desktop/gencode_v19_gene_pos.txt")

#migrate cell name into a column 
panc_annotation$cell <- row.names(panc_annotation)
#switch so cell name is on left (switches order of first and second columns)
panc_annotation <- panc_annotation[,c(2,1)]

#eliminate row and column names
row.names(panc_annotation) <- NULL
colnames(panc_annotation) <- NULL
write.table(panc_annotation, "input_cnv_annotation.txt", row.names = F, sep = "\t", quote = F)


#change a ref group name to either null or use a comparison file

infercnv_obj = CreateInfercnvObject(raw_counts_matrix="panc_combined_raw_matrix.txt",
                                    annotations_file="input_cnv_annotation.txt",
                                    delim="\t",
                                    gene_order_file="/Users/maitralab/Desktop/gencode_v19_gene_pos.txt",
                                    ref_group_names=NULL)


#making HMM false smooths the analysis over and makes things much faster 
#use 0.1 as cutoff for 10x data
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=.1,
                             out_dir="/Users/maitralab/Box/Single_Cell_10x_Data/Cell_Lines/3prime/panc2d3d_rna_122019",  
                             cluster_by_groups=T, 
                             denoise=T,
                             HMM=F)

##migrated to cluster and submitted analysis job there so as not to incapacitate computer. 

