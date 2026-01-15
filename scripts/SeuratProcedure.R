#SEURAT VERSION 4

library(Seurat)
library(dplyr)
#load matrix and create Seurat object
e125abcdTDT.data <- 
  Read10X(data.dir=
            "C:/Users/Jane/Documents/Jaworski Lab/TDTaggr/filtered_feature_bc_matrix/")
e125abcdTDT <- CreateSeuratObject(counts=e125abcdTDT.data, project="E12.5 ABCD TDT", 
                                  min.cells=3, min.features=200)
#Add percentMT and libraryID to metadata
e125abcdTDT[["percent.mt"]] <- PercentageFeatureSet(e125abcdTDT, pattern="^mt-")
libID <- sapply(strsplit(rownames(e125abcdTDT@meta.data), split="-"), "[[", 2)
e125abcdTDT <- AddMetaData(object = e125abcdTDT, 
                           metadata = data.frame(
                             libID=libID, 
                             row.names = rownames(e125abcdTDT@meta.data)))
e125abcdTDTUnfilt <- e125abcdTDT
saveRDS(e125abcdTDTUnfilt, file="~/Documents/Jaworski Lab/TDTaggr/e125abcdTDTUnfilt_Object.RData")
#Initial QC and filtering bad quality cells
vUnfilt <- VlnPlot(e125abcdTDT, 
                   features=c("nCount_RNA","nFeature_RNA","percent.mt"), 
                   group.by="libID", 
                   ncol=3)
#Filter out low-quality/dead cells and TDT-neg cells
e125abcdTDT <- subset(e125abcdTDT, subset = percent.mt<10 & nFeature_RNA >1000)
e125abcdTDT <- subset(e125abcdTDT, subset = tdTomato>0, slot="counts")
vFilt <- VlnPlot(e125abcdTDT,
                 features=c("nCount_RNA","nFeature_RNA","percent.mt"),
                 group.by="libID",
                 ncol=3, pt.size=0.2)
#Add sex based on Xist levels
sex <- factor(ifelse(e125abcdTDT$libID %in% c("1", "2"), "M", "F"))
e125abcdTDT <- AddMetaData(object = e125abcdTDT,
                           metadata = data.frame(
                             sex=sex,
                             row.names=rownames(e125abcdTDT@meta.data)))
#Normalize, find most variable features, scale
e125abcdTDT <- NormalizeData(e125abcdTDT)
e125abcdTDT <- FindVariableFeatures(e125abcdTDT)
top10 <- head(VariableFeatures(e125abcdTDT),10)
varfeat <- LabelPoints(plot=VariableFeaturePlot(e125abcdTDT),
                       points=top10,repel=TRUE)
e125abcdTDT <- ScaleData(e125abcdTDT, features=rownames(e125abcdTDT))
#PCA analysis and determine dimensionality
e125abcdTDT <- RunPCA(e125abcdTDT)
elbow <- ElbowPlot(e125abcdTDT)
#Clustering
e125abcdTDT <- FindNeighbors(e125abcdTDT, dims=1:16)
e125abcdTDT <- FindClusters(e125abcdTDT, resolution = 0.5)
#Run the commented line below only the first time you do this on your machine
#reticulate::py_install(packages='umap-learn')
e125abcdTDT <- RunUMAP(e125abcdTDT, dims = 1:20)
uPlot <- DimPlot(e125abcdTDT, reduction = "umap")
saveRDS(e125abcdTDT, file = "~/Documents/Jaworski Lab/TDTaggr/e125abcdTDT_Object.RData")



