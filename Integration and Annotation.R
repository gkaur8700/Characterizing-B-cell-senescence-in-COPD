#Integrate the data

library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(harmony)

# Read the Data for Non Smokers_1
NS2 <- Read10X(data.dir = "....path_to_folder/Sample_IR_RL_021_2/.../filtered_feature_bc_matrix")

# Make Seurat object
NS2 <- CreateSeuratObject(counts = NS2, project = "nonsmoker", assay = "RNA", min.features = 200, min.cells = 3)

#mitochondrial DNA
NS2 <- PercentageFeatureSet(NS2, pattern = "^MT-", col.name = "percent.mt")

#Filter out bad counts
NS2 <- subset(NS2, subset = nFeature_RNA > 200 & percent.mt < 30 & nCount_RNA < 1e+05)

#Standard Seurat pipeline
NS2 <- NormalizeData(NS2)
NS2 <- FindVariableFeatures(NS2)
NS2 <- ScaleData(NS2)
NS2 <- RunPCA(NS2)
NS2 <- FindNeighbors(object = NS2, dims = 1:20)
NS2 <- FindClusters(object = NS2)
NS2 <- RunUMAP(NS2, dims = 1:20, reduction = 'pca')

#Doublet Finder
Idents(NS2)
sce <- scDblFinder(GetAssayData(NS2, slot = "counts"), clusters = Idents (NS2))
NS2$scDblFinder.class <- sce$scDblFinder.class
table(NS2@meta.data$scDblFinder.class)

#visualize the Doublets 
DimPlot(NS2, reduction = 'umap', group.by = 'scDblFinder.class')
DimPlot(NS2, reduction = 'umap', split.by = 'scDblFinder.class')

#Remove doublets
Idents(NS2) <- NS2@meta.data$scDblFinder.class
Idents(NS2)
VlnPlot(NS2, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
NS2 <- subset(NS2, idents = "singlet")
DimPlot(NS2, reduction = 'umap', label = TRUE)

#Standard Workflow of Seurat
NS2 <- NormalizeData(NS2)
NS2 <- FindVariableFeatures(NS2)
NS2 <- ScaleData(NS2)
NS2 <- RunPCA(NS2)
NS2 <- FindNeighbors(object = NS2, dims = 1:20)
NS2 <- FindClusters(object = NS2)
NS2 <- RunUMAP(NS2, dims = 1:20, reduction = 'pca')
DimPlot(NS2, reduction = 'umap', label = TRUE)

#Save the rds file for NS2
saveRDS(NS2, "path_to_folder/NS2.rds")

#Read data for Non-smoker_2
NS4 <- Read10X(data.dir = "....path_to_folder/Sample_IR_RL_028_2/.../filtered_feature_bc_matrix")
# Make Seurat object
NS4 <- CreateSeuratObject(counts = NS4, project = "nonsmoker", assay = "RNA", min.features = 200, min.cells = 3)

#mitochondrial DNA
NS4 <- PercentageFeatureSet(NS4, pattern = "^MT-", col.name = "percent.mt")
VlnPlot(NS4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Filter out bad counts
NS4 <- subset(NS4, subset = nFeature_RNA > 200 & percent.mt < 25 & nCount_RNA < 1e+05)
VlnPlot(NS4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


#Standard Seurat pipeline
NS4 <- NormalizeData(NS4)
NS4 <- FindVariableFeatures(NS4)
NS4 <- ScaleData(NS4)
NS4 <- RunPCA(NS4)
NS4 <- FindNeighbors(object = NS4, dims = 1:20)
NS4 <- FindClusters(object = NS4)
NS4 <- RunUMAP(NS4, dims = 1:20, reduction = 'pca')

#Doublet Finder
Idents(NS4)
sce1 <- scDblFinder(GetAssayData(NS4, slot = "counts"), clusters = Idents (NS4))
NS4$scDblFinder.class <- sce1$scDblFinder.class
table(NS4@meta.data$scDblFinder.class)

#visualize the Doublets 
DimPlot(NS4, reduction = 'umap', group.by = 'scDblFinder.class')
DimPlot(NS4, reduction = 'umap', split.by = 'scDblFinder.class')

#Remove doublets
Idents(NS4) <- NS4@meta.data$scDblFinder.class
Idents(NS4)
VlnPlot(NS4, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
NS4 <- subset(NS4, idents = "singlet")
DimPlot(NS4, reduction = 'umap', label = TRUE)

#Standard Workflow of Seurat
NS4 <- NormalizeData(NS4)
NS4 <- FindVariableFeatures(NS4)
NS4 <- ScaleData(NS4)
NS4 <- RunPCA(NS4)
NS4 <- FindNeighbors(object = NS4, dims = 1:20)
NS4 <- FindClusters(object = NS4)
NS4 <- RunUMAP(NS4, dims = 1:20, reduction = 'pca')
DimPlot(NS4, reduction = 'umap', label = TRUE)

#Save the rds file for NS2
saveRDS(NS4, "path_to_folder/NS4.rds")


#Read data Non-smoker_3
NS5 <- Read10X(data.dir = "....path_to_folder/Sample_IR_RL_031/.../filtered_feature_bc_matrix")
#Make Seurat object
NS5 <- CreateSeuratObject(counts = NS5, project = "nonsmoker", assay = "RNA", min.features = 200, min.cells = 3)

#mitochondrial DNA
NS5 <- PercentageFeatureSet(NS5, pattern = "^MT-", col.name = "percent.mt")
VlnPlot(NS5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Filter out bad counts
NS5 <- subset(NS5, subset = nFeature_RNA > 200 & percent.mt < 20 & nCount_RNA < 1e+05)
VlnPlot(NS5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


#Standard Seurat pipeline
NS5 <- NormalizeData(NS5)
NS5 <- FindVariableFeatures(NS5)
NS5 <- ScaleData(NS5)
NS5 <- RunPCA(NS5)
NS5 <- FindNeighbors(object = NS5, dims = 1:20)
NS5 <- FindClusters(object = NS5)
NS5 <- RunUMAP(NS5, dims = 1:20, reduction = 'pca')
DimPlot(NS5, reduction = 'umap', label = TRUE)

#Doublet Finder
Idents(NS5)
sce2 <- scDblFinder(GetAssayData(NS5, slot = "counts"), clusters = Idents (NS5))
NS5$scDblFinder.class <- sce2$scDblFinder.class
table(NS5@meta.data$scDblFinder.class)

#visualize the Doublets 
DimPlot(NS5, reduction = 'umap', group.by = 'scDblFinder.class')
DimPlot(NS5, reduction = 'umap', split.by = 'scDblFinder.class')

#Remove doublets
Idents(NS5) <- NS5@meta.data$scDblFinder.class
Idents(NS5)
VlnPlot(NS5, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
NS5 <- subset(NS5, idents = "singlet")
DimPlot(NS5, reduction = 'umap', label = TRUE)

#Standard Workflow of Seurat
NS5 <- NormalizeData(NS5)
NS5 <- FindVariableFeatures(NS5)
NS5 <- ScaleData(NS5)
NS5 <- RunPCA(NS5)
NS5 <- FindNeighbors(object = NS5, dims = 1:20)
NS5 <- FindClusters(object = NS5)
NS5 <- RunUMAP(NS5, dims = 1:20, reduction = 'pca')
DimPlot(NS4, reduction = 'umap', label = TRUE)

#Save the rds file for NS5
saveRDS(NS5, "path_to_folder/NS5.rds")

#Read Non-Smoker_4
NS6 <- Read10X(data.dir = "....path_to_folder/IR_RL_017/.../filtered_feature_bc_matrix_Sample17")

#Make Seurat object
NS6 <- CreateSeuratObject(counts = NS6, project = "nonsmoker", assay = "RNA", min.features = 200, min.cells = 3)

#mitochondrial DNA
NS6 <- PercentageFeatureSet(NS6, pattern = "^MT-", col.name = "percent.mt")
VlnPlot(NS6, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Filter out bad counts
NS6 <- subset(NS6, subset = nFeature_RNA > 200 & percent.mt < 30 & nCount_RNA < 1e+05)
VlnPlot(NS6, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Standard Seurat pipeline
NS6 <- NormalizeData(NS6)
NS6 <- FindVariableFeatures(NS6)
NS6 <- ScaleData(NS6)
NS6 <- RunPCA(NS6)
NS6 <- FindNeighbors(object = NS6, dims = 1:20)
NS6 <- FindClusters(object = NS6)
NS6 <- RunUMAP(NS6, dims = 1:20, reduction = 'pca')
DimPlot(NS6, reduction = 'umap', label = TRUE)

#Doublet Finder
Idents(NS6)
sce3 <- scDblFinder(GetAssayData(NS6, slot = "counts"), clusters = Idents (NS6))
NS6$scDblFinder.class <- sce3$scDblFinder.class
table(NS6@meta.data$scDblFinder.class)

#visualize the Doublets 
DimPlot(NS6, reduction = 'umap', group.by = 'scDblFinder.class')
DimPlot(NS6, reduction = 'umap', split.by = 'scDblFinder.class')

#Remove doublets
Idents(NS6) <- NS6@meta.data$scDblFinder.class
Idents(NS6)
VlnPlot(NS6, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
NS6 <- subset(NS6, idents = "singlet")
DimPlot(NS6, reduction = 'umap', label = TRUE)

#Standard Workflow of Seurat
NS6 <- NormalizeData(NS6)
NS6 <- FindVariableFeatures(NS6)
NS6 <- ScaleData(NS6)
NS6 <- RunPCA(NS6)
NS6 <- FindNeighbors(object = NS6, dims = 1:20)
NS6 <- FindClusters(object = NS6)
NS6 <- RunUMAP(NS6, dims = 1:20, reduction = 'pca')
DimPlot(NS6, reduction = 'umap', label = TRUE)

#Save the rds file for NS6
saveRDS(NS6, "path_to_folder/NS6.rds")

# Read the Data for Smokers
#Read Data for smoker_1
Sm1 <- Read10X(data.dir = "....path_to_folder/Sample_IR_LL_032/.../filtered_feature_bc_matrix")

#Make Seurat Object
Sm1 <- CreateSeuratObject(counts = Sm1, project = "smoker", assay = "RNA", min.features = 200, min.cells = 3)

#mitochondrial DNA
Sm1 <- PercentageFeatureSet(Sm1, pattern = "^MT-", col.name = "percent.mt")
VlnPlot(Sm1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Filter out bad counts
Sm1 <- subset(Sm1, subset = nFeature_RNA > 200 & percent.mt < 10 & nCount_RNA < 2e+05)
VlnPlot(Sm1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Standard Seurat pipeline
Sm1 <- NormalizeData(Sm1)
Sm1 <- FindVariableFeatures(Sm1)
Sm1 <- ScaleData(Sm1)
Sm1 <- RunPCA(Sm1)
Sm1 <- FindNeighbors(object = Sm1, dims = 1:20)
Sm1 <- FindClusters(object = Sm1)
Sm1 <- RunUMAP(Sm1, dims = 1:20, reduction = 'pca')
DimPlot(Sm1, reduction = 'umap', label = TRUE)

#Doublet Finder
Idents(Sm1)
sce4 <- scDblFinder(GetAssayData(Sm1, slot = "counts"), clusters = Idents (Sm1))
Sm1$scDblFinder.class <- sce4$scDblFinder.class
table(Sm1@meta.data$scDblFinder.class)

#visualize the Doublets 
DimPlot(Sm1, reduction = 'umap', group.by = 'scDblFinder.class')
DimPlot(Sm1, reduction = 'umap', split.by = 'scDblFinder.class')

#Remove doublets
Idents(Sm1) <- Sm1@meta.data$scDblFinder.class
Idents(Sm1)
VlnPlot(Sm1, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
Sm1 <- subset(Sm1, idents = "singlet")
DimPlot(Sm1, reduction = 'umap', label = TRUE)

#Standard Workflow of Seurat
Sm1 <- NormalizeData(Sm1)
Sm1 <- FindVariableFeatures(Sm1)
Sm1 <- ScaleData(Sm1)
Sm1 <- RunPCA(Sm1)
Sm1 <- FindNeighbors(object = Sm1, dims = 1:20)
Sm1 <- FindClusters(object = Sm1)
Sm1 <- RunUMAP(Sm1, dims = 1:20, reduction = 'pca')
DimPlot(Sm1, reduction = 'umap', label = TRUE)

#Save the rds file for NS6
saveRDS(Sm1, "path_to_folder/Sm1.rds")

#Read Smoker_2
Sm2 <- Read10X(data.dir = "....path_to_folder/Sample_IR_RL_034_2/.../filtered_feature_bc_matrix")

#Create a Seurat object
Sm2 <- CreateSeuratObject(counts = Sm2, project = "smoker", assay = "RNA", min.features = 200, min.cells = 3)

#mitochondrial DNA
Sm2 <- PercentageFeatureSet(Sm2, pattern = "^MT-", col.name = "percent.mt")
VlnPlot(Sm2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Filter out bad counts
Sm2 <- subset(Sm2, subset = nFeature_RNA > 200 & percent.mt < 25 & nCount_RNA < 1e+05)
VlnPlot(Sm2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Standard Seurat pipeline
Sm2 <- NormalizeData(Sm2)
Sm2 <- FindVariableFeatures(Sm2)
Sm2 <- ScaleData(Sm2)
Sm2 <- RunPCA(Sm2)
Sm2 <- FindNeighbors(object = Sm2, dims = 1:20)
Sm2 <- FindClusters(object = Sm2)
Sm2 <- RunUMAP(Sm2, dims = 1:20, reduction = 'pca')
DimPlot(Sm2, reduction = 'umap', label = TRUE)

#Doublet Finder
Idents(Sm2)
sce5 <- scDblFinder(GetAssayData(Sm2, slot = "counts"), clusters = Idents (Sm2))
Sm2$scDblFinder.class <- sce5$scDblFinder.class
table(Sm2@meta.data$scDblFinder.class)

#visualize the Doublets 
DimPlot(Sm2, reduction = 'umap', group.by = 'scDblFinder.class')
DimPlot(Sm2, reduction = 'umap', split.by = 'scDblFinder.class')

#Remove doublets
Idents(Sm2) <- Sm2@meta.data$scDblFinder.class
Idents(Sm2)
VlnPlot(Sm2, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
Sm2 <- subset(Sm2, idents = "singlet")
DimPlot(Sm2, reduction = 'umap', label = TRUE)

#Standard Workflow of Seurat
Sm2 <- NormalizeData(Sm2)
Sm2 <- FindVariableFeatures(Sm2)
Sm2 <- ScaleData(Sm2)
Sm2 <- RunPCA(Sm2)
Sm2 <- FindNeighbors(object = Sm2, dims = 1:20)
Sm2 <- FindClusters(object = Sm2)
Sm2 <- RunUMAP(Sm2, dims = 1:20, reduction = 'pca')
DimPlot(Sm2, reduction = 'umap', label = TRUE)

#Save the rds file for NS6
saveRDS(Sm2, "path_to_folder/Sm2.rds")

#Read the data Smoker_3
Sm3 <- Read10X(data.dir = "....path_to_folder/IR_LL_014/.../filtered_feature_bc_matrix")

#Create a Seurat object
Sm3 <- CreateSeuratObject(counts = Sm3, project = "smoker", assay = "RNA", min.features = 200, min.cells = 3)

#mitochondrial DNA
Sm3 <- PercentageFeatureSet(Sm3, pattern = "^MT-", col.name = "percent.mt")
VlnPlot(Sm3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Filter out bad counts
Sm3 <- subset(Sm3, subset = nFeature_RNA > 200 & percent.mt < 10 & nCount_RNA < 1e+05)
VlnPlot(Sm3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Standard Seurat pipeline
Sm3 <- NormalizeData(Sm3)
Sm3 <- FindVariableFeatures(Sm3)
Sm3 <- ScaleData(Sm3)
Sm3 <- RunPCA(Sm3)
Sm3 <- FindNeighbors(object = Sm3, dims = 1:20)
Sm3 <- FindClusters(object = Sm3)
Sm3 <- RunUMAP(Sm3, dims = 1:20, reduction = 'pca')
DimPlot(Sm3, reduction = 'umap', label = TRUE)

#Doublet Finder
Idents(Sm3)
sce6 <- scDblFinder(GetAssayData(Sm3, slot = "counts"), clusters = Idents (Sm3))
Sm3$scDblFinder.class <- sce6$scDblFinder.class
table(Sm3@meta.data$scDblFinder.class)

#visualize the Doublets 
DimPlot(Sm3, reduction = 'umap', group.by = 'scDblFinder.class')
DimPlot(Sm3, reduction = 'umap', split.by = 'scDblFinder.class')

#Remove doublets
Idents(Sm3) <- Sm3@meta.data$scDblFinder.class
Idents(Sm3)
VlnPlot(Sm3, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
Sm3 <- subset(Sm3, idents = "singlet")
DimPlot(Sm3, reduction = 'umap', label = TRUE)

#Standard Workflow of Seurat
Sm3 <- NormalizeData(Sm3)
Sm3 <- FindVariableFeatures(Sm3)
Sm3 <- ScaleData(Sm3)
Sm3 <- RunPCA(Sm3)
Sm3 <- FindNeighbors(object = Sm3, dims = 1:20)
Sm3 <- FindClusters(object = Sm3)
Sm3 <- RunUMAP(Sm3, dims = 1:20, reduction = 'pca')
DimPlot(Sm3, reduction = 'umap', label = TRUE)

#Save the rds file for NS6
saveRDS(Sm3, "path_to_folder/Sm3.rds")


#Read the data for Smoker_4
Sm4 <- Read10X(data.dir = "....path_to_folder/IR_RL_015/.../filtered_feature_bc_matrix")

#Create Seurat object
Sm4 <- CreateSeuratObject(counts = Sm4, project = "smoker", assay = "RNA", min.features = 200, min.cells = 3)

#mitochondrial DNA
Sm4 <- PercentageFeatureSet(Sm4, pattern = "^MT-", col.name = "percent.mt")
VlnPlot(Sm4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Filter out bad counts
Sm4 <- subset(Sm4, subset = nFeature_RNA > 200 & percent.mt < 30 & nCount_RNA < 1e+05)
VlnPlot(Sm4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Standard Seurat pipeline
Sm4 <- NormalizeData(Sm4)
Sm4 <- FindVariableFeatures(Sm4)
Sm4 <- ScaleData(Sm4)
Sm4 <- RunPCA(Sm4)
Sm4 <- FindNeighbors(object = Sm4, dims = 1:20)
Sm4 <- FindClusters(object = Sm4)
Sm4 <- RunUMAP(Sm4, dims = 1:20, reduction = 'pca')
DimPlot(Sm4, reduction = 'umap', label = TRUE)

#Doublet Finder
Idents(Sm4)
sce7 <- scDblFinder(GetAssayData(Sm4, slot = "counts"), clusters = Idents (Sm4))
Sm4$scDblFinder.class <- sce7$scDblFinder.class
table(Sm4@meta.data$scDblFinder.class)

#visualize the Doublets 
DimPlot(Sm4, reduction = 'umap', group.by = 'scDblFinder.class')
DimPlot(Sm4, reduction = 'umap', split.by = 'scDblFinder.class')

#Remove doublets
Idents(Sm4) <- Sm4@meta.data$scDblFinder.class
Idents(Sm4)
VlnPlot(Sm4, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
Sm4 <- subset(Sm4, idents = "singlet")
DimPlot(Sm4, reduction = 'umap', label = TRUE)

#Standard Workflow of Seurat
Sm4 <- NormalizeData(Sm4)
Sm4 <- FindVariableFeatures(Sm4)
Sm4 <- ScaleData(Sm4)
Sm4 <- RunPCA(Sm4)
Sm4 <- FindNeighbors(object = Sm4, dims = 1:20)
Sm4 <- FindClusters(object = Sm4)
Sm4 <- RunUMAP(Sm4, dims = 1:20, reduction = 'pca')
DimPlot(Sm4, reduction = 'umap', label = TRUE)

#Save the rds file for NS6
saveRDS(Sm4, "path_to_folder/Sm4.rds")


# Read the Data for COPD
#Read COPD_1
C1 <- Read10X(data.dir = "....path_to_folder/Sample_IR_RL_018/.../filtered_feature_bc_matrix")

#Make Seurat object
C1 <- CreateSeuratObject(counts = C1, project = "COPD", assay = "RNA", min.features = 200, min.cells = 3)

#mitochondrial DNA
C1 <- PercentageFeatureSet(C1, pattern = "^MT-", col.name = "percent.mt")
VlnPlot(C1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Filter out bad counts
C1 <- subset(C1, subset = nFeature_RNA > 200 & percent.mt < 50 & nCount_RNA < 1e+05)
VlnPlot(C1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Standard Seurat pipeline
C1 <- NormalizeData(C1)
C1 <- FindVariableFeatures(C1)
C1 <- ScaleData(C1)
C1 <- RunPCA(C1)
C1 <- FindNeighbors(object = C1, dims = 1:20)
C1 <- FindClusters(object = C1)
C1 <- RunUMAP(C1, dims = 1:20, reduction = 'pca')
DimPlot(C1, reduction = 'umap', label = TRUE)

#Doublet Finder
Idents(C1)
sce8 <- scDblFinder(GetAssayData(C1, slot = "counts"), clusters = Idents (C1))
C1$scDblFinder.class <- sce8$scDblFinder.class
table(C1@meta.data$scDblFinder.class)

#visualize the Doublets 
DimPlot(C1, reduction = 'umap', group.by = 'scDblFinder.class')
DimPlot(C1, reduction = 'umap', split.by = 'scDblFinder.class')

#Remove doublets
Idents(C1) <- C1@meta.data$scDblFinder.class
Idents(C1)
VlnPlot(C1, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
C1 <- subset(C1, idents = "singlet")
DimPlot(C1, reduction = 'umap', label = TRUE)

#Standard Workflow of Seurat
C1 <- NormalizeData(C1)
C1 <- FindVariableFeatures(C1)
C1 <- ScaleData(C1)
C1 <- RunPCA(C1)
C1 <- FindNeighbors(object = C1, dims = 1:20)
C1 <- FindClusters(object = C1)
C1 <- RunUMAP(C1, dims = 1:20, reduction = 'pca')
DimPlot(C1, reduction = 'umap', label = TRUE)

#Save the rds file for NS6
saveRDS(C1, "path_to_folder/C1.rds")

#Read COPD_2
C3 <- Read10X(data.dir = "....path_to_folder/Sample_IR_RL_023/.../filtered_feature_bc_matrix")

#Create Seurat object
C3 <- CreateSeuratObject(counts = C3, project = "COPD", assay = "RNA", min.features = 200, min.cells = 3)

#mitochondrial DNA
C3 <- PercentageFeatureSet(C3, pattern = "^MT-", col.name = "percent.mt")
VlnPlot(C3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Filter out bad counts
C3 <- subset(C3, subset = nFeature_RNA > 200 & percent.mt < 50 & nCount_RNA < 1e+05)
VlnPlot(C3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Standard Seurat pipeline
C3 <- NormalizeData(C3)
C3 <- FindVariableFeatures(C3)
C3 <- ScaleData(C3)
C3 <- RunPCA(C3)
C3 <- FindNeighbors(object = C3, dims = 1:20)
C3 <- FindClusters(object = C3)
C3 <- RunUMAP(C3, dims = 1:20, reduction = 'pca')
DimPlot(C3, reduction = 'umap', label = TRUE)

#Doublet Finder
Idents(C3)
sce9 <- scDblFinder(GetAssayData(C3, slot = "counts"), clusters = Idents (C3))
C3$scDblFinder.class <- sce9$scDblFinder.class
table(C3@meta.data$scDblFinder.class)

#visualize the Doublets 
DimPlot(C3, reduction = 'umap', group.by = 'scDblFinder.class')
DimPlot(C3, reduction = 'umap', split.by = 'scDblFinder.class')

#Remove doublets
Idents(C3) <- C3@meta.data$scDblFinder.class
Idents(C3)
VlnPlot(C3, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
C3 <- subset(C3, idents = "singlet")
DimPlot(C3, reduction = 'umap', label = TRUE)

#Standard Workflow of Seurat
C3 <- NormalizeData(C3)
C3 <- FindVariableFeatures(C3)
C3 <- ScaleData(C3)
C3 <- RunPCA(C3)
C3 <- FindNeighbors(object = C3, dims = 1:20)
C3 <- FindClusters(object = C3)
C3 <- RunUMAP(C3, dims = 1:20, reduction = 'pca')
DimPlot(C3, reduction = 'umap', label = TRUE)

#Save the rds file for NS6
saveRDS(C3, "path_to_folder/C3.rds")

#Read the data for COPD_5

C5 <- Read10X(data.dir = "....path_to_folder/Sample_IR_RL_027_2/.../filtered_feature_bc_matrix")

#Make Seurat object
C5 <- CreateSeuratObject(counts = C5, project = "COPD", assay = "RNA", min.features = 200, min.cells = 3)

#mitochondrial DNA
C5 <- PercentageFeatureSet(C5, pattern = "^MT-", col.name = "percent.mt")
VlnPlot(C5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Filter out bad counts
C5 <- subset(C5, subset = nFeature_RNA > 200 & percent.mt < 50 & nCount_RNA < 1e+05)
VlnPlot(C5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Standard Seurat pipeline
C5 <- NormalizeData(C5)
C5 <- FindVariableFeatures(C5)
C5 <- ScaleData(C5)
C5 <- RunPCA(C5)
C5 <- FindNeighbors(object = C5, dims = 1:20)
C5 <- FindClusters(object = C5)
C5 <- RunUMAP(C5, dims = 1:20, reduction = 'pca')
DimPlot(C5, reduction = 'umap', label = TRUE)

#Doublet Finder
Idents(C5)
sce10 <- scDblFinder(GetAssayData(C5, slot = "counts"), clusters = Idents (C5))
C5$scDblFinder.class <- sce10$scDblFinder.class
table(C5@meta.data$scDblFinder.class)

#visualize the Doublets 
DimPlot(C5, reduction = 'umap', group.by = 'scDblFinder.class')
DimPlot(C5, reduction = 'umap', split.by = 'scDblFinder.class')

#Remove doublets
Idents(C5) <- C5@meta.data$scDblFinder.class
Idents(C5)
VlnPlot(C5, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
C5 <- subset(C5, idents = "singlet")
DimPlot(C5, reduction = 'umap', label = TRUE)

#Standard Workflow of Seurat
C5 <- NormalizeData(C5)
C5 <- FindVariableFeatures(C5)
C5 <- ScaleData(C5)
C5 <- RunPCA(C5)
C5 <- FindNeighbors(object = C5, dims = 1:20)
C5 <- FindClusters(object = C5)
C5 <- RunUMAP(C5, dims = 1:20, reduction = 'pca')
DimPlot(C5, reduction = 'umap', label = TRUE)

#Save the rds file for NS6
saveRDS(C5, "path_to_folder/C5.rds")

#Read the data for COPD_6
C6 <- Read10X(data.dir = "....path_to_folder/Sample_IR_RL_029/.../filtered_feature_bc_matrix")

#Make Seurat object
C6 <- CreateSeuratObject(counts = C6, project = "COPD", assay = "RNA", min.features = 200, min.cells = 3)

#mitochondrial DNA
C6 <- PercentageFeatureSet(C6, pattern = "^MT-", col.name = "percent.mt")
VlnPlot(C6, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Filter out bad counts
C6 <- subset(C6, subset = nFeature_RNA > 200 & percent.mt < 50 & nCount_RNA < 1e+05)
VlnPlot(C6, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Standard Seurat pipeline
C6 <- NormalizeData(C6)
C6 <- FindVariableFeatures(C6)
C6 <- ScaleData(C6)
C6 <- RunPCA(C6)
C6 <- FindNeighbors(object = C6, dims = 1:20)
C6 <- FindClusters(object = C6)
C6 <- RunUMAP(C6, dims = 1:20, reduction = 'pca')
DimPlot(C6, reduction = 'umap', label = TRUE)

#Doublet Finder
Idents(C6)
sce11 <- scDblFinder(GetAssayData(C6, slot = "counts"), clusters = Idents (C6))
C6$scDblFinder.class <- sce11$scDblFinder.class
table(C6@meta.data$scDblFinder.class)

#visualize the Doublets 
DimPlot(C6, reduction = 'umap', group.by = 'scDblFinder.class')
DimPlot(C6, reduction = 'umap', split.by = 'scDblFinder.class')

#Remove doublets
Idents(C6) <- C6@meta.data$scDblFinder.class
Idents(C6)
VlnPlot(C6, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
C6 <- subset(C6, idents = "singlet")
DimPlot(C6, reduction = 'umap', label = TRUE)

#Standard Workflow of Seurat
C6 <- NormalizeData(C6)
C6 <- FindVariableFeatures(C6)
C6 <- ScaleData(C6)
C6 <- RunPCA(C6)
C6 <- FindNeighbors(object = C6, dims = 1:20)
C6 <- FindClusters(object = C6)
C6 <- RunUMAP(C6, dims = 1:20, reduction = 'pca')
DimPlot(C6, reduction = 'umap', label = TRUE)

#Save the rds file for NS6
saveRDS(C6, "path_to_folder/C6.rds")

#Read the data for COPD_7
C7 <- Read10X(data.dir = "....path_to_folder/IR_RL_016/.../filtered_feature_bc_matrix_Sample16")

# Make Seurat object
C7 <- CreateSeuratObject(counts = C7, project = "COPD", assay = "RNA", min.features = 200, min.cells = 3)

#mitochondrial DNA
C7 <- PercentageFeatureSet(C7, pattern = "^MT-", col.name = "percent.mt")
VlnPlot(C7, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Filter out bad counts
C7 <- subset(C7, subset = nFeature_RNA > 200 & percent.mt < 40 & nCount_RNA < 1e+05)
VlnPlot(C7, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Standard Seurat pipeline
C7 <- NormalizeData(C7)
C7 <- FindVariableFeatures(C7)
C7 <- ScaleData(C7)
C7 <- RunPCA(C7)
C7 <- FindNeighbors(object = C7, dims = 1:20)
C7 <- FindClusters(object = C7)
C7 <- RunUMAP(C7, dims = 1:20, reduction = 'pca')
DimPlot(C7, reduction = 'umap', label = TRUE)

#Doublet Finder
Idents(C7)
sce12 <- scDblFinder(GetAssayData(C7, slot = "counts"), clusters = Idents (C7))
C7$scDblFinder.class <- sce12$scDblFinder.class
table(C7@meta.data$scDblFinder.class)

#visualize the Doublets 
DimPlot(C7, reduction = 'umap', group.by = 'scDblFinder.class')
DimPlot(C7, reduction = 'umap', split.by = 'scDblFinder.class')

#Remove doublets
Idents(C7) <- C7@meta.data$scDblFinder.class
Idents(C7)
VlnPlot(C7, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
C7 <- subset(C7, idents = "singlet")
DimPlot(C7, reduction = 'umap', label = TRUE)

#Standard Workflow of Seurat
C7 <- NormalizeData(C7)
C7 <- FindVariableFeatures(C7)
C7 <- ScaleData(C7)
C7 <- RunPCA(C7)
C7 <- FindNeighbors(object = C7, dims = 1:20)
C7 <- FindClusters(object = C7)
C7 <- RunUMAP(C7, dims = 1:20, reduction = 'pca')
DimPlot(C7, reduction = 'umap', label = TRUE)

#Save the rds file for NS6
saveRDS(C7, "path_to_folder/C7.rds")

# Merge all the samples
Sample.list <- list(NS2, NS4, NS5, NS6, Sm1, Sm2, Sm3, Sm4, C1, C3, C5, C6, C7)
merged.samples <- merge(NS2, y = c(NS4, NS5, NS6,Sm1, Sm2, Sm3, Sm4, C1, C3, C5, C6, C7 ), add.cell.ids = c("NS2", "NS4", "NS5", "NS6", "Sm1", "Sm2","Sm3", "Sm4", "C1", "C3", "C5", "C6", "C7"), project = "human")

#Merge Sample of Non-Smokers
NonSmoker <- list(NS2, NS4, NS5, NS6)
merged.NS <- merge(NS2, y = c(NS4, NS5, NS6 ), add.cell.ids = c("NS2", "NS4", "NS5", "NS6"), project = "human")

# Merge Samples of Smokers
Smoker <- list(Sm1, Sm2, Sm3, Sm4)
merged.Sm <- merge(Sm1, y = c(Sm2, Sm3, Sm4), add.cell.ids = c("Sm1", "Sm2","Sm3", "Sm4"), project = "human")

# Merge samples of COPD
COPD <- list(C1, C3, C5, C6, C7)
merged.copds <- merge(C1, y = c(C3, C5, C6, C7), add.cell.ids = c("C1", "C3", "C5", "C6", "C7"), project = "human")

# Create metadata dataframe - This is so additions of new metrics will not interfere with main object

metadata <- merged.samples@meta.data

# Add cell IDs to metadata

metadata$cells <- rownames(metadata)

# Create sample column
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^NS2_"))] <- "NS2"
metadata$sample[which(str_detect(metadata$cells, "^NS4_"))] <- "NS4"
metadata$sample[which(str_detect(metadata$cells, "^NS5_"))] <- "NS5"
metadata$sample[which(str_detect(metadata$cells, "^NS6_"))] <- "NS6"

metadata$sample[which(str_detect(metadata$cells, "^Sm1_"))] <- "Sm1"
metadata$sample[which(str_detect(metadata$cells, "^Sm2_"))] <- "Sm2"
metadata$sample[which(str_detect(metadata$cells, "^Sm3_"))] <- "Sm3"
metadata$sample[which(str_detect(metadata$cells, "^Sm4_"))] <- "Sm4"

metadata$sample[which(str_detect(metadata$cells, "^C1_"))] <- "C1"
metadata$sample[which(str_detect(metadata$cells, "^C3_"))] <- "C3"
metadata$sample[which(str_detect(metadata$cells, "^C5_"))] <- "C5"
metadata$sample[which(str_detect(metadata$cells, "^C6_"))] <- "C6"
metadata$sample[which(str_detect(metadata$cells, "^C7_"))] <- "C7"

merged.samples@meta.data <- metadata

e <- VlnPlot(merged.samples, features = c("nFeature_RNA"), cols = c("orange", "blue", "red", "green"), group.by = "orig.ident")
f <- VlnPlot(merged.samples, features = c("nCount_RNA"), cols = c("orange", "blue", "red", "green"), group.by = "orig.ident")
g <- VlnPlot(merged.samples, features = c("percent.mt"), cols = c("orange", "blue", "red", "green"), group.by = "orig.ident")


# Save the Plots

ggsave(e, filename = 'path_to_folder/nFeatures_postFiltration.tiff', width = 10, height =8)

ggsave(f, filename = 'path_to_folder/nRNA_postFiltration.tiff', width = 10, height =8)

ggsave(g, filename = 'path_tpo_folder/mtPercent_postFiltration.tiff', width = 10, height =8)

#Feature Plot for quality assessment
d <- FeatureScatter(merged.samples, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')

ggsave(d, filename = 'path_to_folder/Feature Scatter Plot_afterdoubletremoval.tiff', width = 10, height =8)

#Follow standard Seurat workflow
# 3. Normalize data ----------

merged.samples <- NormalizeData(merged.samples)

# 4. Identify highly variable features --------------
merged.samples <- FindVariableFeatures(merged.samples)

# 5. Scaling -------------
#all.genes <- rownames(merged.samples_postfilter)
merged.samples <- ScaleData(merged.samples)

# 6. Perform Linear dimensionality reduction --------------
merged.samples <- RunPCA(merged.samples)

ElbowPlot(merged.samples)

merged.samples <- FindNeighbors(object = merged.samples, dims = 1:20)
merged.samples <- FindClusters(object = merged.samples)
merged.samples <- RunUMAP(merged.samples, dims = 1:20, reduction = 'pca')

#Plot before integration
#Plot1 <- DimPlot(merged.samples, reduction = 'umap', split.by = 'orig.ident', label = TRUE)
Plot2 <- DimPlot(merged.samples, reduction = 'umap', group.by = 'RNA_snn_res.0.8', label = TRUE)

ggsave(Plot2, filename = 'path_to_folder/DimensionalityPlot_beforeIntegration.tiff', width = 10, height =8)

saveRDS(merged.samples, "path_to_folder/MergedData_after normalization.rds")

# run Harmony -----------
merged.samples.harmony <- merged.samples %>%
  RunHarmony(group.by.vars = 'orig.ident', plot_convergence = FALSE)

merged.samples.harmony@reductions

merged.samples.harmony.embed <- Embeddings(merged.samples.harmony, "harmony")
merged.samples.harmony.embed[1:10,1:10]

# Do UMAP and clustering using ** Harmony embeddings instead of PCA **
merged.samples.harmony <- merged.samples.harmony %>%
  RunUMAP(reduction = 'harmony', dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.5)

saveRDS(merged.samples.harmony, "path_to_folder/HarmonyIntegratedData.rds") 

library(writexl) 
write_xlsx(merged.samples.harmony.embed, "D:/April12-2024/HarmonyIntegratedData_enbed.xlsx")

i <- DimPlot(merged.samples.harmony, reduction = 'umap', group.by = 'seurat_clusters')

# Run the standard workflow for visualization and clustering

merged.samples.harmony <- RunPCA(merged.samples.harmony, npcs = 50, verbose = TRUE)
merged.samples.harmony <- RunUMAP(merged.samples.harmony, reduction = "pca", dims = 1:50)
merged.samples.harmony <- FindNeighbors(merged.samples.harmony, reduction = "pca", dims = 1:50)
merged.samples.harmony <- FindClusters(merged.samples.harmony, resolution = c(0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.8, 2))

# Visualization

p1 <- DimPlot(merged.samples.harmony, reduction = "umap", group.by = "RNA_snn_res.1", label = TRUE)
ggsave(p1, filename = 'D:/April12-2024/harmonyintegration_resolution1_UMAP.tiff', width = 10, height =8)

p2 <- DimPlot(merged.samples.harmony, reduction = "umap", group.by = "RNA_snn_res.1", split.by = "orig.ident", label = TRUE)
ggsave(p2, filename = 'D:/April12-2024/harmonyintegration_resolution1_UMAP_by groups.tiff', width = 10, height =8)

# Assign identity of clusters

Idents(object = merged.samples.harmony) <- "RNA_snn_res.1"

#Normalize the data once again
merged.samples.harmony <- NormalizeData(merged.samples.harmony)
merged.samples.harmony <- ScaleData(merged.samples.harmony)
merged.samples.harmony <- FindVariableFeatures(merged.samples.harmony, selection.method = "vst", nfeatures = 3000)

#Join the layers
merged.samples.harmony <- JoinLayers(merged.samples.harmony)

#Identify markers
merged.markers <- FindAllMarkers(merged.samples.harmony, only.pos = TRUE)
merged.markers %>%
  group_by(cluster)

write_xlsx(merged.markers, "path_to_folder/harmony_integrate_preannotation.xlsx")
saveRDS(merged.samples.harmony, "path_to_folder/HarmonyIntegratedData.rds")

