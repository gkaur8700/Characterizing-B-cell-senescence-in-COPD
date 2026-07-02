# April 5-2024
#Integrate the data

library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(harmony)


# Read the Data for Non Smokers
NS2 <- Read10X(data.dir = "E:/deliv_Rahman_10X_100422_counts (2).tar/deliv_Rahman_10X_100422_counts (2)/deliv_Rahman_10X_100422_counts/Sample_IR_RL_021_2/outs/filtered_feature_bc_matrix")
NS4 <- Read10X(data.dir = "E:/deliv_Rahman_10X_100422_counts (2).tar/deliv_Rahman_10X_100422_counts (2)/deliv_Rahman_10X_100422_counts/Sample_IR_RL_028_2/outs/filtered_feature_bc_matrix")
NS5 <- Read10X(data.dir = "E:/deliv_Rahman_10X_100422_counts (2).tar/deliv_Rahman_10X_100422_counts (2)/deliv_Rahman_10X_100422_counts/Sample_IR_RL_031/outs/filtered_feature_bc_matrix")
NS6 <- Read10X(data.dir = "F:/Data/filtered_feature_bc_matrix_Sample17")

# Read the Data for Smokers

Sm1 <- Read10X(data.dir = "E:/deliv_Rahman_10X_100422_counts (2).tar/deliv_Rahman_10X_100422_counts (2)/deliv_Rahman_10X_100422_counts/Sample_IR_LL_032/outs/filtered_feature_bc_matrix")
Sm2 <- Read10X(data.dir = "E:/deliv_Rahman_10X_100422_counts (2).tar/deliv_Rahman_10X_100422_counts (2)/deliv_Rahman_10X_100422_counts/Sample_IR_RL_034_2/outs/filtered_feature_bc_matrix")
Sm3 <- Read10X(data.dir = "E:/deliv_Rahman_10X_100422_counts (2).tar/deliv_Rahman_10X_100422_counts (2)/deliv_Rahman_10X_100422_counts/Sample_Liberase_LUL/outs/filtered_feature_bc_matrix")
Sm4 <- Read10X(data.dir = "F:/Data/filtered_feature_bc_matrix_Sample15")

# Read the Data for COPD

C1 <- Read10X(data.dir = "E:/deliv_Rahman_10X_100422_counts (2).tar/deliv_Rahman_10X_100422_counts (2)/deliv_Rahman_10X_100422_counts/Sample_IR_RL_018/outs/filtered_feature_bc_matrix")
C3 <- Read10X(data.dir = "E:/deliv_Rahman_10X_100422_counts (2).tar/deliv_Rahman_10X_100422_counts (2)/deliv_Rahman_10X_100422_counts/Sample_IR_RL_023/outs/filtered_feature_bc_matrix")
C5 <- Read10X(data.dir = "E:/deliv_Rahman_10X_100422_counts (2).tar/deliv_Rahman_10X_100422_counts (2)/deliv_Rahman_10X_100422_counts/Sample_IR_RL_027_2/outs/filtered_feature_bc_matrix")
C6 <- Read10X(data.dir = "E:/deliv_Rahman_10X_100422_counts (2).tar/deliv_Rahman_10X_100422_counts (2)/deliv_Rahman_10X_100422_counts/Sample_IR_RL_029/outs/filtered_feature_bc_matrix")
C7 <- Read10X(data.dir = "F:/Data/filtered_feature_bc_matrix_Sample16")

# Make Seurat object
NS2 <- CreateSeuratObject(counts = NS2, project = "nonsmoker", assay = "RNA", min.features = 200, min.cells = 3)
NS4 <- CreateSeuratObject(counts = NS4, project = "nonsmoker", assay = "RNA", min.features = 200, min.cells = 3)
NS5 <- CreateSeuratObject(counts = NS5, project = "nonsmoker", assay = "RNA", min.features = 200, min.cells = 3)
NS6 <- CreateSeuratObject(counts = NS6, project = "nonsmoker", assay = "RNA", min.features = 200, min.cells = 3)

Sm1 <- CreateSeuratObject(counts = Sm1, project = "smoker", assay = "RNA", min.features = 200, min.cells = 3)
Sm2 <- CreateSeuratObject(counts = Sm2, project = "smoker", assay = "RNA", min.features = 200, min.cells = 3)
Sm3 <- CreateSeuratObject(counts = Sm3, project = "smoker", assay = "RNA", min.features = 200, min.cells = 3)
Sm4 <- CreateSeuratObject(counts = Sm4, project = "smoker", assay = "RNA", min.features = 200, min.cells = 3)

C1 <- CreateSeuratObject(counts = C1, project = "COPD", assay = "RNA", min.features = 200, min.cells = 3)
C3 <- CreateSeuratObject(counts = C3, project = "COPD", assay = "RNA", min.features = 200, min.cells = 3)
C5 <- CreateSeuratObject(counts = C5, project = "COPD", assay = "RNA", min.features = 200, min.cells = 3)
C6 <- CreateSeuratObject(counts = C6, project = "COPD", assay = "RNA", min.features = 200, min.cells = 3)
C7 <- CreateSeuratObject(counts = C7, project = "COPD", assay = "RNA", min.features = 200, min.cells = 3)

# Merge all the samples
Sample.list <- list(NS2, NS4, NS5, NS6, Sm1, Sm2, Sm3, Sm4, C1, C3, C5, C6, C7)
merged.samples <- merge(NS2, y = c(NS4, NS5, NS6,Sm1, Sm2, Sm3, Sm4, C1, C3, C5, C6, C7 ), add.cell.ids = c("NS2", "NS4", "NS5", "NS6", "Sm1", "Sm2","Sm3", "Sm4", "C1", "C3", "C5", "C6", "C7"), project = "human")

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

#Filtering and data cleaning#
# Add mitochondrial gene fraction that will be used to exclude dying cells
#Use capital letters for human data ^MT
merged.samples <- PercentageFeatureSet(merged.samples, pattern = "^MT-", col.name = "percent.mt")

#Use Feature Scatter for further quality assessment
d <- FeatureScatter(merged.samples, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')

#split the list for filtration
split_samples_prefilter <- SplitObject(merged.samples, split.by = "sample")

# Filter out Unwanted cell/debris/multiplets 

split_NS2 <- subset(split_samples_prefilter$NS2, subset = nFeature_RNA > 200 & percent.mt < 30 & nCount_RNA < 1e+05)
split_NS4 <- subset(split_samples_prefilter$NS4, subset = nFeature_RNA > 200 & percent.mt < 30 & nCount_RNA < 1e+05)
split_NS5 <- subset(split_samples_prefilter$NS5, subset = nFeature_RNA > 200 & percent.mt < 30 & nCount_RNA < 1e+05)
split_NS6 <- subset(split_samples_prefilter$NS6, subset = nFeature_RNA > 200 & percent.mt < 30 & nCount_RNA < 1e+05)

split_Sm1 <- subset(split_samples_prefilter$Sm1, subset = nFeature_RNA > 200 & percent.mt < 30 & nCount_RNA < 1e+05)
split_Sm2 <- subset(split_samples_prefilter$Sm2, subset = nFeature_RNA > 200 & percent.mt < 30 & nCount_RNA < 1e+05)
split_Sm3 <- subset(split_samples_prefilter$Sm3, subset = nFeature_RNA > 200 & percent.mt < 30 & nCount_RNA < 1e+05)
split_Sm4 <- subset(split_samples_prefilter$Sm4, subset = nFeature_RNA > 200 & percent.mt < 30 & nCount_RNA < 1e+05)

split_C1 <- subset(split_samples_prefilter$C1, subset = nFeature_RNA > 200 & percent.mt < 30 & nCount_RNA < 1e+05)
split_C3 <- subset(split_samples_prefilter$C3, subset = nFeature_RNA > 200 & percent.mt < 30 & nCount_RNA < 1e+05)
split_C5 <- subset(split_samples_prefilter$C5, subset = nFeature_RNA > 200 & percent.mt < 30 & nCount_RNA < 1e+05)
split_C6 <- subset(split_samples_prefilter$C6, subset = nFeature_RNA > 200 & percent.mt < 30 & nCount_RNA < 1e+05)
split_C7 <- subset(split_samples_prefilter$C7, subset = nFeature_RNA > 200 & percent.mt < 30 & nCount_RNA < 1e+05)

# Merge data after filtering 
merged.samples_postfilter <- merge(x = split_NS2, y = c(split_NS4, split_NS5, split_NS6, split_Sm1, split_Sm2, split_Sm3, split_Sm4, split_C1, split_C3, split_C5, split_C6, split_C7), add.cell.ids = c("NS2", "NS4", "NS5", "NS6", "Sm1", "Sm2", "Sm3", "Sm4", "C1", "C3", "C5", "C6", "C7"))

# Plot the feature plots post filtration

e <- VlnPlot(merged.samples_postfilter, features = c("nFeature_RNA"), cols = c("orange", "blue", "red", "green"), split.by = "sample")
f <- VlnPlot(merged.samples_postfilter, features = c("nCount_RNA"), cols = c("orange", "blue", "red", "green"), split.by = "sample")
g <- VlnPlot(merged.samples_postfilter, features = c("percent.mt"), cols = c("orange", "blue", "red", "green"), split.by = "sample")


#Save the the mitochondrial percentages

ggsave(g, filename = 'D:/Combined Results/mtPercent_postFiltration_mtdna_copd_30.tiff', width = 10, height =8)


# 3. Normalize data ----------

merged.samples_postfilter <- NormalizeData(merged.samples_postfilter)

# 4. Identify highly variable features --------------
merged.samples_postfilter <- FindVariableFeatures(merged.samples_postfilter)

# 5. Scaling -------------
#all.genes <- rownames(merged.samples_postfilter)
merged.samples_postfilter <- ScaleData(merged.samples_postfilter)

# 6. Perform Linear dimensionality reduction --------------
merged.samples_postfilter <- RunPCA(merged.samples_postfilter)

ElbowPlot(merged.samples_postfilter)

merged.samples_postfilter <- RunUMAP(merged.samples_postfilter, dims = 1:20, reduction = 'pca')

#Plot before integration
before <- DimPlot(merged.samples_postfilter, reduction = 'umap', group.by = 'orig.ident')

# run Harmony -----------
merged.samples.harmony <- merged.samples_postfilter %>%
  RunHarmony(group.by.vars = 'orig.ident', plot_convergence = FALSE)

merged.samples.harmony@reductions

merged.samples.harmony.embed <- Embeddings(merged.samples.harmony, "harmony")
merged.samples.harmony.embed[1:10,1:10]

# Do UMAP and clustering using ** Harmony embeddings instead of PCA **
merged.samples.harmony <- merged.samples.harmony %>%
  RunUMAP(reduction = 'harmony', dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.5)

# visualize 
after <- DimPlot(merged.samples.harmony, reduction = 'umap', group.by = 'orig.ident')

h <- before|after

ggsave(h, filename = 'D:/Combined Results/before and after harmony.tiff', width = 10, height =8)

saveRDS(merged.samples.harmony, "D:/Combined Results/HarmonyIntegratedData.rds") 

i <- DimPlot(merged.samples.harmony, reduction = 'umap', group.by = 'seurat_clusters')

ggsave(i, filename = 'D:/Combined Results/UMAP_after_harmony.tiff', width = 10, height =8)




#Integration Method#2 
#SC Tranform
#Split seurat object by condition to perform QC

merged.samples_preintegration <- SplitObject(merged.samples_postfilter, split.by = "orig.ident")
#merged.samples_preintegration1 <- SplitObject(merged.samples_postfilter, split.by = "sample")

#Integrations of dataset#
# normalize and identify variable features for each dataset independently

split_mergedsamples_filtered <- lapply(X = merged.samples_preintegration, FUN = function(x) 
{x <- SCTransform(x, method = "glmGamPoi", variable.features.n = 3000, return.only.var.genes = TRUE)})

#split_mergedsamples_filtered1 <- lapply(X = merged.samples_preintegration1, FUN = function(x) 
#{x <- SCTransform(x, method = "glmGamPoi", variable.features.n = 3000, return.only.var.genes = TRUE)})

# select features that are repeatedly variable across datasets for integration

features <- SelectIntegrationFeatures(object.list = split_mergedsamples_filtered, nfeatures = 3000)
#features1 <- SelectIntegrationFeatures(object.list = split_mergedsamples_filtered1, nfeatures = 3000)

# Prepare the SCT list object for integration

split_samples_prep <- PrepSCTIntegration(object.list = split_mergedsamples_filtered, anchor.features = features)

#split_samples_prep1 <- PrepSCTIntegration(object.list = split_mergedsamples_filtered1, anchor.features = features1)

# Find anchor for integration

split_sample_anchor <- FindIntegrationAnchors(object.list = split_samples_prep1, anchor.features = features1, dims = 1:40, normalization.method = "SCT")

integrated_sample <- IntegrateData(anchorset = split_sample_anchor, dims = 1:40, normalization.method = "SCT")

# Save the Data

saveRDS(integrated_sample, "D:/Combined Results/SCT_integration.rds") 


DefaultAssay(integrated_ecig) <- "integrated"

# Run the standard workflow for visualization and clustering

integrated_sample <- RunPCA(integrated_sample, npcs = 40, verbose = TRUE)
integrated_sample <- RunUMAP(integrated_sample, reduction = "pca", dims = 1:40)
integrated_sample <- FindNeighbors(integrated_sample, reduction = "pca", dims = 1:40)
integrated_sample <- FindClusters(integrated_sample, resolution = c(0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.8, 2))

# Assign identity of clusters

Idents(object = integrated_sample) <- "integrated_snn_res.1"

# Visualization

p1 <- DimPlot(integrated_sample, reduction = "umap", group.by = "integrated_snn_res.1")
ggsave(p1, filename = 'D:/Combined Results/SCTintegration_UMAP.tiff', width = 10, height =8)






