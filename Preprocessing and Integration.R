#Integrate the data

library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(harmony)

# Read the Data for Grp1
Sample1 <- Read10X(data.dir = "....path_to_folder/Saple ID/.../filtered_feature_bc_matrix")

# Make Seurat object
Sample1 <- CreateSeuratObject(counts = Sample1, project = "nonsmoker", assay = "RNA", min.features = 200, min.cells = 3)

#mitochondrial DNA
Sample1 <- PercentageFeatureSet(Sample1, pattern = "^MT-", col.name = "percent.mt")

#Filter out bad counts
Sample1 <- subset(Sample1, subset = nFeature_RNA > 200 & percent.mt < 30 & nCount_RNA < 1e+05)

#Standard Seurat pipeline
Sample1 <- NormalizeData(Sample1)
Sample1 <- FindVariableFeatures(Sample1)
Sample1 <- ScaleData(Sample1)
Sample1 <- RunPCA(Sample1)
Sample1 <- FindNeighbors(object = Sample1, dims = 1:20)
Sample1 <- FindClusters(object = Sample1)
Sample1 <- RunUMAP(Sample1, dims = 1:20, reduction = 'pca')

#Doublet Finder
Idents(Sample1)
sce <- scDblFinder(GetAssayData(Sample1, slot = "counts"), clusters = Idents (Sample1))
Sample1$scDblFinder.class <- sce$scDblFinder.class
table(Sample1@meta.data$scDblFinder.class)

#visualize the Doublets 
DimPlot(Sample1, reduction = 'umap', group.by = 'scDblFinder.class')
DimPlot(Sample1, reduction = 'umap', split.by = 'scDblFinder.class')

#Remove doublets
Idents(Sample1) <- Sample1@meta.data$scDblFinder.class
Idents(Sample1)
VlnPlot(Sample1, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
Sample1 <- subset(Sample1, idents = "singlet")
DimPlot(Sample1, reduction = 'umap', label = TRUE)

#Standard Workflow of Seurat
Sample1 <- NormalizeData(Sample1)
Sample1 <- FindVariableFeatures(Sample1)
Sample1 <- ScaleData(Sample1)
Sample1 <- RunPCA(Sample1)
Sample1 <- FindNeighbors(object = v, dims = 1:20)
Sample1 <- FindClusters(object = Sample1)
Sample1 <- RunUMAP(Sample1, dims = 1:20, reduction = 'pca')
DimPlot(Sample1, reduction = 'umap', label = TRUE)

#Save the rds file for NS2
saveRDS(Sample1, "path_to_folder/Sample1.rds")

# Repeat for all samples

# Merge all the samples
Sample.list <- list(Sample 1,....... Sample 13)
merged.samples <- merge(Sample 1, y = c(Sample 2,....., Sample 13 ), add.cell.ids = c("Sample1", ........, "Sample 13"), project = "human")

#Merge Sample of Non-Smokers
Grp1 <- list(Sample 1, Sample2, Sample3, Sample4)
merged.grp1 <- merge(Sample1, y = c(Sample2, Sample3, Sample4 ), add.cell.ids = c("Sample 1", "Sample2", "Sample3", "Sample4"), project = "human")

# Merge other groups

# Create metadata dataframe - This is so additions of new metrics will not interfere with main object

metadata <- merged.samples@meta.data

# Add cell IDs to metadata

metadata$cells <- rownames(metadata)

# Create sample column
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^Sample1_"))] <- "Grp1"
.
.
.
.
.
metadata$sample[which(str_detect(metadata$cells, "^sample13_"))] <- "Grp3"

merged.samples@meta.data <- metadata

#Feature Plot for quality assessment
d <- FeatureScatter(merged.samples, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')

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

# Run the standard workflow for visualization and clustering

merged.samples.harmony <- RunPCA(merged.samples.harmony, npcs = 50, verbose = TRUE)
merged.samples.harmony <- RunUMAP(merged.samples.harmony, reduction = "pca", dims = 1:50)
merged.samples.harmony <- FindNeighbors(merged.samples.harmony, reduction = "pca", dims = 1:50)
merged.samples.harmony <- FindClusters(merged.samples.harmony, resolution = c(0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.8, 2))

# Visualization

p1 <- DimPlot(merged.samples.harmony, reduction = "umap", group.by = "RNA_snn_res.1", label = TRUE)
p2 <- DimPlot(merged.samples.harmony, reduction = "umap", group.by = "RNA_snn_res.1", split.by = "orig.ident", label = TRUE)

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

```
