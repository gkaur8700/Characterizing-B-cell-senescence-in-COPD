# Visium Final Analyses

#Load libraries
library(Seurat)
library(dplyr)
library(ggplot2)

#Load the samples
ns1 <- Load10X_Spatial(data.dir = "/...path_to_folder/Sample_1_346/outs", slice = "ns1")
ns2 <- Load10X_Spatial(data.dir = "/...path_to_folder/Sample_1_350/outs", slice = "ns2")
ns3 <- Load10X_Spatial(data.dir = "/...path_to_folder/Sample_2_350/outs", slice = "ns3")

sm1 <- Load10X_Spatial(data.dir = "/...path_to_folder/Sample_2_346/outs", slice = "sm1")

c1 <- Load10X_Spatial(data.dir = "/...path_to_folder/Sample_3_346/outs", slice = "c1")
c2 <- Load10X_Spatial(data.dir = "/...path_to_folder/Sample_3_350/outs", slice = "c2")
c3 <- Load10X_Spatial(data.dir = "/...path_to_folder/Sample_4_350/outs", slice = "c3")

#Data Preprocessing

#log normalization
ns1 <- NormalizeData(ns1)
ns2 <- NormalizeData(ns2)
ns3 <- NormalizeData(ns3)

sm1 <- NormalizeData(sm1)

c1 <- NormalizeData(c1)
c2 <- NormalizeData(c2)
c3 <- NormalizeData(c3)

#Merge all the samples

all.merge <- merge(ns1, y = c(ns2, ns3, sm1, c1, c2, c3), add.cell.ids = c("nonsmoker1", "nonsmoker2", "nonsmoker3", "smoker1", "copd1", "copd2", "copd3"))

#add a column with sample name
all.merge$sample <- sapply(strsplit(colnames(all.merge), "_"), `[`, 1)

# Make sure Spatial is the input assay
DefaultAssay(all.merge) <- "Spatial"

# Run SCTransform (this creates the SCT assay)
all.merge <- SCTransform(
  all.merge,
  assay = "Spatial",
  ncells = 3000,
  verbose = FALSE
)

# SCT is now the default
DefaultAssay(all.merge) <- "SCT"

# Dimensional reduction & clustering
all.merge <- RunPCA(all.merge, verbose = FALSE)
all.merge <- FindNeighbors(all.merge, dims = 1:30)
all.merge <- FindClusters(all.merge, resolution = 0.4)
all.merge <- RunUMAP(all.merge, dims = 1:30)

plot4 <- DimPlot(all.merge, reduction = "umap", group.by = c("ident", "sample"))

plot4

#Visualize the spatial plot
plot5 <- SpatialDimPlot(
  all.merge,
  label = TRUE,
  label.size = 3
)

plot5

#Marker Discovery

DefaultAssay(all.merge) <- "Spatial"

all.merge <- JoinLayers(all.merge)
markers <- FindAllMarkers(
  all.merge,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

#View Top markers
markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 5)

#Save the files
saveRDS(all.merge, "/path_to_folder/All_merged_after_integration.rds") 
saveRDS(markers, file = "/path_to_folder/cluster_markersrds")

#read the rds
data <- readRDS("path_to_folder/All_merged_after_integration.rds")

#rename the clusters
#based on the markers from scRNA seq

new_cluster_names <- c(
  "0" = "Endothelial Cells",
  "1" = "Monocytes",
  "2" = "B cells",
  "3" = "Smooth Muscle Cells",
  "4" = "Alveolar Macrophages",
  "5" = "Foam cells",
  "6" = "Reticulocyte",
  "7" = "Ciliated Epithelial Cells",
  "8" = "Transtional-AT2 cells",
  "9" = "AT2"
)

data <- RenameIdents(data, new_cluster_names)

Idents(data)

saveRDS(data, "/path_to_folder/All_merged_after_integration_annotated.rds")

```
