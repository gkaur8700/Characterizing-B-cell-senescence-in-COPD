
#Load Libraries
library(Seurat) 
library(ggplot2)
library(BiocManager)
library(writexl)
library(devtools)
library(Azimuth)

# Using Azimuth annotations
lungref <- readRDS("path_to_folder/Azimuth_Reference.rds")

#check the dimensionality plot of the reference once
DimPlot(lungref, group.by = "ann_finest_level", raster = FALSE, label = TRUE, repel = TRUE) + NoLegend()

#Read the query for harmony integration
lungquery <- readRDS("path_to_folder/HarmonyIntegratedData.rds")

DimPlot(lungquery, raster = FALSE, label = TRUE, repel = TRUE) + NoLegend()

#Run Azimuth on the Harmony integration
lungquery <- RunAzimuth(lungquery, reference = "path_to_folder")

# Dimensionality Plots
Plot1 <- DimPlot(lungquery, group.by = "predicted.ann_finest_level", raster = FALSE, label = TRUE, repel = TRUE) + NoLegend() 
Plot2 <- DimPlot(lungquery, group.by = "predicted.ann_level_4", raster = FALSE, label = TRUE, repel = TRUE) 

Idents(lungquery) <- "predicted.ann_finest_level"
lungquery <- NormalizeData(lungquery)

# Save RDS
saveRDS(lungquery, "path_to_folder/HarmonyIntegrated_annotated_data.rds")

# Find Variable features after renaming
lungquery <- FindVariableFeatures(lungquery, selection.method = "vst", nfeatures = 3000)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones - Takes approximately a minute per cluster depending on the number of clusters.
lungquery <- FindAllMarkers(lungquery, only.pos = TRUE)
lungquery %>%
  group_by(cluster)

write_xlsx(lungquery, "path_to_folder/Marker_annotations.xlsx")

#Found clusters to remove based upon the cluster frequency information.
#These clusters have readings for 1-3 samples 
clusters_to_remove <- c("Fibromyocytes", "Tuft", "Goblet (nasal)", "SMG duct", "SMG serous (bronchial)", "Migratory DCs", "AT2 proliferating", "Multiciliated (nasal)")

#Cannot remove the alveolar M(phi) CCL3+ cell cluster
#Removed these cells
sobj_subset <- subset(x = annotated_data, idents = clusters_to_remove, invert = TRUE)

#Requires more tweaking
#Removing myofibroblasts, subrabasal cells and T-cell proliferating as they are over-represented in only one sample and could be an artifact

clusters_to_remove_1 <- c("Suprabasal", "T cells proliferating", "Myofibroblasts")

sobj_subset <- subset(x = sobj_subset, idents = clusters_to_remove_1, invert = TRUE)

# Find Variable features after renaming
sobj_subset <- FindVariableFeatures(sobj_subset, selection.method = "vst", nfeatures = 3000)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones - Takes approximately a minute per cluster depending on the number of clusters.
All.markers<- FindAllMarkers(sobj_subset, only.pos = TRUE)
All.markers %>%
  group_by(cluster)

write_xlsx(All.markers, "path_to_folder/Final_Marker_annotation.xlsx")

#Calculate the cell cluster frequency again
#This should be final

clusters_freq_2 <- sobj_subset@meta.data %>%
  group_by(sample, predicted.ann_finest_level) %>%
  summarise(n=n()) %>%
  mutate(relative_freq = n/sum(n))
clusters_freq_2$predicted.ann_finest_level <- factor(clusters_freq_2$predicted.ann_finest_level)

clusters_freq_2

write.csv(as.data.frame (clusters_freq_2),file="path_to_folder/Cluster_frequency_sample_wise.csv")

#save rds file

saveRDS(sobj_subset, "path_to_folder/HarmonyIntegrated_annotated_data.rds")
saveRDS(All.markers, "path_to_folder/annotated_markers.rds")
