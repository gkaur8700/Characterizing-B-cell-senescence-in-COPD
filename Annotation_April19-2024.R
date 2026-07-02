
library(Seurat) 
library(ggplot2)
library(tidyverse)
library(BiocManager)
library(writexl)
library(devtools)
#Make sure the TBFSTools are updated
# https://github.com/satijalab/seurat/issues/8202
#install Azimuth
devtools::install_github("satijalab/azimuth","seurat5", quiet = TRUE)

library(Azimuth)

# Using Azimuth annotations
lungref <- readRDS("D:/Azimuth Reference/ref.rds")

#check the dimensionality plot of the reference once
DimPlot(lungref, group.by = "ann_finest_level", raster = FALSE, label = TRUE, repel = TRUE) + NoLegend()

#Read the query for harmony integration

lungquery <- readRDS("D:/April12-2024/HarmonyIntegratedData.rds")

DimPlot(lungquery, raster = FALSE, label = TRUE, repel = TRUE) + NoLegend()

#Run Azimuth on the Harmony integration
lungquery <- RunAzimuth(lungquery, reference = "D:/Azimuth Reference/")

#Check the Dimensionality Plots
Plot1 <- DimPlot(lungquery, group.by = "predicted.ann_finest_level", raster = FALSE, label = TRUE, repel = TRUE) + NoLegend() 

Plot2 <- DimPlot(lungquery, group.by = "predicted.ann_level_4", raster = FALSE, label = TRUE, repel = TRUE) 

ggsave(Plot1, filename = 'D:/April12-2024/annotated_finestlevel.tiff', width = 10, height =8)

ggsave(Plot2, filename = 'D:/April12-2024/annotationlevel_4.tiff', width = 10, height =8)

# Set Idents
Idents(lungquery)

DimPlot(lungquery)

Idents(lungquery) <- "predicted.ann_finest_level"

DimPlot(lungquery)

# Azimuth normalizes data before mapping, but does not return the results
# normalize the data here before visualization.
lungquery <- NormalizeData(lungquery)

FeaturePlot(lungquery, features = c ("EPCAM", "PTPRC", "CLDN5", "COL1A2"))

# Save RDS
saveRDS(lungquery, "D:/April12-2024/HarmonyIntegrated_annotated_data.rds")

# Find Variable features after renaming
lungquery <- FindVariableFeatures(lungquery, selection.method = "vst", nfeatures = 3000)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones - Takes approximately a minute per cluster depending on the number of clusters.
lungquery <- FindAllMarkers(lungquery, only.pos = TRUE)
lungquery %>%
  group_by(cluster)

write_xlsx(lungquery, "D:/April12-2024/Annotated_after_Azimuth_annotations.xlsx")

clusters_freq <- annotated_data@meta.data %>%
  group_by(orig.ident, predicted.ann_finest_level) %>%
  summarise(n=n()) %>%
  mutate(relative_freq = n/sum(n))
clusters_freq_$predicted.ann_finest_level <- factor(clusters_freq$predicted.ann_finest_level)

clusters_freq

write.csv(as.data.frame (clusters_freq),file="D:/April12-2024/Cluster_frequency_detailed.csv")

clusters_freq_2 <- annotated_data@meta.data %>%
  group_by(sample, predicted.ann_finest_level) %>%
  summarise(n=n()) %>%
  mutate(relative_freq = n/sum(n))
clusters_freq_2$predicted.ann_finest_level <- factor(clusters_freq_2$predicted.ann_finest_level)

clusters_freq_2

write.csv(as.data.frame (clusters_freq_2),file="D:/April12-2024/Cluster_frequency_sample_detailed.csv")
