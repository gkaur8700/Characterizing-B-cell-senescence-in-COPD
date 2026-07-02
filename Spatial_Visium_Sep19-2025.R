# Sep 18 2025

# Visium Final Analyses

#Load libraries
install.packages(c("SeuratObject", "Seurat")) # To install/update both

library(Seurat)
library(dplyr)
library(ggplot2)


#Load the samples
ns1 <- Load10X_Spatial(data.dir = "/Volumes/My Passport/deliv_Rahman_Spatial_121321_counts/Sample_1_346/outs", slice = "ns1")
ns2 <- Load10X_Spatial(data.dir = "/Volumes/My Passport/deliv_Rahman_Spatial_121321_counts/Sample_1_350/outs", slice = "ns2")
ns3 <- Load10X_Spatial(data.dir = "/Volumes/My Passport/deliv_Rahman_Spatial_121321_counts/Sample_2_350/outs", slice = "ns3")

sm1 <- Load10X_Spatial(data.dir = "/Volumes/My Passport/deliv_Rahman_Spatial_121321_counts/Sample_2_346/outs", slice = "sm1")

c1 <- Load10X_Spatial(data.dir = "/Volumes/My Passport/deliv_Rahman_Spatial_121321_counts/Sample_3_346/outs", slice = "c1")
c2 <- Load10X_Spatial(data.dir = "/Volumes/My Passport/deliv_Rahman_Spatial_121321_counts/Sample_3_350/outs", slice = "c2")
c3 <- Load10X_Spatial(data.dir = "/Volumes/My Passport/deliv_Rahman_Spatial_121321_counts/Sample_4_350/outs", slice = "c3")

#Data Preprocessing

#log normalization
ns1 <- NormalizeData(ns1)
ns2 <- NormalizeData(ns2)
ns3 <- NormalizeData(ns3)

sm1 <- NormalizeData(sm1)

c1 <- NormalizeData(c1)
c2 <- NormalizeData(c2)
c3 <- NormalizeData(c3)

# get the valcano plot first
plot1 <- VlnPlot(ns1, features = "nCount_Spatial", pt.size = 0.1) 
plot2 <- VlnPlot(ns2, features = "nCount_Spatial", pt.size = 0.1) 
plot3 <- VlnPlot(ns3, features = "nCount_Spatial", pt.size = 0.1) 

plot1|plot2|plot3

plot4 <- VlnPlot(c1, features = "nCount_Spatial", pt.size = 0.1)
plot5 <- VlnPlot(c2, features = "nCount_Spatial", pt.size = 0.1) 
plot6 <- VlnPlot(c3, features = "nCount_Spatial", pt.size = 0.1) 

plot4|plot5|plot6

plot7 <- VlnPlot(sm1, features = "nCount_Spatial", pt.size = 0.1)

plot7

ggsave(plot1, filename = '~/Desktop/Spatial_visium/QC_ns1.tiff', width = 10, height =8)
ggsave(plot2, filename = '~/Desktop/Spatial_visium/QC_ns2.tiff', width = 10, height =8)
ggsave(plot3, filename = '~/Desktop/Spatial_visium/QC_ns3.tiff', width = 10, height =8)

ggsave(plot4, filename = '~/Desktop/Spatial_visium/QC_c1.tiff', width = 10, height =8)
ggsave(plot5, filename = '~/Desktop/Spatial_visium/QC_c2.tiff', width = 10, height =8)
ggsave(plot6, filename = '~/Desktop/Spatial_visium/QC_c3.tiff', width = 10, height =8)

ggsave(plot7, filename = '~/Desktop/Spatial_visium/QC_sm1.tiff', width = 10, height =8)

# Spatial Feature Plot

plot8 <- SpatialFeaturePlot(ns1, features = "nCount_Spatial") + theme(legend.position = "right")
plot9 <- SpatialFeaturePlot(ns2, features = "nCount_Spatial") + theme(legend.position = "right")
plot10 <- SpatialFeaturePlot(ns3, features = "nCount_Spatial") + theme(legend.position = "right")

plot11 <- SpatialFeaturePlot(c1, features = "nCount_Spatial") + theme(legend.position = "right")
plot12 <- SpatialFeaturePlot(c2, features = "nCount_Spatial") + theme(legend.position = "right")
plot13 <- SpatialFeaturePlot(c3, features = "nCount_Spatial") + theme(legend.position = "right")

plot14 <- SpatialFeaturePlot(sm1, features = "nCount_Spatial") + theme(legend.position = "right")

#save the plots
ggsave(plot8, filename = '~/Desktop/Spatial_visium/SP_ns1.tiff', width = 10, height =8)
ggsave(plot9, filename = '~/Desktop/Spatial_visium/SP_ns2.tiff', width = 10, height =8)
ggsave(plot10, filename = '~/Desktop/Spatial_visium/SP_ns3.tiff', width = 10, height =8)

ggsave(plot11, filename = '~/Desktop/Spatial_visium/SP_c1.tiff', width = 10, height =8)
ggsave(plot12, filename = '~/Desktop/Spatial_visium/SP_c2.tiff', width = 10, height =8)
ggsave(plot13, filename = '~/Desktop/Spatial_visium/SP_c3.tiff', width = 10, height =8)

ggsave(plot14, filename = '~/Desktop/Spatial_visium/SP_sm1.tiff', width = 10, height =8)

#Doing scTranform and compare the two

ns1 <- SCTransform (ns1, assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE)
ns2 <- SCTransform (ns2, assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE)
ns3 <- SCTransform (ns3, assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE)

sm1 <- SCTransform (sm1, assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE)

c1 <- SCTransform (c1, assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE)
c2 <- SCTransform (c2, assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE)
c3 <- SCTransform (c3, assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE)

# Correlate scTransform with log normalization and plot

ns1 <- GroupCorrelation(ns1, group.assay = "Spatial", assay = "Spatial", slot = "data", do.plot = FALSE)
ns1 <- GroupCorrelation(ns1, group.assay = "Spatial", assay = "SCT", slot = "scale.data", do.plot = FALSE)

p1 <- GroupCorrelationPlot(ns1, assay = "Spatial", cor = "nCount_Spatial_cor") + ggtitle("Log Normalization") +
  theme(plot.title = element_text(hjust = 0.5))
p2 <- GroupCorrelationPlot(ns1, assay = "SCT", cor = "nCount_Spatial_cor") + ggtitle("SCTransform Normalization") +
  theme(plot.title = element_text(hjust = 0.5))
plot15 <- p1 + p2

ns2 <- GroupCorrelation(ns2, group.assay = "Spatial", assay = "Spatial", slot = "data", do.plot = FALSE)
ns2 <- GroupCorrelation(ns2, group.assay = "Spatial", assay = "SCT", slot = "scale.data", do.plot = FALSE)

p3 <- GroupCorrelationPlot(ns2, assay = "Spatial", cor = "nCount_Spatial_cor") + ggtitle("Log Normalization") +
  theme(plot.title = element_text(hjust = 0.5))
p4 <- GroupCorrelationPlot(ns2, assay = "SCT", cor = "nCount_Spatial_cor") + ggtitle("SCTransform Normalization") +
  theme(plot.title = element_text(hjust = 0.5))
plot16 <- p3 + p4

ns3 <- GroupCorrelation(ns3, group.assay = "Spatial", assay = "Spatial", slot = "data", do.plot = FALSE)
ns3 <- GroupCorrelation(ns3, group.assay = "Spatial", assay = "SCT", slot = "scale.data", do.plot = FALSE)

p5 <- GroupCorrelationPlot(ns3, assay = "Spatial", cor = "nCount_Spatial_cor") + ggtitle("Log Normalization") +
  theme(plot.title = element_text(hjust = 0.5))
p6 <- GroupCorrelationPlot(ns3, assay = "SCT", cor = "nCount_Spatial_cor") + ggtitle("SCTransform Normalization") +
  theme(plot.title = element_text(hjust = 0.5))
plot17 <- p5 + p6

sm1 <- GroupCorrelation(sm1, group.assay = "Spatial", assay = "Spatial", slot = "data", do.plot = FALSE)
sm1 <- GroupCorrelation(sm1, group.assay = "Spatial", assay = "SCT", slot = "scale.data", do.plot = FALSE)

p7 <- GroupCorrelationPlot(sm1, assay = "Spatial", cor = "nCount_Spatial_cor") + ggtitle("Log Normalization") +
  theme(plot.title = element_text(hjust = 0.5))
p8 <- GroupCorrelationPlot(sm1, assay = "SCT", cor = "nCount_Spatial_cor") + ggtitle("SCTransform Normalization") +
  theme(plot.title = element_text(hjust = 0.5))
plot18 <- p7 + p8

c1 <- GroupCorrelation(c1, group.assay = "Spatial", assay = "Spatial", slot = "data", do.plot = FALSE)
c1 <- GroupCorrelation(c1, group.assay = "Spatial", assay = "SCT", slot = "scale.data", do.plot = FALSE)

p9 <- GroupCorrelationPlot(c1, assay = "Spatial", cor = "nCount_Spatial_cor") + ggtitle("Log Normalization") +
  theme(plot.title = element_text(hjust = 0.5))
p10 <- GroupCorrelationPlot(c1, assay = "SCT", cor = "nCount_Spatial_cor") + ggtitle("SCTransform Normalization") +
  theme(plot.title = element_text(hjust = 0.5))
plot19 <- p9 + p10

c2 <- GroupCorrelation(c2, group.assay = "Spatial", assay = "Spatial", slot = "data", do.plot = FALSE)
c2 <- GroupCorrelation(c2, group.assay = "Spatial", assay = "SCT", slot = "scale.data", do.plot = FALSE)

p11 <- GroupCorrelationPlot(c2, assay = "Spatial", cor = "nCount_Spatial_cor") + ggtitle("Log Normalization") +
  theme(plot.title = element_text(hjust = 0.5))
p12 <- GroupCorrelationPlot(c2, assay = "SCT", cor = "nCount_Spatial_cor") + ggtitle("SCTransform Normalization") +
  theme(plot.title = element_text(hjust = 0.5))
plot20 <- p11 + p12

c3 <- GroupCorrelation(c3, group.assay = "Spatial", assay = "Spatial", slot = "data", do.plot = FALSE)
c3 <- GroupCorrelation(c3, group.assay = "Spatial", assay = "SCT", slot = "scale.data", do.plot = FALSE)

p13 <- GroupCorrelationPlot(c3, assay = "Spatial", cor = "nCount_Spatial_cor") + ggtitle("Log Normalization") +
  theme(plot.title = element_text(hjust = 0.5))
p14 <- GroupCorrelationPlot(c3, assay = "SCT", cor = "nCount_Spatial_cor") + ggtitle("SCTransform Normalization") +
  theme(plot.title = element_text(hjust = 0.5))
plot21 <- p13 + p14

#save the plots
ggsave(plot15, filename = '~/Desktop/Spatial_visium/nor_ns1.tiff', width = 10, height =8)
ggsave(plot16, filename = '~/Desktop/Spatial_visium/nor_ns2.tiff', width = 10, height =8)
ggsave(plot17, filename = '~/Desktop/Spatial_visium/nor_ns3.tiff', width = 10, height =8)

ggsave(plot19, filename = '~/Desktop/Spatial_visium/nor_c1.tiff', width = 10, height =8)
ggsave(plot20, filename = '~/Desktop/Spatial_visium/nor_c2.tiff', width = 10, height =8)
ggsave(plot21, filename = '~/Desktop/Spatial_visium/nor_c3.tiff', width = 10, height =8)

ggsave(plot18, filename = '~/Desktop/Spatial_visium/nor_sm1.tiff', width = 10, height =8)

#merge ns together

ns.merge <- merge(ns1, y = c(ns2, ns3), add.cell.ids = c("nonsmoker1", "nonsmoker2", "nonsmoker3"))
#add a column with sample name
ns.merge$sample <- sapply(strsplit(colnames(ns.merge), "_"), `[`, 1)

c.merge <- merge(c1, y = c(c2, c3), add.cell.ids = c("copd1", "copd2", "copd3"))
c.merge$sample <- sapply(strsplit(colnames(c.merge), "_"), `[`, 1)

#make scTransform as the default assay

DefaultAssay(ns.merge)
DefaultAssay(c.merge)
DefaultAssay(sm1)

VariableFeatures(ns.merge) <- c(VariableFeatures(ns1), VariableFeatures(ns2), VariableFeatures(ns3))
ns.merge <- RunPCA(ns.merge, verbose = FALSE)
ns.merge <- FindNeighbors(ns.merge, dims = 1:30)
ns.merge <- FindClusters(ns.merge, verbose = FALSE)
ns.merge <- RunUMAP(ns.merge, dims = 1:30)

VariableFeatures(c.merge) <- c(VariableFeatures(c1), VariableFeatures(c2), VariableFeatures(c3))
c.merge <- RunPCA(c.merge, verbose = FALSE)
c.merge <- FindNeighbors(c.merge, dims = 1:30)
c.merge <- FindClusters(c.merge, verbose = FALSE)
c.merge <- RunUMAP(c.merge, dims = 1:30)

VariableFeatures(sm1) <- c(VariableFeatures(sm1))
sm1 <- RunPCA(sm1, verbose = FALSE)
sm1 <- FindNeighbors(sm1, dims = 1:30)
sm1 <- FindClusters(sm1, verbose = FALSE)
sm1 <- RunUMAP(sm1, dims = 1:30)

#Visualize the merged files

plot22 <- DimPlot(ns.merge, reduction = "umap", group.by = c("ident", "sample"))

ggsave(plot22, filename = '~/Desktop/Spatial_visium/Merged_nonsmoker.tiff', width = 10, height =8)

plot23 <- DimPlot(c.merge, reduction = "umap", group.by = c("ident", "sample"))

ggsave(plot23, filename = '~/Desktop/Spatial_visium/Merged_copd.tiff', width = 10, height =8)

plot24 <- DimPlot(sm1, reduction = "umap", group.by = c("ident"))

ggsave(plot24, filename = '~/Desktop/Spatial_visium/Smoker.tiff', width = 10, height =8)

#Spatial plots without integration

plot25 <- SpatialDimPlot(ns.merge)
plot26 <- SpatialDimPlot(c.merge)
plot27 <- SpatialDimPlot(sm1)

ggsave(plot25, filename = '~/Desktop/Spatial_visium/SpatialDimplot_nonsmoker.tiff', width = 10, height =8)
ggsave(plot26, filename = '~/Desktop/Spatial_visium/SpatialDimplot_copd.tiff', width = 10, height =8)
ggsave(plot27, filename = '~/Desktop/Spatial_visium/SpatialDimplot_smoker.tiff', width = 10, height =8)

#Integrate the datasets

#first normalise the merged dataset
ns.merge <- SCTransform(ns.merge, assay = "Spatial", ncells = 3000, verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)

plot28 <- DimPlot(ns.merge, reduction = "umap", group.by = c("ident", "sample"))
ggsave(plot28, filename = '~/Desktop/Spatial_visium/SpatialDimplot_nonsmoker_afterintegration.tiff', width = 10, height =8)

c.merge <- SCTransform(c.merge, assay = "Spatial", ncells = 3000, verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)

plot29 <- DimPlot(c.merge, reduction = "umap", group.by = c("ident", "sample"))
ggsave(plot29, filename = '~/Desktop/Spatial_visium/SpatialDimplot_copd_afterintegration.tiff', width = 10, height =8)

all.merge <- merge(ns.merge, y = c(c.merge, sm1), add.cell.ids = c("nonsmoker", "copd", "smoker"))

all.merge$group <- sapply(strsplit(colnames(all.merge), "_"), `[`, 1)

VariableFeatures(all.merge) <- c(VariableFeatures(ns.merge), VariableFeatures(c.merge), VariableFeatures(sm1))
all.merge <- RunPCA(all.merge, verbose = FALSE)
all.merge <- FindNeighbors(all.merge, dims = 1:30)
all.merge <- FindClusters(all.merge, verbose = FALSE)
all.merge <- RunUMAP(all.merge, dims = 1:30)

#Visualize the merged NS files

plot30 <- DimPlot(all.merge, reduction = "umap", group.by = c("ident", "group"))

ggsave(plot30, filename = '~/Desktop/Spatial_visium/All_merged.tiff', width = 10, height =8)

saveRDS(all.merge, "~/Desktop/Spatial_visium/All_merged_after_integration.rds") 
