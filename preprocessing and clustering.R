library(Seurat) 
library(ggplot2)
library(tidyverse)
library(writexl)
library(DESeq2)

data <- readRDS("/Volumes/My Passport/Final Data/HarmonyIntegrated_annotated_data_final.rds")

data$CellTypes <- Idents(data)

p1 <- DimPlot(data,
              reduction = "umap",
              label = TRUE,
              label.size = 5, fontface = "bold", pt.size = 0.5, repel = TRUE
) + NoLegend()

# Make the fonts Bold

p2 <- p1 + theme(
  text = element_text(face = "bold"),
  legend.text = element_text(face = "bold"),
  axis.title = element_text(face = "bold"),
  axis.text = element_text(face = "bold")
)

# Plot the UMAP of before and after to compare
Plot1 <- DimPlot(data, split.by = "orig.ident",
                reduction = "umap",
                label = FALSE,
                label.size = 5, pt.size = 0.5, repel = TRUE
) + NoLegend()

#Save the Plot

ggsave(Plot1, filename = '~/Desktop/Group_split.tiff', width = 10, height =8)

Plot2 <- DimPlot(data,
                reduction = "umap",
                label = TRUE,
                label.size = 5, pt.size = 0.5, repel = TRUE
)

Plot3 <- DimPlot(data,
                reduction = "umap", pt.size = 0.8,
                label = FALSE,
)  + NoLegend()

ggsave(p2, filename = '~/Desktop/Annotated_umap_Plot_general_final.tiff', width = 10, height =8)
ggsave(Plot2, filename = '~/Desktop/AnnotatedPlotwithlabels_final.tiff', width = 12, height =10)
ggsave(Plot3, filename = '~/Desktop/AnnotatedPlotwithoutlabels_final.tiff', width = 12, height =10)

#Marker file
markers <- readRDS("/Volumes/My Passport/Final Data/annotated_markers.rds")

write.csv(as.data.frame (markers),file="~/Desktop/Marker_genes.csv")
#Select clusters 

select_marker <- read.csv("~/Desktop/select marker.csv")

selected_clusters <- c("AT1", "AT2", "Alveolar macrophages", "Multiciliated (non-nasal)", "B cells", "CD4 T cells", "CD8 T cells", "Classical monocytes", "Non-classical monocytes","DC2", "Mast cells", "NK cells", "Plasma cells")
marker_genes <- c("MS4A4E", "APOBEC3A", "SLC24A4", "ADGRE2", "SLC8A1", "KYNU", "FCAR", "FAM49A", "CHST15", "FCN1", "IGKV1-39", "IGHV1-24", 
"IGHV3-48", "IGHV4-59", "IGKV1-5", "IGHG1", "IGHV3-23", "IGKV3-15", "IGKV4-1", "IGHG3", "CD8B", "CCL5", "CD3G", "TRGC2", "GZMH", "EOMES", 
"CD8A", "KLRC4", "CD3D", "SLA2", "S100A12", "S100A8", "S100A9", "FCN1", "ASGR2", "RNASE2", "IL1R2", "VCAN", "EREG", "FCAR", "TCL1A", "MS4A1", "VPREB3",
"FCRLA", "PAX5", "BANK1", "IGLC6", "BLK", "CR2", "FCRL1", "LGI3", "SFTPC", "PLA2G3", "WFDC12", "SFTPA2", "SFTPA1", "PLA2G1B", "PGC", "NAPSA", "FGG", "LGALS9C", "LAIR2", "FGFBP2", "LGALS9B", "LIM2",
"TRDC", "CX3CR1", "NKG7", "GZMA", "GZMB", "CD40LG", "CCR4", "CTLA4", "THEMIS", "LTB", "CD28", "ICOS", "CAMK4", "LEF1", "NELL2", "IRX4", "HOXA11", "CD1E", "C1QL1", "DLX1", "FCER1A", "SIX2", "CA9", "TFAP2A", "L1CAM",
"AGER", "ANKRD1", "RTKN2", "CEACAM5", "BDNF", "NTM", "UPK3B", "SERTM1", "CEACAM6", "ZBED2", "TPSB2", "CPA3", "TPSAB1", "TPSG1", "MS4A2", "SIGLEC6", "RD3", "CTSG", "TPSD1", "PTGDR2", "DNASE2B", "FABP3",
"RBP4", "PARAL1", "GPD1", "INHBA", "AGRP", "APOC2", "LSAMP", "OR6N1")


subset_obj <- subset(data, idents = selected_clusters)


# For each cluster, get the top 5 variable genes based on the average log fold change
top5_genes_per_cluster <- markers %>%
  group_by(cluster) %>%
  top_n(5, avg_log2FC) %>%
  pull(gene)

write.csv(as.data.frame (top10_genes_per_cluster),file="~/Desktop/Top_10_gene_per_cluster.csv")

top10_genes_expression <- FetchData(data, vars = unique(top10_genes_per_cluster))

DoHeatmap(subset_obj, features = marker_genes) +
  scale_fill_gradientn(colors = c("navy", "white", "firebrick3"))

write.csv(as.data.frame (top10_genes_expression),file="~/Desktop/Top_10_gene_expression.csv")
