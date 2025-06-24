library(Seurat)
library(sctransform)
library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(Matrix)

PATH_results = "./output/"

harmonized <- readRDS("../hdrg_proteomics/data/multi-omics/DRG_nonneurons_human.rds")

#-------------------------------------------------------------------------------

Seurat::Assays(harmonized)
Seurat::Idents(harmonized) <- "Atlas_annotation"
# 
# Seurat::DotPlot(harmonized, 
#                 assay ="RNA",
#                 features = c("Raet1a", "Raet1b", "Raet1c", "Raet1d", "Raet1e", "Raet1f", "Raet1g", "Raet1l", 
#                              "Ulbp1", "Ulbp2", "Ulbp3", "Mill2")) + 
#   theme(axis.text.x = element_text(angle = 45, hjust=1)) 
# 
# 
# features = c("KLRC1", "KIT", "KLRD1", "CD56", "CXCR3", "ICAM1", "NCR1")
# 
# FeaturePlot(harmonized, features = features)
g <- DimPlot(harmonized) + theme_bw()

pdf(paste0(PATH_results, "immune_NKs_fullumap.pdf"), height = 4, width = 5)
print(g)
dev.off()


#-------------------------------------------------------------------------------

seu.immune <- subset(harmonized, subset = Atlas_annotation == "Immune")
Seurat::Idents(seu.immune) <- "Atlas_annotation"

rm(harmonized)

# rebuild clean
counts_mat <- GetAssayData(seu.immune, assay = "RNA", layer = "counts")
meta_df    <- seu.immune@meta.data
seu.immune <- CreateSeuratObject(counts = counts_mat, meta.data = meta_df)

# cluster immune subtypes
seu.immune <- NormalizeData(seu.immune)
seu.immune <- FindVariableFeatures(seu.immune, selection.method = "vst", nfeatures = 2000)
seu.immune <- ScaleData(seu.immune)
seu.immune <- RunPCA(seu.immune)
ElbowPlot(seu.immune)  # Check how many PCs to use

seu.immune <- FindNeighbors(seu.immune, dims = 1:10)
seu.immune <- FindClusters(seu.immune, resolution = 1)

seu.immune <- RunUMAP(seu.immune, dims = 1:10)
DimPlot(seu.immune, reduction = "umap", label = TRUE)

# check for non-neuronal clusters after neuron-specific clustering
features <- c(
  "Acap1", "Tmc8", "Cd3e", "Cd8a",                  # T cells
  "Myh11", "Acta2",                                 # SMCs
  "Fabp7", "Apoe", "Serpina5", "Aatk", "Lgi4", "Gnai1",  # SCGs (note: "GIAI" likely meant "GNAI1")
  "Tinagl1", "Notch3",                              # Pericytes
  "Mpz", "Mbp", "Ncam1", "Scn7a", "Mlip", "Prx",    # Schwann
  "Snap25", "Prph", "Tubb2a", "Rbfox3",             # Neurons
  "Plcb2", "Csf1r",                                 # Monocytes
  "Col1a1", "Podn", "Dcn", "Pi16", "Pdgfra",        # Fibroblasts
  "Vwf", "Egfl7", "Cd34", "Tek",                    # Endothelial
  "Plin4", "Plin1", "Cidea", "Adipoq", "Fabp4",     # Adipocytes
  "Cd14", "Cd163", "Ptprc",                         # Immune
  "Tmsb10", "Tmsb4x"                                # Misc / junk
)

# NK features
features <- c("Klrd1", "Nkg7", "Gnly", "Gzmb",   
  "Prf1", "Klrf1", "Fcgr3", "Klrc1", "Klrc2", "Klrc3", "Klrc4", "Il2rb", "Klrk1")

DotPlot(seu.immune, scale = TRUE,
        features = features, 
        assay ="RNA") + theme(axis.text.x = element_text(angle = 45, hjust = 1))

FeaturePlot(
  object = seu.immune,
  features = c("Klrd1", "Klrc1"),
  #blend = TRUE,
  #blend.threshold = 0.01,
  pt.size = 2
)

DimPlot(seu.immune) + theme_bw()


df1 <- FetchData(seu.immune, vars = c("umap_1", "umap_2", "Klrd1"))
df2 <- FetchData(seu.immune, vars = c("umap_1", "umap_2", "Klrk1"))

# Background: all cells in grey
g1 <- ggplot(df1, aes(x = umap_1, y = umap_2)) +
  geom_point(color = "lightgrey", size = 0.5) +
  # Foreground: only cells expressing Klrd1
  geom_point(data = subset(df1, Klrd1 > 0), 
             aes(color = Klrd1), size = 1.5) +
  scale_color_viridis_c(option = "plasma") +
  theme_bw()

# Background: all cells in grey
g2 <- ggplot(df2, aes(x = umap_1, y = umap_2)) +
  geom_point(color = "lightgrey", size = 0.5) +
  # Foreground: only cells expressing Klrc1
  geom_point(data = subset(df2, Klrk1 > 0), 
             aes(color = Klrk1), size = 1.5) +
  scale_color_viridis_c(option = "plasma") +
  theme_bw()

g3 <- DimPlot(seu.immune) + theme_bw()


pdf(paste0(PATH_results, "immune_NKs.pdf"), height = 4, width = 13)
g1 + g2 + g3
dev.off()

pdf(paste0(PATH_results, "immune_NKs_dots.pdf"), height = 4, width = 6)
DotPlot(seu.immune, scale = TRUE,
        features = features, 
        assay ="RNA") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
