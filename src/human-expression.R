library(Seurat)
library(sctransform)
library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(Matrix)

PATH_results = "./output/"

#-------------------------------------------------------------------------------

#Available human data
# hDRG bulk - Wangzhou
# hDRG bulk - ray
# hDRG spatial - diana, lumbar
# hDRG harmonized sn - mixed level
# hDRG sc - Ish. Published???

#-------------------------------------------------------------------------------

load("../hdrg_immune/processing/Ray2023_gsea.RData")
ray_colData <- read.csv("../hdrg_immune/processing/ray2023-colData.csv")
ray_colData$sample_id <- paste0("X", ray_colData$sample_id) #add X (rownames can't start with a number)

head(tcounts)

ensembl   <- biomaRt::useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
gene_info <- biomaRt::getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'), 
                   filters = 'ensembl_gene_id', 
                   values = tcounts$gene_id, 
                   mart = ensembl)

tcounts <- merge(tcounts, gene_info, by.x = "gene_id", by.y = "ensembl_gene_id", all.x = TRUE)


gene_data <- tcounts %>%
  filter(hgnc_symbol %in% c("RAET1A", "RAET1B", "RAET1C", "RAET1D", "RAET1E", "RAET1F", "RAET1G", "RAET1L",
                     "ULBP1", "ULBP2", "ULBP3", "ULBP4", "ULBP5", "ULBP6",
                     "MICA", "MICB"))

# Create the plot with ggplot
gene_data$pain.state <- as.factor(gene_data$pain.state)
gene_data$hgnc_symbol <- as.factor(gene_data$hgnc_symbol)

# Create the main plot with `interaction(state, gene)` on the x-axis
g <- ggplot(gene_data, aes(x = interaction(pain.state, hgnc_symbol), y = expression, fill = pain.state)) +
  geom_boxplot(width = 0.6, position = position_dodge(width = 0.7)) +  # Narrow boxplots and adjust spacing
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.7)) +
  theme_bw() +
  labs(title = "", y = "TPM") +
  theme(
    axis.text.x = element_blank(),       # Hide default x-axis labels
    axis.ticks.x = element_blank(),      # Hide x-axis ticks
    panel.grid.major.x = element_blank(),# Remove vertical grid lines
    plot.margin = margin(0, 0, 30, 0)    # Add margin at the bottom for labels
  )

# Calculate x-axis positions for each gene by finding the midpoint of each group of boxes
gene_positions <- gene_data %>%
  group_by(hgnc_symbol) %>%
  summarize(x_position = mean(as.numeric(interaction(pain.state, hgnc_symbol)))) %>%
  ungroup()

# Add gene labels using `annotate()`, centered below each group of boxes
for (i in seq_len(nrow(gene_positions))) {
  g <- g + annotate(
    "text",
    x = gene_positions$x_position[i],
    y = min(gene_data$expression) - 1, # Adjust y-position as needed
    label = gene_positions$hgnc_symbol[i],
    size = 4, hjust = 0.5, vjust = 1, angle = 45
  )
}

print(g)

pdf(paste0(PATH_results, "ray_boxplot.pdf"), height = 5, width = 7)
print(g)
dev.off()

#-------------------------------------------------------------------------------

harmonized <- readRDS("../hdrg_proteomics/data/multi-omics/DRG_neurons_human.rds")

Seurat::Assays(harmonized)
Seurat::Idents(harmonized) <- "Atlas_annotation"

Seurat::VlnPlot(harmonized, 
                assay ="RNA",
                features = c("Raet1a", "Raet1b", "Raet1c", "Raet1d", "Raet1e", "Raet1f", "Raet1g", "Raet1l", 
                             "Ulbp1", "Ulbp2", "Ulbp3", "Mill2")) + 
                theme(axis.text.x = element_text(angle = 45, hjust=1)) 

Seurat::DotPlot(harmonized, 
                assay ="RNA",
                features = c("Raet1a", "Raet1b", "Raet1c", "Raet1d", "Raet1e", "Raet1f", "Raet1g", "Raet1l", 
                             "Ulbp1", "Ulbp2", "Ulbp3", "Mill2")) + 
                theme(axis.text.x = element_text(angle = 45, hjust=1)) 


pdf(paste0(PATH_results, "harmonized_dots.pdf"), height = 6, width = 5)
Seurat::DotPlot(harmonized, 
                assay ="RNA",
                features = c("Raet1a", "Raet1b", "Raet1c", "Raet1d", "Raet1e", "Raet1f", "Raet1g", "Raet1l", 
                             "Ulbp1", "Ulbp2", "Ulbp3", "Mill2")) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) 
dev.off()

#-------------------------------------------------------------------------------


harmonized <- readRDS("../hdrg_proteomics/data/multi-omics/DRG_nonneurons_human.rds")

Seurat::Assays(harmonized)
Seurat::Idents(harmonized) <- "Atlas_annotation"

Seurat::DotPlot(harmonized, 
                assay ="RNA",
                features = c("Raet1a", "Raet1b", "Raet1c", "Raet1d", "Raet1e", "Raet1f", "Raet1g", "Raet1l", 
                             "Ulbp1", "Ulbp2", "Ulbp3", "Mill2")) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) 


pdf(paste0(PATH_results, "harmonized_non-neuronal_dots.pdf"), height = 6, width = 5)
Seurat::DotPlot(harmonized, 
                assay ="RNA",
                features = c("Raet1a", "Raet1b", "Raet1c", "Raet1d", "Raet1e", "Raet1f", "Raet1g", "Raet1l", 
                             "Ulbp1", "Ulbp2", "Ulbp3", "Mill2")) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) 
dev.off()

#-------------------------------------------------------------------------------



visium.lumbar <- readRDS("../pig_deconvolution/data/drg.combined.human_forCONOS.rds")

Seurat::Assays(visium.lumbar)

Seurat::DotPlot(visium.lumbar, 
                assay ="RNA",
                features = c("RAET1A", "RAET1B", "RAET1C", "RAET1D", "RAET1E", "RAET1F", "RAET1G", "RAET1L",
                             "ULBP1", "ULBP2", "ULBP3", "ULBP4", "ULBP5", "ULBP6",
                             "MICA", "MICB")) + theme(axis.text.x = element_text(angle = 45, hjust=1)) 


pdf(paste0(PATH_results, "visium.lumbar_dots.pdf"), height = 5, width = 6)
Seurat::DotPlot(visium.lumbar, 
                assay ="RNA",
                features = c("RAET1A", "RAET1B", "RAET1C", "RAET1D", "RAET1E", "RAET1F", "RAET1G", "RAET1L",
                             "ULBP1", "ULBP2", "ULBP3", "ULBP4", "ULBP5", "ULBP6",
                             "MICA", "MICB")) + theme(axis.text.x = element_text(angle = 45, hjust=1)) 
dev.off()


#-------------------------------------------------------------------------------

c2.sn <- readRDS("../hdrg_c2study/data/merged_FINAL_clean_neuronsubset_102524.rds") #not for publishing

Seurat::Assays(c2.sn)

Seurat::DotPlot(c2.sn, 
                assay ="RNA",
                group.by = "Condition",
                features = c("RAET1A", "RAET1B", "RAET1C", "RAET1D", "RAET1E", "RAET1F", "RAET1G", "RAET1L",
                             "ULBP1", "ULBP2", "ULBP3", "ULBP4", "ULBP5", "ULBP6",
                             "MICA", "MICB")) + theme(axis.text.x = element_text(angle = 45, hjust=1)) 


pdf(paste0(PATH_results, "c2.sn_dots.pdf"), height = 5, width = 6)
Seurat::DotPlot(c2.sn, 
                assay ="RNA",
                group.by = "Condition",
                features = c("RAET1A", "RAET1B", "RAET1C", "RAET1D", "RAET1E", "RAET1F", "RAET1G", "RAET1L",
                             "ULBP1", "ULBP2", "ULBP3", "ULBP4", "ULBP5", "ULBP6",
                             "MICA", "MICB")) + theme(axis.text.x = element_text(angle = 45, hjust=1)) 
dev.off()

#-------------------------------------------------------------------------------

wangzhou      <- read.csv("../hdrg_proteomics/data/multi-omics/wangzhou2020_tpm.csv")
wangzhou_meta <- read.csv("../hdrg_proteomics/data/multi-omics/wangzhou2020_metadata.csv")

gene_data <- wangzhou %>%
  filter(gene %in% c("RAET1A", "RAET1B", "RAET1C", "RAET1D", "RAET1E", "RAET1F", "RAET1G", "RAET1L",
                     "ULBP1", "ULBP2", "ULBP3", "ULBP4", "ULBP5", "ULBP6",
                     "MICA", "MICB", "P2RX3"))

gene_data <- gene_data %>%
  pivot_longer(cols = starts_with("hDIV") | starts_with("hDRG"), 
               names_to = "sample", 
               values_to = "expression")

gene_data$sample <- gsub("\\.", "-", gene_data$sample)

gene_data <- gene_data %>%
  left_join(wangzhou_meta, by = c("sample" = "sample"))

# Create the plot with ggplot
gene_data$state <- as.factor(gene_data$state)
gene_data$gene <- as.factor(gene_data$gene)

# Create the main plot with `interaction(state, gene)` on the x-axis
g <- ggplot(gene_data, aes(x = interaction(state, gene), y = expression, fill = state)) +
  geom_boxplot(width = 0.6, position = position_dodge(width = 0.7)) +  # Narrow boxplots and adjust spacing
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.7)) +
  theme_bw() +
  labs(title = "", y = "TPM") +
  theme(
    axis.text.x = element_blank(),       # Hide default x-axis labels
    axis.ticks.x = element_blank(),      # Hide x-axis ticks
    panel.grid.major.x = element_blank(),# Remove vertical grid lines
    plot.margin = margin(0, 0, 30, 0)    # Add margin at the bottom for labels
  )

# Calculate x-axis positions for each gene by finding the midpoint of each group of boxes
gene_positions <- gene_data %>%
  group_by(gene) %>%
  summarize(x_position = mean(as.numeric(interaction(state, gene)))) %>%
  ungroup()

# Add gene labels using `annotate()`, centered below each group of boxes
for (i in seq_len(nrow(gene_positions))) {
  g <- g + annotate(
    "text",
    x = gene_positions$x_position[i],
    y = min(gene_data$expression) - 1, # Adjust y-position as needed
    label = gene_positions$gene[i],
    size = 4, hjust = 0.5, vjust = 1, angle = 45
  )
}

print(g)

pdf(paste0(PATH_results, "wangzhou_boxplot.pdf"), height = 5, width = 7)
print(g)
dev.off()

#-------------------------------------------------------------------------------

library(data.table)
library(sctransform)
library(Seurat)

lu_data <- fread("./data/GSE249746_Expression_matrix_raw_counts.csv.gz")

lu_data <- as.data.frame(lu_data)
rownames(lu_data) <- lu_data$V1
lu_data$V1 <- NULL

# ss.data <- CreateSeuratObject(counts = lu_data, project = "lu_singlesomma")
# ss.data[["percent.mt"]] <- PercentageFeatureSet(ss.data, pattern = "^MT-")
# VlnPlot(ss.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# 
# ss.data <- SCTransform(ss.data, vars.to.regress = "percent.mt", verbose = FALSE)
# 
# ss.data <- RunPCA(ss.data, verbose = FALSE)
# ss.data <- RunUMAP(ss.data, dims = 1:30, verbose = FALSE)
# 
# ss.data <- FindNeighbors(ss.data, dims = 1:30, verbose = FALSE)
# ss.data <- FindClusters(ss.data, verbose = FALSE)
# ss.data <- DimPlot(ss.data, label = TRUE)

donor_ids <- unique(sub("^(\\w+-\\w+)-.*", "\\1", colnames(lu_data)))

# Loop through each donor ID to create a Seurat object
seurat_objects <- list()
for (donor in donor_ids) {
  donor_data <- lu_data[, grepl(paste0("^", donor, "-"), colnames(lu_data))]
  seurat_obj <- CreateSeuratObject(counts = donor_data)

  seurat_obj <- NormalizeData(seurat_obj)
  #seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 4500)

  seurat_objects[[donor]] <- seurat_obj
}

ss.data <- FindIntegrationAnchors(object.list = seurat_objects,
  scale=TRUE, 
  normalization.method="LogNormalize",
  reduction="cca", 
  l2.norm = TRUE,
  dims=1:30, 
  k.anchor=5, 
  k.filter=200, 
  k.score=30, 
  max.features=200,
  nn.method="annoy", 
  n.trees=50, 
  eps=0)

ss.data <- IntegrateData(ss.data,
  normalization.method="LogNormalize", 
  features=NULL, 
  features.to.integrate=NULL, 
  dims=1:30, 
  k.weight=80, 
  weight.reduction=NULL, 
  sd.weight=1, 
  sample.tree=NULL, 
  preserve.order=FALSE, 
  eps=0) 

ss.data <- ScaleData(ss.data)
ss.data <- RunPCA(ss.data, npcs=50)

ss.data <- RunUMAP(ss.data, reduction="pca", dims=1:25)
ss.data <- FindNeighbors(ss.data, reduction="pca", dims=1:25) 
ss.data <- FindClusters(ss.data, resolution=3.4)

saveRDS(ss.data, "./data/lu-singlesoma.RDS")

DimPlot(ss.data, reduction = "umap")

table(Idents(ss.data))

# Highly similar clusters without clearly distinguishable markers
# were merged to produce the final 16 clusters.

load("./data/human_meta_final_cluster.Rdata")

ss.data@meta.data <- human_meta_final_cluster
ss.data <- SetIdent(ss.data, value = ss.data@meta.data$cl.conserv_final)

library(ggplot2)

# Seurat::FeaturePlot(ss.data,
#                 features = c("TRPM8", "TRPV1", "CHRNA7", "SST", 
#                              "MRGPRX1", "PIEZO2", "KIT", "NTRK2",
#                              "TRPA1", "NTRK3", "PVALB")) + theme(axis.text.x = element_text(angle = 45, hjust=1)) 

Seurat::DotPlot(ss.data,
                assay ="RNA",
                features = c("TRPM8", "TRPV1", "CHRNA7", "SST", 
                             "MRGPRX1", "PIEZO2", "KIT", "NTRK2",
                             "TRPA1", "NTRK3", "PVALB")) + theme(axis.text.x = element_text(angle = 45, hjust=1)) 


PATH_results = "./output/"

pdf(paste0(PATH_results, "ss_dots_named.pdf"), height = 5, width = 6)
Seurat::DotPlot(ss.data, 
                assay ="RNA",
                features = c("RAET1A", "RAET1B", "RAET1C", "RAET1D", "RAET1E", "RAET1F", "RAET1G", "RAET1L",
                             "ULBP1", "ULBP2", "ULBP3", "ULBP4", "ULBP5", "ULBP6",
                             "MICA", "MICB")) + theme(axis.text.x = element_text(angle = 45, hjust=1)) 
dev.off()





