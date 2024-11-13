library(ggrepel)
library(ggplot2)
library(dplyr)

PATH_results = "./output/"

sig   <- read.csv("./data/RE_ NK in human sequencing data/DE_analysis_IPSDSN_vs_iPSC.csv", row.names = 1)
names <- read.csv("./data/RE_ NK in human sequencing data/genes_names.csv", row.names = 1)

sig$ensembl_gene_id <- rownames(sig)

input <- merge(sig, names, by = "ensembl_gene_id", all.x = TRUE)
input$refseq_mrna <- NULL

input <- input[!duplicated(input),]

g <- ggplot(input, aes(x = log2FoldChange, y = -log10(pvalue))) 
g <- g + geom_point(color = "#B63679ff", size = 4) 
g <- g + theme(aspect.ratio=1) + theme_bw() +
  ggrepel::geom_text_repel(aes(label=symbol), size=4, segment.alpha= 0.2, force =2, max.overlaps=16) +
  labs(title = "Volcano Plot",
       x = "Log2 Fold Change",
       y = "-log10(pvalue)") +
  theme(plot.title = element_text(hjust = 0.5))

plot(g)

pdf(paste0(PATH_results, "ipsc_volcano.pdf"), height = 4, width = 4)
plot(g)
dev.off()

#-------------------------------------------------------------------------------

mat   <- read.csv("./data/RE_ NK in human sequencing data/genes_eset_alex.csv", row.names = 1)
rownames(mat) <- mat$symbol

col_names <- colnames(mat)

# need to check what the metadata actually is; placeholder here.
metadata <- data.frame(
  sample    = col_names
)

metadata <- metadata %>%
  mutate(group = case_when(
    grepl("^N", sample) ~ "iPSC-SN",            # Names starting with "N"
    grepl("^ON_", sample) ~ "old-SN",      # "ON1" or "ON2"
    grepl("^ONS", sample) ~ "old-SN-serine",    # Names starting with "ONS"
    grepl("^DRG", sample) ~ "DRG",              # Names starting with "DRG"
    grepl("^i", sample) ~ "iPSC",               # Names starting with "i"
    TRUE ~ NA_character_                        # Default case if no match
  ))

head(metadata)

gene_data <- mat %>%
  pivot_longer(cols = !starts_with("symbol"), 
    names_to = "sample", 
    values_to = "expression")

gene_data <- gene_data %>%
  left_join(metadata, by = c("sample" = "sample"))

metadata$group     <- as.factor(metadata$group)

gene_data <- gene_data[gene_data$group %in% c("iPSC", "iPSC-SN"), ] #, "DRG"

ggplot(gene_data, aes(x = symbol, y = expression, fill = group)) +
  geom_boxplot() +
  labs(title = "",
       x = "",
       y = "Expression") + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

pdf(paste0(PATH_results, "ipsc_boxplot.pdf"), height = 5, width = 7)
ggplot(gene_data, aes(x = symbol, y = expression, fill = group)) +
  geom_boxplot() +
  labs(title = "",
       x = "",
       y = "Expression") + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
dev.off()
