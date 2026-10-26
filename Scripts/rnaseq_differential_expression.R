library(DESeq2)
library(GEOquery)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggplot2)
library(pheatmap)
library(clusterProfiler)
library(enrichplot)

project_dir <- file.path(
  "/Users/nathanperkins/Documents/Programming",
  "RNAseq-Differential-Expression"
)

setwd(project_dir)

dir.create("Results", showWarnings = FALSE)
dir.create("Figures", showWarnings = FALSE)

raw_data <- read.csv(
  "Data/Raw/GSE91061_BMS038109Sample.hg19KnownGene.raw.csv.gz",
  check.names = FALSE
)

rownames(raw_data) <- raw_data[[1]]
raw_data <- raw_data[, -1, drop = FALSE]

counts_raw <- as.matrix(raw_data)
storage.mode(counts_raw) <- "numeric"

geo_data <- getGEO("GSE91061", GSEMatrix = TRUE)
geo_metadata <- pData(geo_data[[1]])

geo_match <- match(colnames(counts_raw), geo_metadata$title)
stopifnot(!anyNA(geo_match))

analysis_metadata <- data.frame(
  sample_id = colnames(counts_raw),
  response = geo_metadata[geo_match, "response:ch1"],
  timepoint = geo_metadata[geo_match, "visit (pre or on treatment):ch1"],
  row.names = colnames(counts_raw)
)

keep_samples <- analysis_metadata$timepoint == "Pre" &
  analysis_metadata$response %in% c("PD", "PRCR")

pre_metadata <- analysis_metadata[keep_samples, , drop = FALSE]
pre_counts <- counts_raw[, keep_samples, drop = FALSE]

keep_genes <- rowSums(pre_counts >= 10) >= 3
filtered_counts <- pre_counts[keep_genes, , drop = FALSE]

pre_metadata$response <- factor(
  pre_metadata$response,
  levels = c("PD", "PRCR")
)

dds <- DESeqDataSetFromMatrix(
  countData = filtered_counts,
  colData = pre_metadata,
  design = ~ response
)

dds <- DESeq(dds)

res <- results(
  dds,
  contrast = c("response", "PRCR", "PD"),
  alpha = 0.05
)

res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)

sig_genes <- res_df[
  !is.na(res_df$padj) &
    res_df$padj < 0.05 &
    abs(res_df$log2FoldChange) >= 1,
]

up_in_prcr <- sig_genes[sig_genes$log2FoldChange > 0, ]
up_in_pd <- sig_genes[sig_genes$log2FoldChange < 0, ]

write.csv(
  res_df,
  "Results/all_differential_expression_results.csv",
  row.names = FALSE
)

write.csv(
  sig_genes,
  "Results/significant_genes.csv",
  row.names = FALSE
)

write.csv(
  up_in_prcr,
  "Results/genes_higher_in_PRCR.csv",
  row.names = FALSE
)

write.csv(
  up_in_pd,
  "Results/genes_higher_in_PD.csv",
  row.names = FALSE
)

gene_annotation <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(sig_genes$gene_id),
  keytype = "ENTREZID",
  columns = c("SYMBOL", "GENENAME")
)

sig_genes_annotated <- merge(
  sig_genes,
  gene_annotation,
  by.x = "gene_id",
  by.y = "ENTREZID",
  all.x = TRUE
)

write.csv(
  sig_genes_annotated,
  "Results/significant_genes_annotated.csv",
  row.names = FALSE
)

library_size_data <- data.frame(
  sample_id = colnames(pre_counts),
  response = pre_metadata$response,
  library_size = colSums(pre_counts)
)

library_size_plot <- ggplot(
  library_size_data,
  aes(x = response, y = library_size, fill = response)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.65) +
  geom_jitter(width = 0.12, alpha = 0.75, size = 2) +
  labs(
    title = "Library Sizes by Response Group",
    x = "Response",
    y = "Total Counts"
  ) +
  guides(fill = "none") +
  theme_classic()

ggsave(
  "Figures/library_sizes_by_response.png",
  plot = library_size_plot,
  width = 8,
  height = 6,
  dpi = 300
)


volcano_data <- res_df[
  !is.na(res_df$padj) & !is.na(res_df$log2FoldChange),
]

volcano_data$significant <- volcano_data$padj < 0.05 &
  abs(volcano_data$log2FoldChange) >= 1

volcano_plot <- ggplot(
  volcano_data,
  aes(x = log2FoldChange, y = -log10(padj), color = significant)
) +
  geom_point(alpha = 0.7, size = 1.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  scale_color_manual(
    values = c("FALSE" = "grey70", "TRUE" = "firebrick"),
    labels = c("FALSE" = "Not significant", "TRUE" = "Significant")
  ) +
  labs(
    title = "Differential Expression: PRCR vs PD",
    subtitle = "Positive values indicate higher expression in PRCR; negative values indicate higher expression in PD",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = NULL
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 9),
    legend.position = "right"
  )

ggsave(
  "Figures/volcano_plot.png",
  plot = volcano_plot,
  width = 10,
  height = 6,
  dpi = 300
)


png(
  "Figures/ma_plot.png",
  width = 2400,
  height = 1800,
  res = 300
)

plotMA(
  res,
  alpha = 0.05,
  main = "MA Plot: PRCR vs PD"
)

dev.off()

vsd <- vst(dds, blind = FALSE)

pca_data <- plotPCA(
  vsd,
  intgroup = "response",
  returnData = TRUE
)

percent_variance <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(
  pca_data,
  aes(x = PC1, y = PC2, color = response)
) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    title = "PCA of Pretreatment Melanoma Samples",
    x = paste0("PC1: ", percent_variance[1], "% variance"),
    y = paste0("PC2: ", percent_variance[2], "% variance"),
    color = "Response"
  ) +
  theme_classic()

ggsave(
  "Figures/pca_response.png",
  plot = pca_plot,
  width = 8,
  height = 6,
  dpi = 300
)


gene_variances <- apply(assay(vsd), 1, var)

top_variable_genes <- names(
  sort(gene_variances, decreasing = TRUE)
)[1:50]

heatmap_matrix <- assay(vsd)[top_variable_genes, , drop = FALSE]
heatmap_matrix_scaled <- t(scale(t(heatmap_matrix)))

sample_annotation <- data.frame(
  Response = pre_metadata$response,
  row.names = rownames(pre_metadata)
)

sample_annotation <- sample_annotation[
  colnames(heatmap_matrix_scaled),
  ,
  drop = FALSE
]

annotation_colors <- list(
  Response = c(
    PD = "#F8766D",
    PRCR = "#00BFC4"
  )
)

png(
  "Figures/top_variable_genes_heatmap.png",
  width = 3600,
  height = 2400,
  res = 300
)

pheatmap(
  heatmap_matrix_scaled,
  annotation_col = sample_annotation,
  annotation_colors = annotation_colors,
  show_rownames = FALSE,
  main = "Top 50 Most Variable Genes",
  fontsize_col = 8
)

dev.off()


pd_genes <- unique(as.character(up_in_pd$gene_id))
prcr_genes <- unique(as.character(up_in_prcr$gene_id))
background_genes <- unique(as.character(rownames(dds)))

go_pd <- enrichGO(
  gene = pd_genes,
  universe = background_genes,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.20,
  readable = TRUE
)

go_prcr <- enrichGO(
  gene = prcr_genes,
  universe = background_genes,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.20,
  readable = TRUE
)

go_pd_df <- as.data.frame(go_pd)
go_prcr_df <- as.data.frame(go_prcr)

write.csv(
  go_pd_df,
  "Results/GO_BP_higher_in_PD.csv",
  row.names = FALSE
)

write.csv(
  go_prcr_df,
  "Results/GO_BP_higher_in_PRCR.csv",
  row.names = FALSE
)

go_pd_plot <- dotplot(
  go_pd,
  showCategory = 10,
  title = "GO Biological Processes: Genes Higher in PD"
)

ggsave(
  "Figures/GO_BP_higher_in_PD.png",
  plot = go_pd_plot,
  width = 10,
  height = 7,
  dpi = 300
)

writeLines(
  capture.output(sessionInfo()),
  "Results/session_info.txt"
)
