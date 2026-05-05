library(DESeq2)
library(ggplot2)

counts_file   <- snakemake@input[["counts"]]
metadata_file <- snakemake@input[["metadata"]]

results_file  <- snakemake@output[["results"]]
pca_file      <- snakemake@output[["pca"]]
volcano_file  <- snakemake@output[["volcano"]]

# Read counts
counts <- read.csv(counts_file, check.names = FALSE, row.names = 1)
counts <- as.matrix(counts)
storage.mode(counts) <- "integer"

#Read metadata
meta <- read.csv(metadata_file, stringsAsFactors = FALSE)

#Keep one row per GSM
meta <- meta[!duplicated(meta$gsm), ]

# Keep only samples present in counts
meta <- meta[meta$gsm %in% colnames(counts), ]
counts <- counts[, meta$gsm, drop = FALSE]
meta <- meta[match(colnames(counts), meta$gsm), , drop = FALSE]

stopifnot(all(meta$gsm == colnames(counts)))

# Condition factor
meta$condition <- factor(meta$condition, levels = c("control", "PDAC"))
rownames(meta) <- meta$gsm

# Low ount filtering
keep <- rowSums(counts >= 10) >= 3
counts <- counts[keep, , drop = FALSE]

# DESeq2
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = meta,
  design    = ~ condition
)

dds <- estimateSizeFactors(dds, type = "poscounts")
dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "PDAC", "control"), alpha = 0.05)

# shrinkage
res_out <- as.data.frame(res)

res_shrunk <- lfcShrink(dds, coef = "condition_PDAC_vs_control", type = "apeglm")
res_out <- as.data.frame(res_shrunk)


res_out$miRNA <- rownames(res_out)
res_out <- res_out[order(res_out$padj, na.last = TRUE), ]

write.csv(res_out, results_file, row.names = FALSE)

# PCA plot
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pct_var <- round(100 * attr(pca_data, "percentVar"))

pdf(pca_file, width = 8, height = 6)
print(
  ggplot(pca_data, aes(PC1, PC2, color = condition, label = name)) +
    geom_point(size = 3, alpha = 0.9) +
    geom_text(vjust = -0.7, size = 2.5, check_overlap = TRUE) +
    labs(
      title = "PCA plot",
      x = paste0("PC1: ", pct_var[1], "% variance"),
      y = paste0("PC2: ", pct_var[2], "% variance")
    ) +
    theme_minimal(base_size = 12)
)
dev.off()

# Volcano plot
volc <- res_out
volc$neglog10padj <- -log10(volc$padj)
volc$group <- "NS"
volc$group[!is.na(volc$padj) & volc$padj < 0.05 & volc$log2FoldChange > 0] <- "Up in PDAC"
volc$group[!is.na(volc$padj) & volc$padj < 0.05 & volc$log2FoldChange < 0] <- "Down in PDAC"

pdf(volcano_file, width = 8, height = 6)
print(
  ggplot(volc, aes(log2FoldChange, neglog10padj, color = group)) +
    geom_point(alpha = 0.7, size = 1.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "solid") +
    labs(
      title = "Volcano plot",
      x = "log2 fold change (PDAC vs control)",
      y = expression(-log[10](adjusted~p))
    ) +
    theme_minimal(base_size = 12)
)
dev.off()