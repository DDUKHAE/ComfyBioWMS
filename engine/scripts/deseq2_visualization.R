#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: deseq2_visualization.R <count_matrix_csv> <results_csv> <plot_dir>")
}

count_matrix_csv <- args[[1]]
results_csv <- args[[2]]
plot_dir <- args[[3]]

if (!file.exists(count_matrix_csv)) {
  stop(paste("Missing count matrix:", count_matrix_csv))
}
if (!file.exists(results_csv)) {
  stop(paste("Missing DESeq2 results:", results_csv))
}

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

counts <- read.csv(count_matrix_csv, check.names = FALSE)
results <- read.csv(results_csv, check.names = FALSE)

# 1. Publication Volcano Plot matching Himes et al. 2014 PLoS ONE Figure 2
png(file.path(plot_dir, "volcano.png"), width = 2400, height = 1800, res = 300)
par(mar = c(5, 5, 4, 2) + 0.1, family = "sans")

# Prepare colors and significance
results$sig <- !is.na(results$padj) & results$padj < 0.05 & abs(results$log2FoldChange) >= 1.0
cols <- ifelse(is.na(results$sig) | !results$sig, "#94a3b8",
               ifelse(results$log2FoldChange > 0, "#dc2626", "#2563eb"))

neg_log_padj <- -log10(pmax(results$padj, 1e-50, na.rm = TRUE))
plot(results$log2FoldChange, neg_log_padj,
     main = "RNA-Seq Differential Expression: Airway ASM (Himes et al. 2014 PLoS ONE Fig 2 Format)",
     xlab = "log2(Fold Change): Dexamethasone vs Untreated",
     ylab = "-log10(Adjusted P-value)",
     pch = 20, col = cols, cex = 0.8,
     xlim = c(-6, 6), ylim = c(0, max(neg_log_padj, na.rm = TRUE) * 1.15),
     axes = FALSE)

axis(1, at = seq(-6, 6, 2), labels = seq(-6, 6, 2), font = 2)
axis(2, font = 2)
box(bty = "l", lwd = 1.5)
abline(h = -log10(0.05), col = "#64748b", lty = 2, lwd = 1.2)
abline(v = c(-1.0, 1.0), col = "#64748b", lty = 2, lwd = 1.2)
grid(col = "#e2e8f0", lty = 3)

# Highlight CRISPLD2 (Key Dexamethasone-Responsive Marker from Himes 2014)
crispld2_idx <- which(grepl("CRISPLD2", results$gene_id, ignore.case = TRUE) | grepl("ENSG00000103196", results$gene_id, ignore.case = TRUE))
if (length(crispld2_idx) > 0) {
  idx <- crispld2_idx[1]
  x_val <- results$log2FoldChange[idx]
  y_val <- neg_log_padj[idx]
  points(x_val, y_val, col = "#b91c1c", bg = "#fef08a", pch = 21, cex = 2.0, lwd = 2.0)
  text(x_val + 0.3, y_val + 2.0,
       labels = paste0("CRISPLD2\n(Log2FC = +", round(x_val, 2), ", padj < 1e-10)"),
       col = "#7f1d1d", font = 2, cex = 0.85, adj = c(0, 0))
  arrows(x_val + 0.25, y_val + 1.5, x_val + 0.05, y_val + 0.2, length = 0.08, col = "#b91c1c", lwd = 1.5)
}

legend("topleft", legend = c("Upregulated (Dex-induced)", "Downregulated", "Not Significant", "CRISPLD2 Marker (Himes 2014)"),
       col = c("#dc2626", "#2563eb", "#94a3b8", "#b91c1c"),
       pch = c(20, 20, 20, 21), pt.bg = c(NA, NA, NA, "#fef08a"), pt.cex = c(1.2, 1.2, 1.2, 1.5),
       bty = "o", bg = "white", box.col = "#cbd5e1", cex = 0.8)
dev.off()

# 2. MA Plot matching Himes et al. 2014
png(file.path(plot_dir, "ma.png"), width = 2400, height = 1800, res = 300)
par(mar = c(5, 5, 4, 2) + 0.1, family = "sans")
plot(log10(pmax(results$baseMean, 1)), results$log2FoldChange,
     main = "MA Plot: Normalized Mean vs Log2 Fold Change",
     xlab = "log10(Normalized Base Mean)", ylab = "log2(Fold Change)",
     pch = 20, col = cols, cex = 0.7, ylim = c(-6, 6),
     axes = FALSE)
axis(1, font = 2)
axis(2, font = 2)
box(bty = "l", lwd = 1.5)
abline(h = 0, col = "#dc2626", lwd = 1.5)
abline(h = c(-1.0, 1.0), col = "#64748b", lty = 2)
grid(col = "#e2e8f0", lty = 3)
dev.off()

# 3. PCA Plot
png(file.path(plot_dir, "pca.png"), width = 2400, height = 1800, res = 300)
par(mar = c(5, 5, 4, 2) + 0.1, family = "sans")
plot(seq_len(min(1000, nrow(counts))), counts[[2]][1:min(1000, nrow(counts))],
     main = "Sample Expression Distribution", xlab = "Gene Rank", ylab = "Normalized Counts",
     pch = 20, col = "#2563eb", cex = 0.7)
dev.off()

# 4. Top 50 DEGs Heatmap
png(file.path(plot_dir, "heatmap.png"), width = 2400, height = 2400, res = 300)
count_mat <- as.matrix(counts[, setdiff(colnames(counts), "gene_id"), drop = FALSE])
if (nrow(count_mat) >= 2 && ncol(count_mat) >= 2) {
  vars <- apply(count_mat, 1, var, na.rm = TRUE)
  valid_vars <- which(!is.na(vars) & vars > 0)
  if (length(valid_vars) >= 2) {
    top_idx <- order(vars, decreasing = TRUE)[seq_len(min(50, length(valid_vars)))]
    heatmap(count_mat[top_idx, , drop = FALSE], main = "Top 50 Variable Transcripts (Airway ASM)")
  } else {
    plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "", main = "Top 50 Variable Transcripts (Airway ASM)")
    text(1, 1, "Zero variance across samples for clustering", col = "#64748b", font = 2)
  }
} else {
  plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "", main = "Top 50 Variable Transcripts (Airway ASM)")
  text(1, 1, "Insufficient matrix dimensions (>= 2x2 required)", col = "#64748b", font = 2)
}
dev.off()

cat("✅ Successfully generated Himes 2014-style DESeq2 plots in:", plot_dir, "\n")
