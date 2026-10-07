#################################################################################
# Identifying Candidate Biomarkers of Psoriasis
# Convergence across differential expression, cross-validation selection
# stability, and machine-learning (SHAP) feature importance.
#
# Author : Cailey Coghlan  (MSc Data Science, University of Sheffield, 2026)
# Data   : GEO GSE54456 (Li et al., 2014): 92 psoriatic + 82 normal skin biopsies
#
# HOW TO RUN
#   1. Open psoriasis-biomarker-discovery.Rproj in RStudio (this sets the
#      working directory to the repository root), or setwd() to the repo root.
#   2. Run R/00_install_packages.R once.
#   3. Run this script top to bottom: source("R/01_analysis.R")
#
# Every figure is written to figures/ and every table to tables/.
# Runtime: about 30-60 minutes on a laptop. The SVM Kernel SHAP step (Step 28)
# is the slowest.
#
# STEP 1 (package installation) lives in R/00_install_packages.R.
#################################################################################

################################################
# STEP 2: LOAD LIBRARIES
################################################

library(nestedcv)
library(limma)
library(glmnet)
library(ranger)
library(xgboost)
library(caret)
library(e1071)
library(shapviz)
library(pROC)
library(ggplot2)
library(ggrepel)
library(matrixStats)
library(mltools)
library(kernlab)
library(treeshap)
library(kernelshap)

################################################
# OUTPUT HELPER
# Routes every saved file into the right folder:
#   .png -> figures/     .csv -> tables/
################################################

dir.create("figures", showWarnings = FALSE)
dir.create("tables",  showWarnings = FALSE)

out <- function(filename) {
  folder <- if (grepl("\\.png$", filename)) "figures" else "tables"
  file.path(folder, filename)
}

################################################
# STEP 3: SET RANDOM SEED
################################################

set.seed(42)

################################################
# STEP 4: LOAD DATA
################################################

# Load RPKM expression matrix
# Genes are rows, samples are columns
rpkm <- read.delim(
  file.path("data", "raw", "RPKM_Values.txt"),
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)

# Load sample ID and tissue type CSV table
tissue_type_and_sample_ID <- read.csv(
  file.path("data", "raw", "sample_metadata.csv"),
  header = TRUE,
  stringsAsFactors = FALSE
)

# Confirm dimensions
cat("RPKM matrix dimensions (genes x samples):", dim(rpkm), "\n")
cat("sample ID and tissue type CSV table dimensions:", dim(tissue_type_and_sample_ID), "\n")

################################################
# STEP 5: VERIFY SAMPLE NAMES MATCH
################################################

matched <- sum(colnames(rpkm) %in% tissue_type_and_sample_ID$Sample.ID)
cat("Samples in RPKM matrix:", ncol(rpkm), "\n")
cat("Samples in sample ID and tissue type table:", nrow(tissue_type_and_sample_ID), "\n")
cat("Samples that match:", matched, "\n")

################################################
# STEP 6: TRANSPOSE AND MERGE
################################################

# Transpose so rows are samples and columns are genes
rpkm_transpose <- as.data.frame(t(rpkm))

# Move row names into a column for merging
rpkm_transpose$Sample.ID <- rownames(rpkm_transpose)

# Merge with tissue_type_and_sample_ID table by sample ID
genes_table <- merge(
  rpkm_transpose,
  tissue_type_and_sample_ID[, c("Sample.ID", "Tissue.Type")],
  by = "Sample.ID",
  all = FALSE
)

# Convert Tissue Type to factor with Normal as reference level
genes_table$Tissue.Type <- factor(genes_table$Tissue.Type,
                            levels = c("Normal", "Psoriasis"))

# Confirm final dimensions and Tissue Type distribution
cat("\nMerged dataset dimensions:", nrow(genes_table), "x", ncol(genes_table), "\n")
head(genes_table[, c(1:6, ncol(genes_table))])
cat("\nLabel distribution:\n")
print(table(genes_table$Tissue.Type))

################################################
# STEP 7: CHECK FOR NAs AND ZEROS
################################################

# Check for NAs in Sample.ID and Tissue.Type first
na_sample_id <- sum(is.na(genes_table$Sample.ID))
na_tissue_type <- sum(is.na(genes_table$Tissue.Type))
cat("NAs in Sample.ID:", na_sample_id, "\n")
cat("NAs in Tissue.Type:", na_tissue_type, "\n")

# --- Check for duplicate sample IDs ---
# Duplicates would inflate the sample count and could place the same
# sample in both training and test sets. Confirm each ID is unique.
n_dup <- sum(duplicated(genes_table$Sample.ID))
cat("Duplicate Sample.IDs:", n_dup, "\n")
if (n_dup > 0) {
  cat("WARNING - the following Sample.IDs are duplicated:\n")
  print(genes_table$Sample.ID[duplicated(genes_table$Sample.ID)])
  stop("Resolve duplicate sample IDs before proceeding\n")
} else {
  cat("All Sample.IDs are unique — no duplicates\n")
}

# Extract just the gene expression columns
# Exclude Sample.ID and Tissue.Type
gene_cols <- genes_table[, !colnames(genes_table) %in% c("Sample.ID", "Tissue.Type")]

# Count NAs in gene expression matrix
total_nas <- sum(is.na(gene_cols))
cat("Total NAs in gene expression matrix:", total_nas, "\n")

# Count zeros in gene expression matrix
total_zeros <- sum(gene_cols == 0)
total_values <- nrow(gene_cols) * ncol(gene_cols)
prop_zeros <- round((total_zeros / total_values) * 100, 2)
cat("Total zeros:", total_zeros, "\n")
cat("Proportion of zeros:", prop_zeros, "%\n")

# Remove genes with zero expression across all samples
# These contain no information and their removal is
# independent of sample partitioning
all_zero_mask <- colSums(gene_cols) == 0
cat("\nGenes with zero expression across all samples:", sum(all_zero_mask), "\n")
gene_cols <- gene_cols[, !all_zero_mask]
cat("Genes remaining after removing all-zero genes:", ncol(gene_cols), "\n")

# Rebuild genes_table with zero genes removed
genes_table <- cbind(
  Sample.ID   = genes_table$Sample.ID,
  gene_cols,
  Tissue.Type = genes_table$Tissue.Type
)
genes_table$Tissue.Type <- factor(genes_table$Tissue.Type,
                                  levels = c("Normal", "Psoriasis"))

cat("genes_table dimensions after cleaning:",
    nrow(genes_table), "x", ncol(genes_table), "\n")

# Check for samples with zero expression across all genes
all_zero_samples <- rowSums(gene_cols) == 0
cat("Samples with zero expression across all genes:", sum(all_zero_samples), "\n")

if (sum(all_zero_samples) > 0) {
  cat("WARNING - the following samples have zero expression across all genes:\n")
  print(genes_table$Sample.ID[all_zero_samples])
  stop("Remove or investigate these samples before proceeding\n")
} else {
  cat("No samples with zero expression across all genes — data looks clean\n")
}

################################################
# STEP 8: DESCRIPTIVE STATISTICS ON RAW RPKM
################################################

# Overall statistics across entire expression matrix
all_values <- unlist(gene_cols)

cat("\n--- Overall Expression Statistics (RPKM) ---\n")
cat("Mean:  ", round(mean(all_values), 4), "\n")
cat("Median:", round(median(all_values), 4), "\n")
cat("SD:    ", round(sd(all_values), 4), "\n")
cat("Min:   ", round(min(all_values), 4), "\n")
cat("Max:   ", round(max(all_values), 4), "\n")

# Class balance
cat("\n--- Class Balance ---\n")
print(table(genes_table$Tissue.Type))
cat("Proportion Psoriasis:",
    round(sum(genes_table$Tissue.Type == "Psoriasis") /
            nrow(genes_table) * 100, 1), "%\n")
cat("Proportion Normal:",
    round(sum(genes_table$Tissue.Type == "Normal") /
            nrow(genes_table) * 100, 1), "%\n")

################################################
# STEP 9: EDA — UNTRANSFORMED RPKM GENE EXPRESSION DATA
################################################

# Reshape gene_cols to long format for ggplot
# This creates one row per gene per sample
rpkm_long <- data.frame(
  Sample.ID   = rep(genes_table$Sample.ID, times = ncol(gene_cols)),
  Tissue.Type = rep(genes_table$Tissue.Type, times = ncol(gene_cols)),
  Expression  = unlist(gene_cols)
)

# Order samples by tissue type then by median expression
# This groups Normal and Psoriasis samples together on the plot
sample_order <- genes_table$Sample.ID[order(genes_table$Tissue.Type)]
rpkm_long$Sample.ID <- factor(rpkm_long$Sample.ID, levels = sample_order)

rpkm_long$Tissue.Type <- factor(rpkm_long$Tissue.Type,
                                levels = c("Normal", "Psoriasis"))

# --- Density plot of untransformed RPKM distribution ---
rpkm_long$Tissue.Type <- factor(rpkm_long$Tissue.Type,
                                levels = c("Psoriasis", "Normal"))
density_raw <- ggplot(rpkm_long, aes(x = Expression,
                                     fill = Tissue.Type,
                                     colour = Tissue.Type)) +
  geom_density(alpha     = 0.4,
               linewidth = 0.8) +
  scale_fill_manual(
    values = c("Normal"    = "#0072B2",
               "Psoriasis" = "#E69F00"),
    breaks = c("Normal", "Psoriasis"),
    name   = "Tissue Type"
  ) +
  scale_colour_manual(
    values = c("Normal"    = "#0072B2",
               "Psoriasis" = "#E69F00"),
    breaks = c("Normal", "Psoriasis"),
    name   = "Tissue Type"
  ) +
  coord_cartesian(xlim = c(0, 490)) +
  scale_x_continuous(labels = scales::comma) +
  labs(
    title    = "Untransformed RPKM Gene Expression Density Distribution",
    subtitle = "Extreme right skew confirms the need for log2(x+1) transformation.\nX-axis limited to 0-500 for readability — values extend to 35,651.",
    x        = "Untransformed RPKM Gene Expression Value",
    y        = "Density",
    caption  = "All 21,099 genes across 174 samples"
  ) +
  theme_classic(base_size   = 12,
                base_family = "Arial") +
  theme(
    plot.title         = element_text(face = "bold", size = 14),
    plot.subtitle      = element_text(size = 12, colour = "grey40"),
    legend.position    = "top",
    legend.title       = element_text(face = "bold"),
    panel.grid.major.y = element_line(colour = "grey90"),
    plot.caption       = element_text(size = 12, colour = "grey50")
  )

# Save density plot
ggsave(out("density_rpkm_untransformed.png"),
       plot   = density_raw,
       width  = 10,
       height = 6,
       dpi    = 300,
       bg     = "white")
cat("Untransformed RPKM density plot saved\n")

################################################
# STEP 10: LOG2(x+1) TRANSFORMATION
################################################

# Apply log2(x+1) to all gene expression columns only
gene_cols_log2 <- log2(gene_cols + 1)

# Rebuild genes_table with log2 transformed values
# Preserve Sample.ID and Tissue.Type columns
genes_table_log2 <- cbind(
  Sample.ID   = genes_table$Sample.ID,
  gene_cols_log2,
  Tissue.Type = genes_table$Tissue.Type
)
genes_table_log2$Tissue.Type <- factor(genes_table_log2$Tissue.Type,
                                       levels = c("Normal", "Psoriasis"))

cat("Log2(x+1) transformation applied\n")
cat("Transformed dataset dimensions:",
    nrow(genes_table_log2), "x", ncol(genes_table_log2), "\n")

###############################################################
# STEP 11: DESCRIPTIVE STATISTICS ON TRANSFORMED DATA
###############################################################

# Overall statistics on transformed data
all_values_log2 <- unlist(gene_cols_log2)

cat("\n--- Overall Expression Statistics (Log2 Transformed) ---\n")
cat("Mean:  ", round(mean(all_values_log2), 4), "\n")
cat("Median:", round(median(all_values_log2), 4), "\n")
cat("SD:    ", round(sd(all_values_log2), 4), "\n")
cat("Min:   ", round(min(all_values_log2), 4), "\n")
cat("Max:   ", round(max(all_values_log2), 4), "\n")

# Per-sample statistics on transformed data
per_sample_stats_log2 <- data.frame(
  Sample.ID   = genes_table_log2$Sample.ID,
  Tissue.Type = genes_table_log2$Tissue.Type,
  Mean        = round(rowMeans(gene_cols_log2), 4),
  Median      = round(apply(gene_cols_log2, 1, median), 4),
  SD          = round(apply(gene_cols_log2, 1, sd), 4),
  Min         = round(apply(gene_cols_log2, 1, min), 4),
  Max         = round(apply(gene_cols_log2, 1, max), 4)
)

# Compare before and after transformation
cat("\n--- Skew Comparison: Before vs After Transformation ---\n")
cat("Untransformed — Mean:", round(mean(unlist(gene_cols)), 4),
    "| Median:", round(median(unlist(gene_cols)), 4), "\n")
cat("Log2 transformed — Mean:", round(mean(all_values_log2), 4),
    "| Median:", round(median(all_values_log2), 4), "\n")
cat("Mean-median gap before:", 
    round(mean(unlist(gene_cols)) - median(unlist(gene_cols)), 4), "\n")
cat("Mean-median gap after:", 
    round(mean(all_values_log2) - median(all_values_log2), 4), "\n")

################################################
# STEP 12A: EDA — LOG2 TRANSFORMED DATA
################################################

# Reshape transformed data to long format
log2_long <- data.frame(
  Sample.ID   = rep(genes_table_log2$Sample.ID,
                    times = ncol(gene_cols_log2)),
  Tissue.Type = rep(genes_table_log2$Tissue.Type,
                    times = ncol(gene_cols_log2)),
  Expression  = unlist(gene_cols_log2)
)

# Use same sample order as untransformed plots for direct comparison
log2_long$Sample.ID <- factor(log2_long$Sample.ID,
                              levels = sample_order)

log2_long$Tissue.Type <- factor(log2_long$Tissue.Type,
                                levels = c("Normal", "Psoriasis"))

# --- Boxplot of log2 transformed values ---
boxplot_log2 <- ggplot(log2_long, aes(x = Sample.ID,
                                      y = Expression,
                                      fill = Tissue.Type)) +
  geom_boxplot(outlier.size  = 0.2,
               outlier.alpha = 0.2,
               linewidth     = 0.3) +
  scale_fill_manual(
    values = c("Normal"    = "#0072B2",
               "Psoriasis" = "#E69F00"),
    name = "Tissue Type"
  ) +
  facet_wrap(~ Tissue.Type,
             scales = "free_x",
             ncol   = 2) +
  labs(
    title    = "Log2(x+1) Transformed RPKM Gene Expression Distribution Per Sample",
    subtitle = "Consistent medians confirm transformation success.\nReduced spread compared to untransformed data.",
    x        = "Sample",
    y        = "Log2(x+1) RPKM Gene Expression Value",
    caption  = "Each box represents one sample. Normal n = 82, Psoriasis n = 92."
  ) +
  theme_classic(base_size   = 12,
                base_family = "Arial") +
  theme(
    axis.text.x        = element_blank(),
    axis.ticks.x       = element_blank(),
    plot.title         = element_text(face = "bold", size = 14),
    plot.subtitle      = element_text(size = 12, colour = "grey40"),
    legend.position    = "none",
    strip.text         = element_text(face = "bold", size = 12),
    strip.background   = element_rect(fill = "grey95", colour = "grey70"),
    panel.grid.major.y = element_line(colour = "grey90"),
    plot.caption       = element_text(size = 12, colour = "grey50")
  )

# Save log2 boxplot
ggsave(out("boxplot_log2_transformed.png"),
       plot   = boxplot_log2,
       width  = 16,
       height = 7,
       dpi    = 300,
       bg     = "white")
cat("Log2 transformed boxplot saved\n")

# --- Density plot of log2 transformed values ---

# Draw Psoriasis first so Normal appears on top
log2_long$Tissue.Type <- factor(log2_long$Tissue.Type,
                                levels = c("Psoriasis", "Normal"))

density_log2 <- ggplot(log2_long, aes(x = Expression,
                                      fill = Tissue.Type,
                                      colour = Tissue.Type)) +
  geom_density(alpha     = 0.4,
               linewidth = 0.8) +
  scale_fill_manual(
    values = c("Normal"    = "#0072B2",
               "Psoriasis" = "#E69F00"),
    breaks = c("Normal", "Psoriasis"),
    name   = "Tissue Type"
  ) +
  scale_colour_manual(
    values = c("Normal"    = "#0072B2",
               "Psoriasis" = "#E69F00"),
    breaks = c("Normal", "Psoriasis"),
    name   = "Tissue Type"
  ) +
  labs(
    title    = "Log2(x+1) Transformed RPKM Gene Expression Density Distribution",
    subtitle = "Extreme right skew of untransformed RPKM data successfully reduced by log2(x+1) transformation",
    x        = "Log2(x+1) RPKM Gene Expression Value",
    y        = "Density",
    caption  = "All 21,099 genes across 174 samples"
  ) +
  theme_classic(base_size   = 12,
                base_family = "Arial") +
  theme(
    plot.title         = element_text(face = "bold", size = 14),
    plot.subtitle      = element_text(size = 12, colour = "grey40"),
    legend.position    = "top",
    legend.title       = element_text(face = "bold"),
    panel.grid.major.y = element_line(colour = "grey90"),
    plot.caption       = element_text(size = 12, colour = "grey50")
  )

# Save log2 density plot
ggsave(out("density_log2_transformed.png"),
       plot   = density_log2,
       width  = 10,
       height = 6,
       dpi    = 300,
       bg     = "white")
cat("Log2 transformed density plot saved\n")

################################################
# STEP 12B: SAMPLE CORRELATION OUTLIER CHECK
# Relational outlier check: does each sample's overall
# expression pattern resemble the other 173? Complements
# the boxplot, which only judges each sample's distribution
# in isolation.
# Log2 data used — Pearson correlation is distorted by the
# raw RPKM skew, so log2 gives a truer similarity measure.
################################################

# --- Compute pairwise sample correlations ---
sample_cor <- cor(t(as.matrix(gene_cols_log2)))
rownames(sample_cor) <- colnames(sample_cor) <- genes_table_log2$Sample.ID

cat("\nSample correlation matrix dimensions:", dim(sample_cor), "\n")
cat("Correlation range:", round(min(sample_cor), 3), "to",
    round(max(sample_cor), 3), "\n")

# --- Score each sample by how typical it is ---
# Self-correlation is always 1 and meaningless for outlier
# detection, so blank the diagonal before averaging.
diag_backup <- sample_cor
diag(diag_backup) <- NA

# Mean correlation with all other samples. Low = resembles
# nothing else here, a candidate outlier.
mean_cor_per_sample <- rowMeans(diag_backup, na.rm = TRUE)

outlier_summary <- data.frame(
  Sample.ID   = genes_table_log2$Sample.ID,
  Tissue.Type = genes_table_log2$Tissue.Type,
  Mean.Correlation = round(mean_cor_per_sample, 4)
)
outlier_summary <- outlier_summary[order(outlier_summary$Mean.Correlation), ]

cat("\n--- 10 samples with lowest average correlation to all others ---\n")
print(head(outlier_summary, 10))

cat("\n--- Overall range of average correlation scores ---\n")
cat("Minimum:", round(min(outlier_summary$Mean.Correlation), 4), "\n")
cat("Maximum:", round(max(outlier_summary$Mean.Correlation), 4), "\n")
cat("Mean:   ", round(mean(outlier_summary$Mean.Correlation), 4), "\n")
cat("\nInspect the table above for any sample showing a marked gap\n")
cat("from the rest of the distribution, rather than a gradual decline\n")

# Save the full table as supporting evidence for the
# dissertation
write.csv(outlier_summary, out("sample_correlation_summary.csv"), row.names = FALSE)
cat("Sample correlation summary saved to sample_correlation_summary.csv\n")

################################################
# STEP 12C: PCA — MULTIVARIATE STRUCTURE CHECK
# Purpose: a third view that the
# per-sample boxplot and pairwise correlation check cannot
# provide. Projects all 174 samples onto their two largest
# axes of variation and colours by tissue type, to confirm
# (a) that psoriatic and normal samples separate, and
# (b) that no large unexplained structure dominates the data.
# Log2 data used.
################################################

pca <- prcomp(gene_cols_log2, center = TRUE, scale. = FALSE)

# Percentage of total variance captured by each component
var_explained <- pca$sdev^2 / sum(pca$sdev^2) * 100
pc1_var <- round(var_explained[1], 1)
pc2_var <- round(var_explained[2], 1)

cat("\n--- PCA Variance Explained ---\n")
cat("PC1:", pc1_var, "%\n")
cat("PC2:", pc2_var, "%\n")
cat("PC1 + PC2:", round(pc1_var + pc2_var, 1), "%\n")

# Build a scores table: one row per sample, its tissue type,
# and its coordinates on the first two components
pca_scores <- data.frame(
  Sample.ID   = genes_table_log2$Sample.ID,
  Tissue.Type = factor(genes_table_log2$Tissue.Type,
                       levels = c("Normal", "Psoriasis")),
  PC1         = pca$x[, 1],
  PC2         = pca$x[, 2]
)

# --- PCA scores plot ---
pca_plot <- ggplot(pca_scores, aes(x = PC1,
                                   y = PC2,
                                   colour = Tissue.Type)) +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "grey80", linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey80", linewidth = 0.4) +
  stat_ellipse(linewidth = 0.7, linetype = "dashed") +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_text_repel(
    data = subset(pca_scores, Sample.ID == "M3171"),
    aes(label = Sample.ID),
    colour             = "black",
    size               = 3.2,
    fontface           = "bold",
    box.padding        = 0.6,
    point.padding      = 0.4,
    segment.colour     = "grey50",
    segment.size       = 0.3,
    min.segment.length = 0,
    show.legend        = FALSE
  ) +
  scale_colour_manual(
    values = c("Normal"    = "#0072B2",
               "Psoriasis" = "#E69F00"),
    breaks = c("Normal", "Psoriasis"),
    name   = "Tissue Type"
  ) +
  labs(
    title    = "PCA Graph of Log2(x+1) RPKM Gene Expression Across All Samples",
    subtitle = "Separation along PC1 indicates the disease signal dominates the variation.",
    x        = paste0("PC1 (", pc1_var, "% of variance)"),
    y        = paste0("PC2 (", pc2_var, "% of variance)"),
    caption  = "Each point is one sample (n = 174; 82 Normal, 92 Psoriasis) across 21,099 genes.\nDashed ellipses: 95% confidence region per group. M3171 is an atypical but retained psoriasis sample."
  ) +
  theme_classic(base_size   = 12,
                base_family = "Arial") +
  theme(
    plot.title         = element_text(face = "bold", size = 14),
    plot.subtitle      = element_text(size = 12, colour = "grey40"),
    legend.position    = "top",
    legend.title       = element_text(face = "bold"),
    panel.grid.major.y = element_line(colour = "grey90"),
    plot.caption       = element_text(size = 12, colour = "grey50")
  )

# Save PCA plot
ggsave(out("pca_plot.png"),
       plot   = pca_plot,
       width  = 9,
       height = 7,
       dpi    = 300,
       bg     = "white")
cat("\nPCA plot saved to pca_plot.png\n")

# --- Identify the most extreme sample on each PC for cross-checking ---
# Names the visible outliers on the PCA plot so they can be looked up
# against the Step 12B correlation table.

cat("\n--- Most extreme samples on PC2 (the visible outlier axis) ---\n")
print(head(pca_scores[order(pca_scores$PC2), c("Sample.ID", "Tissue.Type", "PC1", "PC2")], 3))

cat("\n--- Most extreme samples on PC1 ---\n")
print(head(pca_scores[order(pca_scores$PC1), c("Sample.ID", "Tissue.Type", "PC1", "PC2")], 3))

################################################
# STEP 13: LIMMA DIFFERENTIAL EXPRESSION ANALYSIS
# Purpose: Answer RQ1
################################################

# Prepare expression matrix for limma
# Limma requires genes as rows and samples as columns
expr_matrix <- t(as.matrix(gene_cols_log2))

# Confirm dimensions. Should be 21099 genes x 174 samples
cat("Expression matrix dimensions (genes x samples):",
    dim(expr_matrix), "\n")

# Create design matrix
# Normal = 0 (reference level), Psoriasis = 1
group <- genes_table_log2$Tissue.Type
design <- model.matrix(~ group)

cat("\nDesign matrix dimensions:", dim(design), "\n")
cat("First few rows of design matrix:\n")
print(head(design))

# Fit linear model to each gene
fit <- lmFit(expr_matrix, design)
cat("\nLinear model fitted to", nrow(fit$coefficients), "genes\n")

# Apply empirical Bayes shrinkage
# Stabilises variance estimates across genes
# Essential for small sample sizes relative to gene number
fit <- eBayes(fit)
cat("Empirical Bayes shrinkage applied\n")

# Extract results for all genes using topTable
# coef = 2 extracts the Psoriasis vs Normal comparison
# (second coefficient in the design matrix)
# number = Inf returns all 21,099 genes not just the top ones
# adjust.method = "BH" applies Benjamini-Hochberg FDR correction
# to control false discovery rate across all genes
# sort.by = "P" sorts all genes by FDR-adjusted p-value
# most statistically significant genes appear first
# Significant genes will be re-sorted by log2FC afterwards
limma_results <- topTable(
  fit,
  coef          = 2,
  number        = Inf,
  adjust.method = "BH",
  sort.by       = "P"
)

# Add gene names as a column
limma_results$gene <- rownames(limma_results)

# Confirm output
cat("\nLimma results dimensions:", dim(limma_results), "\n")

# Summary of results
cat("\n--- Limma Results Summary ---\n")
cat("Total genes tested:", nrow(limma_results), "\n")

# Apply significance thresholds
# Benjamini-Hochberg FDR-adjusted p-value < 0.05
# And absolute log2FC > 1 (at least twofold change)
sig_genes <- limma_results[
  limma_results$adj.P.Val < 0.05 &
    abs(limma_results$logFC) > 1, ]

# Split into upregulated and downregulated
# Order each by log2FC — most dramatic changes first
sig_genes_up <- sig_genes[sig_genes$logFC > 1, ]
sig_genes_up <- sig_genes_up[order(sig_genes_up$logFC,
                                   decreasing = TRUE), ]

sig_genes_down <- sig_genes[sig_genes$logFC < -1, ]
sig_genes_down <- sig_genes_down[order(sig_genes_down$logFC,
                                       decreasing = FALSE), ]

# RQ1 summary statement
cat("\n--- RQ1 Summary Statement ---\n")
# Check how many genes have adjusted p-value below 0.05
# regardless of log2FC threshold
sum(limma_results$adj.P.Val < 0.05)
cat("Genes with adjusted p-value < 0.05:",
    sum(limma_results$adj.P.Val < 0.05), "\n")
cat("Out of total genes tested:", nrow(limma_results), "\n")
cat("Percentage:",
    round(sum(limma_results$adj.P.Val < 0.05) /
            nrow(limma_results) * 100, 1), "%\n")
cat("Differential expression analysis using limma identified",
    nrow(sig_genes),
    "significantly differentially expressed genes\n")
cat("meeting the threshold of FDR-adjusted p-value < 0.05",
    "and |log2FC| > 1.\n")
cat(sum(sig_genes$logFC > 1),
    "genes were upregulated in psoriatic skin relative to normal skin.\n")
cat(sum(sig_genes$logFC < -1),
    "genes were downregulated in psoriatic skin relative to normal skin.\n")
cat("\nMost dramatically upregulated gene:",
    sig_genes_up$gene[1],
    "| log2FC =", round(sig_genes_up$logFC[1], 2), "\n")
cat("Most dramatically downregulated gene:",
    sig_genes_down$gene[1],
    "| log2FC =", round(sig_genes_down$logFC[1], 2), "\n")

# Print top 10 most upregulated
cat("\nTop 10 most upregulated genes in psoriasis (by log2FC):\n")
temp_up <- head(sig_genes_up[, c("gene", "logFC", "adj.P.Val")], 10)
colnames(temp_up)[colnames(temp_up) == "logFC"] <- "log2FC"
print(temp_up)
rm(temp_up)

# Print top 10 most downregulated
cat("\nTop 10 most downregulated genes in psoriasis (by log2FC):\n")
temp_down <- head(sig_genes_down[, c("gene", "logFC", "adj.P.Val")], 10)
colnames(temp_down)[colnames(temp_down) == "logFC"] <- "log2FC"
print(temp_down)
rm(temp_down)

# Save significant gene names for later SHAP cross-reference
# Contains all 1131 significant genes — not just top 10
sig_gene_names <- sig_genes$gene
cat("\nSignificant gene list saved for SHAP cross-reference\n")
cat("Total significant genes:", length(sig_gene_names), "\n")


# SAVE CSV FILES

# Full results table — ALL 21,099 genes regardless of significance
# Includes genes with any log2FC and any adjusted p-value
# Use this file to look up any specific gene's log2FC and adjusted p-value
limma_results_export <- limma_results
colnames(limma_results_export)[
  colnames(limma_results_export) == "logFC"] <- "log2FC"
write.csv(limma_results_export,
          out("limma_results_all_genes.csv"),
          row.names = FALSE)
cat("\nFull limma results saved to limma_results_all_genes.csv\n")
cat("Contains all", nrow(limma_results_export),
    "genes regardless of statistical significance\n")

# Significantly upregulated genes only
# Filtered to: FDR-adjusted p-value < 0.05 AND log2FC > 1
# Sorted by log2FC — most dramatically upregulated genes first
# These genes show at least twofold higher expression in psoriasis
sig_genes_up_export <- sig_genes_up
colnames(sig_genes_up_export)[
  colnames(sig_genes_up_export) == "logFC"] <- "log2FC"
write.csv(sig_genes_up_export,
          out("limma_sig_genes_upregulated.csv"),
          row.names = FALSE)
cat("\nUpregulated genes saved to limma_sig_genes_upregulated.csv\n")
cat("Contains", nrow(sig_genes_up_export),
    "genes upregulated in psoriasis",
    "(FDR-adjusted p-value < 0.05, log2FC > 1)\n")

# Significantly downregulated genes only
# Filtered to: FDR-adjusted p-value < 0.05 AND log2FC < -1
# Sorted by log2FC — most dramatically downregulated genes first
# These genes show at least twofold lower expression in psoriasis
sig_genes_down_export <- sig_genes_down
colnames(sig_genes_down_export)[
  colnames(sig_genes_down_export) == "logFC"] <- "log2FC"
write.csv(sig_genes_down_export,
          out("limma_sig_genes_downregulated.csv"),
          row.names = FALSE)
cat("\nDownregulated genes saved to limma_sig_genes_downregulated.csv\n")
cat("Contains", nrow(sig_genes_down_export),
    "genes downregulated in psoriasis",
    "(FDR-adjusted p-value < 0.05, log2FC < -1)\n")

################################################
# STEP 14: VOLCANO PLOT
# VISUALISATION FOR ANSWERING RQ1
################################################

# Create a column marking significance status for colouring
limma_results$significance <- "Not Significant"
limma_results$significance[
  limma_results$adj.P.Val < 0.05 &
    limma_results$logFC > 1] <- "Upregulated"
limma_results$significance[
  limma_results$adj.P.Val < 0.05 &
    limma_results$logFC < -1] <- "Downregulated"

cat("Volcano plot significance breakdown:\n")
print(table(limma_results$significance))

# Confirm how many genes pass p-value threshold alone
# vs both thresholds together — important methodological point
cat("\nGenes with adjusted p-value < 0.05 (any log2FC):",
    sum(limma_results$adj.P.Val < 0.05), "\n")
cat("Genes meeting BOTH thresholds (adj p < 0.05 AND |log2FC| > 1):",
    nrow(sig_genes), "\n")

# Genes to label: literature review genes plus top 10 up/down by log2FC
literature_genes <- c("S100A8", "S100A9", "DEFB4A", "CXCL8",
                      "SERPINB4", "WIF1", "BTC")

top10_up_genes <- head(sig_genes_up$gene, 10)
top10_down_genes <- head(sig_genes_down$gene, 10)

known_genes <- unique(c(literature_genes,
                        top10_up_genes,
                        top10_down_genes))

known_genes_data <- limma_results[
  limma_results$gene %in% known_genes, ]

# Build the volcano plot
volcano_plot <- ggplot(limma_results,
                       aes(x = logFC,
                           y = -log10(adj.P.Val),
                           colour = significance)) +
  geom_point(alpha = 0.5,
             size  = 0.8) +
  scale_colour_manual(
    values = c("Upregulated"     = "#E69F00",
               "Downregulated"   = "#0072B2",
               "Not Significant" = "grey70"),
    breaks = c("Upregulated", "Downregulated", "Not Significant"),
    name   = "Differential Expression"
  ) +
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4, 6, 8, 10),
                     limits = c(NA, 12)
  ) +
  geom_vline(xintercept = c(-1, 1),
             linetype   = "dashed",
             colour     = "grey40",
             linewidth  = 0.4) +
  geom_hline(yintercept = -log10(0.05),
             linetype   = "dashed",
             colour     = "grey40",
             linewidth  = 0.4) +
  annotate("text",
           x      = -1,
           y      = max(-log10(limma_results$adj.P.Val)) * 0.97,
           label  = "Log2FC = -1",
           angle  = 90,
           vjust  = -0.5,
           size   = 4.5,
           colour = "grey30") +
  annotate("text",
           x      = 1,
           y      = max(-log10(limma_results$adj.P.Val)) * 0.97,
           label  = "Log2FC = 1",
           angle  = 90,
           vjust  = -0.5,
           size   = 4.5,
           colour = "grey30") +
  annotate("text",
           x      = max(limma_results$logFC) * 0.95,
           y      = -log10(0.05),
           label  = "-Log10(FDR-Adjusted P-Value = 0.05)",
           vjust  = -0.5,
           hjust  = 1,
           size   = 4.5,
           colour = "grey30") +
  geom_point(data   = known_genes_data,
             aes(x  = logFC, y = -log10(adj.P.Val)),
             colour = "black",
             size   = 1.8,
             shape  = 21,
             fill   = "white",
             stroke = 0.8) +
  geom_text_repel(data  = known_genes_data,
                  aes(x = logFC, y = -log10(adj.P.Val), label = gene),
                  colour             = "black",
                  fontface           = "bold",
                  size               = 3.2,
                  box.padding        = 0.5,
                  point.padding      = 0.4,
                  segment.colour     = "grey50",
                  segment.size       = 0.3,
                  segment.alpha      = 0.5,
                  max.overlaps       = Inf,
                  force              = 1,
                  force_pull         = 0.5,
                  max.iter           = 20000,
                  max.time           = 2,
                  min.segment.length = 0) +
  labs(
    title    = "Volcano Plot: Differentially Expressed Genes in Psoriatic vs Normal Skin",
    subtitle = paste0("1,131 of 21,099 genes met both significance criteria ",
                      "(FDR-adjusted p < 0.05, |log2FC| > 1).\n",
                      "Labelled genes: literature review markers plus top 10 ",
                      "most upregulated and downregulated by log2FC."),
    x        = "Log2 Fold Change (Psoriasis vs Normal)",
    y        = "-Log10(FDR-Adjusted P-Value)",
    caption  = paste0("All 21,099 genes tested. 15,466 significant by p-value alone. ",
                      nrow(sig_genes), " significant by both thresholds (",
                      sum(sig_genes$logFC > 1), " up, ",
                      sum(sig_genes$logFC < -1), " down)")
  ) +
  theme_classic(base_size   = 12,
                base_family = "Arial") +
  theme(
    plot.title         = element_text(face = "bold", size = 14),
    plot.subtitle      = element_text(size = 12, colour = "grey40"),
    legend.position    = "top",
    legend.title       = element_text(face = "bold"),
    panel.grid.major.y = element_line(colour = "grey90"),
    plot.caption       = element_text(size = 12, colour = "grey50")
  )

# Save volcano plot
ggsave(out("volcano_plot.png"),
       plot   = volcano_plot,
       width  = 12,
       height = 9,
       dpi    = 300,
       bg     = "white")
cat("\nVolcano plot saved to volcano_plot.png\n")

##############################################################
# STEP 15: STRATIFIED TRAIN/TEST SPLIT
# 80/20 split stratified by Tissue.Type.
# train_index is reused in Step 17 to build a raw
# (untransformed) version for Random Forest and XGBoost. Row
# numbers stay valid because row order is unchanged.
###############################################################

set.seed(42)

# Create stratified 80/20 split 
# Stratification preserves 53:47 class ratio in both sets
train_index <- createDataPartition(
  genes_table_log2$Tissue.Type,
  p    = 0.80,
  list = FALSE
)

# Use these row numbers to split the LOG2 dataset into
# training and test sets
train_data_log2 <- genes_table_log2[train_index, ]
test_data_log2  <- genes_table_log2[-train_index, ]

# Confirm dimensions
cat("Training set dimensions:", nrow(train_data_log2), "x", ncol(train_data_log2), "\n")
cat("Test set dimensions:", nrow(test_data_log2), "x", ncol(test_data_log2), "\n")

# Confirm class ratio preserved in both sets
cat("\nTraining set class distribution:\n")
print(table(train_data_log2$Tissue.Type))
cat("Training proportion Psoriasis:",
    round(sum(train_data_log2$Tissue.Type == "Psoriasis") /
            nrow(train_data_log2) * 100, 1), "%\n")
cat("Training proportion Normal:",
    round(sum(train_data_log2$Tissue.Type == "Normal") /
            nrow(train_data_log2) * 100, 1), "%\n")

cat("\nTest set class distribution:\n")
print(table(test_data_log2$Tissue.Type))
cat("Test proportion Psoriasis:",
    round(sum(test_data_log2$Tissue.Type == "Psoriasis") /
            nrow(test_data_log2) * 100, 1), "%\n")
cat("Test proportion Normal:",
    round(sum(test_data_log2$Tissue.Type == "Normal") /
            nrow(test_data_log2) * 100, 1), "%\n")

################################################################
# STEP 16: REMOVE ZERO AND LOW-VARIANCE GENES
# Variance computed on log2 training data only, then the same
# mask is applied to both branches (SVM log2, RF/XGBoost untransformed).
# log2(1+x) is monotonic, so near-constant genes are the same
# on unstransformed or log2 data — one mask works for both.
#################################################################

train_gene_cols_log2 <- train_data_log2[, !colnames(train_data_log2) %in% c("Sample.ID", "Tissue.Type")]
test_gene_cols_log2  <- test_data_log2[,  !colnames(test_data_log2)  %in% c("Sample.ID", "Tissue.Type")]

# Per-gene variance, Training samples only
gene_vars <- matrixStats::colVars(as.matrix(train_gene_cols_log2))
names(gene_vars) <- colnames(train_gene_cols_log2)


var_threshold <- 0.01 

var_hist_df <- data.frame(variance = gene_vars[gene_vars > 0])

# Zero-variance genes (removed below)
n_zero_var_train <- sum(gene_vars == 0)
cat("\n--- Zero-Variance Gene Check (Training) ---\n")
cat("Genes with exactly zero variance (training):", n_zero_var_train, "\n")
cat("Genes with variance in (0, 0.01] (low but non-zero):",
    sum(gene_vars > 0 & gene_vars <= var_threshold), "\n")

# Variance histogram to justify the threshold (positive variances only, for the log scale)
variance_histogram <- ggplot(var_hist_df, aes(x = variance)) +
  geom_histogram(bins = 100, fill = "#0072B2", colour = "white", linewidth = 0.1) +
  geom_vline(xintercept = var_threshold,
             linetype = "dashed", colour = "grey30", linewidth = 0.5) +
  annotate("text", x = var_threshold, y = Inf,
           label = paste0("Threshold = ", var_threshold),
           angle = 90, vjust = -0.4, hjust = 1, size = 4, colour = "grey30") +
  scale_x_log10(
    breaks = c(0.00001, 0.0001, 0.001, 0.01, 0.02, 0.1, 1, 10),
    labels = c("0.00001", "0.0001", "0.001", "0.01", "0.02", "0.1", "1", "10"),
    minor_breaks = NULL
  ) +
  labs(
    title    = "Distribution of Per-Gene Variance (Training Set)",
    subtitle = "A variance threshold (dashed line) removes the low-variance shoulder; the informative peak is retained.",
    x        = "Per-Gene Variance (log10 scale)",
    y        = "Number of Genes",
    caption  = paste0("Variance computed on log2(x+1) RPKM gene expression, ", nrow(train_gene_cols_log2),
                      " training samples across ", ncol(train_gene_cols_log2), " genes.")
  ) +
  theme_classic(base_size = 12, base_family = "Arial") +
  theme(
    plot.title    = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 12, colour = "grey40"),
    plot.caption  = element_text(size = 12, colour = "grey50")
  )

ggsave(out("variance_histogram_training.png"), plot = variance_histogram,
       width = 10, height = 6, dpi = 300, bg = "white")
cat("Variance distribution histogram saved\n")

# Variance quantiles (to check the threshold is defensible)
cat("\n--- Per-gene variance quantiles (training) ---\n")
print(round(quantile(gene_vars, probs = c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1)), 5))

# Apply the training-derived mask to both sets
 low_var_mask <- gene_vars <= var_threshold | is.na(gene_vars)
cat("\n--- Low-Variance Gene Removal Summary ---\n")
cat("Variance threshold used:", var_threshold, "\n")
cat("Genes removed (zero or low variance):", sum(low_var_mask), "\n")
cat("Genes remaining:", sum(!low_var_mask),
    "(", round(sum(!low_var_mask) / length(low_var_mask) * 100, 1), "% kept )\n")
train_gene_cols_log2 <- train_gene_cols_log2[, !low_var_mask]
test_gene_cols_log2  <- test_gene_cols_log2[,  !low_var_mask]

cat("Genes remaining in training set:", ncol(train_gene_cols_log2), "\n")
cat("Genes remaining in test set:",     ncol(test_gene_cols_log2), "\n")
cat("Both sets have the same number of genes:",
    ncol(train_gene_cols_log2) == ncol(test_gene_cols_log2), "\n")

################################################
# STEP 17: BUILD MATRICES FOR MODELLING
# Two unscaled versions built here: log2(1+x) and untransformed.
# Neither is standardised. Standardisation is deferred to inside
# nested CV, where each fold refits its own scaling, to avoid the
# leakage of computing scaling parameters before the folds exist.
################################################

# Log2 branch — used by Elastic Net and SVM (both scale/skew sensitive)
x_train_log2 <- as.matrix(train_gene_cols_log2)
x_test_log2  <- as.matrix(test_gene_cols_log2)

# Untransformed branch. Random Forest and XGBoost only. Reuses the same
# train_index (Step 15) and low_var_mask (Step 16) on the untransformed data.
train_gene_cols_untransformed <- gene_cols[train_index, ][, !low_var_mask]
test_gene_cols_untransformed  <- gene_cols[-train_index, ][, !low_var_mask]
x_train_untransformed <- as.matrix(train_gene_cols_untransformed)
x_test_untransformed  <- as.matrix(test_gene_cols_untransformed)

# Labels as separate vectors (the modelling functions expect x and y apart).
# Only one version needed. Tissue.Type was never transformed.
y_train <- train_data_log2$Tissue.Type
y_test  <- test_data_log2$Tissue.Type

# Attach Sample.ID as row names for traceability (e.g. SHAP waterfall plots),
# without exposing it to the models as a predictor.
rownames(x_train_log2)          <- train_data_log2$Sample.ID
rownames(x_test_log2)           <- test_data_log2$Sample.ID
rownames(x_train_untransformed) <- train_data_log2$Sample.ID
rownames(x_test_untransformed)  <- test_data_log2$Sample.ID

cat("\nx_train_log2 (EN/SVM):", dim(x_train_log2), "\n")
cat("x_train_untransformed (RF):", dim(x_train_untransformed), "\n")
cat("Same gene count in both branches:",
    ncol(x_train_log2) == ncol(x_train_untransformed), "\n")

################################################
# STEP 18: SHARED OUTER FOLDS
# One fixed 5-fold split of the training samples, stratified by
# Tissue.Type, reused across all three models so each is evaluated
# on identical held-out samples per fold. Needed for fair comparison.
################################################

set.seed(42)

outer_folds <- caret::createFolds(
  y_train,
  k = 5,
  list = TRUE,
  returnTrain = FALSE
)

cat("\nShared outer folds created\n")
cat("Number of folds:", length(outer_folds), "\n")
cat("Samples per fold:", sapply(outer_folds, length), "\n")

#################################################################
# STEP 19A: ALPHA DIAGNOSTIC
# Run once on full training set to determine CV-optimal alpha.
# Shows CV error and gene count at each alpha value.
# Evidence for alpha=0.8 selection used in en_filter (Step 19B).
##################################################################

cat("\n--- Step 19A: Elastic net alpha diagnostic ---\n")

set.seed(1000)

alpha_diagnostic <- do.call(rbind, lapply(seq(0.1, 1, by = 0.1), function(a) {
  fit <- cv.glmnet(
    log2(x_train_untransformed + 1),
    y_train,
    family = "binomial",
    alpha  = a,
    nfolds = 5
  )
  n_genes <- sum(coef(fit, s = "lambda.1se") != 0) - 1
  data.frame(
    Alpha    = a,
    N_Genes  = n_genes,
    CV_Error = round(fit$cvm[fit$index["1se", ]], 5)
  )
}))

cat("Alpha | Genes | CV Error\n")
print(alpha_diagnostic)

best_alpha <- alpha_diagnostic$Alpha[which.min(alpha_diagnostic$CV_Error)]
cat("\nCV-optimal alpha:", best_alpha, "\n")
cat("Genes at alpha=0.8:",
    alpha_diagnostic$N_Genes[which.min(abs(alpha_diagnostic$Alpha - 0.8))], "\n")

write.csv(alpha_diagnostic,
          out("elastic_net_alpha_diagnostic.csv"),
          row.names = FALSE)
cat("Diagnostic saved to elastic_net_alpha_diagnostic.csv\n")

################################################
# STEP 19B: SHARED ELASTIC NET FILTER FUNCTION
# Called fresh inside every fold, for all three classifiers.
# Elastic Net is used only as a feature filter (not a classifier)
# to reduce the gene count going into each model.
#
# already_log tells the function whether x is already log2-transformed,
# so it always fits on log2 exactly once:
#   - RF / XGBoost pass raw data  -> already_log = FALSE -> log2 internally
#   - SVM passes log2 data        -> already_log = TRUE  -> used as is
# This prevents a double log2 transform on the SVM branch.
#
# Alpha fixed at 0.8 (Step 19A diagnostic on full training set):
# alpha=0.2 was CV-optimal (384 genes, CV error 0.02571) but 0.8 was
# chosen for interpretability (64 genes, CV error 0.02856; Δ=0.00285,
# negligible). set.seed(1000) gives identical selection for the
# raw and log2 calls on the same fold.
################################################

en_filter <- function(y, x, already_log = FALSE, ...) {
  
  x_log <- if (already_log) x else log2(x + 1)
  
# Fixed seed ensures identical gene selection across model
# families for the same fold's data
  set.seed(1000)
  
  cvfit <- cv.glmnet(x_log, y,
                     family = "binomial",
                     alpha  = 0.8,
                     nfolds = 5)
  
  sel <- which(coef(cvfit, s = "lambda.1se")[-1, ] != 0)
  if (length(sel) == 0) {
    sel <- which(coef(cvfit, s = "lambda.min")[-1, ] != 0)
  }
  sel
}

cat("\nElastic Net filter function defined (alpha=0.8)\n")
cat("Expected gene count per fold: ~35-55 genes\n")
cat("Expected final model gene count: ~64 genes\n")

#####################################################################
# STEP 19C: PER-FOLD GENE SELECTION (SHARED ACROSS MODELS)

# The elastic-net filter is deterministic (set.seed(1000)) and uses
# the shared outer_folds, so gene selection is identical for all
# three models. Computed once here, before the models.
######################################################################

# Run the filter on each outer fold's training portion
per_fold_genes <- lapply(seq_along(outer_folds), function(f) {
  train_idx <- setdiff(seq_len(nrow(x_train_untransformed)), outer_folds[[f]])
  sel <- en_filter(y = y_train[train_idx],
                   x = x_train_untransformed[train_idx, ],
                   already_log = FALSE)
  colnames(x_train_untransformed)[sel]
})
names(per_fold_genes) <- paste0("Fold", seq_along(per_fold_genes))

# CSV 1: gene names listed under each fold column, with a per-fold total row
max_len   <- max(sapply(per_fold_genes, length))
fold_cols <- lapply(per_fold_genes, function(g) c(sort(g), rep("", max_len - length(g))))
fold_df   <- as.data.frame(fold_cols, stringsAsFactors = FALSE, check.names = FALSE)
totals    <- sapply(per_fold_genes, length)
fold_df   <- rbind(fold_df, as.character(totals))
fold_df   <- cbind(Row = c(seq_len(max_len), "TOTAL_GENES"), fold_df)
write.csv(fold_df, out("genes_per_fold.csv"), row.names = FALSE)
cat("Genes per fold:", totals, "\n")

# CSV 2: final gene panel (filter on the full training set)
final_sel   <- en_filter(y = y_train, x = x_train_untransformed, already_log = FALSE)
final_genes <- sort(colnames(x_train_untransformed)[final_sel])
write.csv(data.frame(Gene = final_genes), out("final_gene_panel.csv"), row.names = FALSE)
cat("Final panel:", length(final_genes), "genes\n")

# CSV 3: genes selected in all five folds 
stable_genes <- Reduce(intersect, per_fold_genes)
write.csv(data.frame(Gene = stable_genes), out("genes_in_all_folds.csv"), row.names = FALSE)
cat("Genes in all", length(outer_folds), "folds:", length(stable_genes), "\n")

################################################
# STEP 20: XGBOOST NESTED CROSS-VALIDATION
# Classifier 1 of 3. Trees are scale-invariant, so XGBoost uses the
# untransformed branch. The elastic-net filter is refit inside every
# fold (already_log = FALSE); shared outer_folds give all three models
# identical held-out samples. Tuned on AUC.
################################################

xgb_grid <- expand.grid(
  nrounds          = c(50, 100, 150),
  max_depth        = c(2, 3),
  eta              = c(0.1, 0.3),
  gamma            = 0,
  colsample_bytree = c(0.6, 0.8),
  min_child_weight = 1,
  subsample        = 0.8
)
xgb_ctrl <- caret::trainControl(method = "cv", number = 5,
                                classProbs = TRUE, summaryFunction = caret::twoClassSummary)
set.seed(42)
xgb_fit <- nestcv.train(
  y = y_train, x = x_train_untransformed, method = "xgbTree",
  filterFUN = en_filter, filter_options = list(already_log = FALSE),
  outer_folds = outer_folds, n_inner_folds = 5,
  metric = "ROC", trControl = xgb_ctrl, tuneGrid = xgb_grid, cv.cores = 1
)

###################################################
# STEP 20B: XGBOOST RESULTS EXTRACTION 
###################################################

xgb_n_folds <- length(xgb_fit$outer_result)

# Per-fold AUC (from each fold's held-out predictions)
xgb_fold_aucs <- sapply(seq_len(xgb_n_folds), function(f) {
  fold <- xgb_fit$outer_result[[f]]
  as.numeric(pROC::auc(fold$preds$testy, fold$preds$predyp,
                       levels = c("Normal", "Psoriasis"), direction = "<"))
})
write.csv(data.frame(Fold = seq_len(xgb_n_folds), AUC = round(xgb_fold_aucs, 4)),
          out("xgb_per_fold_auc.csv"), row.names = FALSE)

# Pooled outer-fold predictions (nested CV, not resubstitution)
xgb_true_class <- factor(xgb_fit$output$testy, levels = c("Normal", "Psoriasis"))
xgb_pred_class <- factor(xgb_fit$output$predy, levels = c("Normal", "Psoriasis"))
xgb_prob_pos   <- xgb_fit$output$predyp
xgb_cm_cv      <- caret::confusionMatrix(xgb_pred_class, xgb_true_class, positive = "Psoriasis")
xgb_pooled_auc <- as.numeric(pROC::auc(xgb_true_class, xgb_prob_pos,
                                       levels = c("Normal", "Psoriasis"), direction = "<"))
xgb_pooled_metrics <- data.frame(
  Metric = c("AUC (pooled)", "Accuracy", "Sensitivity", "Specificity", "F1", "MCC"),
  Value  = c(round(xgb_pooled_auc, 4),
             round(as.numeric(xgb_cm_cv$overall["Accuracy"]), 4),
             round(as.numeric(xgb_cm_cv$byClass["Sensitivity"]), 4),
             round(as.numeric(xgb_cm_cv$byClass["Specificity"]), 4),
             round(as.numeric(xgb_cm_cv$byClass["F1"]), 4),
             round(mltools::mcc(preds = xgb_pred_class, actuals = xgb_true_class), 4))
)

write.csv(xgb_fit$final_fit$bestTune, out("xgb_final_hyperparameters.csv"), row.names = FALSE)

############################################################
# STEP 21: RANDOM FOREST NESTED CROSS-VALIDATION
# Classifier 2 of 3. untransformed branch (scale-invariant). 
############################################################

rf_grid <- expand.grid(mtry = c(5, 10, 20, 30), splitrule = "gini", min.node.size = c(5, 10, 15))
rf_ctrl <- caret::trainControl(method = "cv", number = 5,
                               classProbs = TRUE, summaryFunction = caret::twoClassSummary)
set.seed(42)

rf_fit <- nestcv.train(
  y = y_train, x = x_train_untransformed, method = "ranger",
  filterFUN = en_filter, filter_options = list(already_log = FALSE),
  outer_folds = outer_folds, n_inner_folds = 5,
  metric = "ROC", trControl = rf_ctrl, tuneGrid = rf_grid,
  num.trees = 1000, importance = "impurity", cv.cores = 1
)

################################################
# STEP 21B: RANDOM FOREST RESULTS EXTRACTION
################################################

rf_n_folds <- length(rf_fit$outer_result)
rf_fold_aucs <- sapply(seq_len(rf_n_folds), function(f) {
  fold <- rf_fit$outer_result[[f]]
  as.numeric(pROC::auc(fold$preds$testy, fold$preds$predyp,
                       levels = c("Normal", "Psoriasis"), direction = "<"))
})
write.csv(data.frame(Fold = seq_len(rf_n_folds), AUC = round(rf_fold_aucs, 4)),
          out("rf_per_fold_auc.csv"), row.names = FALSE)

rf_true_class <- factor(rf_fit$output$testy, levels = c("Normal", "Psoriasis"))
rf_pred_class <- factor(rf_fit$output$predy, levels = c("Normal", "Psoriasis"))
rf_prob_pos   <- rf_fit$output$predyp
rf_cm_cv      <- caret::confusionMatrix(rf_pred_class, rf_true_class, positive = "Psoriasis")
rf_pooled_auc <- as.numeric(pROC::auc(rf_true_class, rf_prob_pos,
                                      levels = c("Normal", "Psoriasis"), direction = "<"))
rf_pooled_metrics <- data.frame(
  Metric = c("AUC (pooled)", "Accuracy", "Sensitivity", "Specificity", "F1", "MCC"),
  Value  = c(round(rf_pooled_auc, 4),
             round(as.numeric(rf_cm_cv$overall["Accuracy"]), 4),
             round(as.numeric(rf_cm_cv$byClass["Sensitivity"]), 4),
             round(as.numeric(rf_cm_cv$byClass["Specificity"]), 4),
             round(as.numeric(rf_cm_cv$byClass["F1"]), 4),
             round(mltools::mcc(preds = rf_pred_class, actuals = rf_true_class), 4))
)

write.csv(rf_fit$final_fit$bestTune, out("rf_final_hyperparameters.csv"), row.names = FALSE)

########################################################
# STEP 22: SVM (LINEAR KERNEL) NESTED CROSS-VALIDATION
# Classifier 3 of 3. 
########################################################

svm_grid <- expand.grid(C = c(0.25, 0.5, 1, 2, 4))
svm_ctrl <- caret::trainControl(method = "cv", number = 5,
                                classProbs = TRUE, summaryFunction = caret::twoClassSummary)
set.seed(42)

svm_fit <- nestcv.train(
  y = y_train, x = x_train_log2, method = "svmLinear",
  filterFUN = en_filter, filter_options = list(already_log = TRUE),
  outer_folds = outer_folds, n_inner_folds = 5,
  metric = "ROC", trControl = svm_ctrl, tuneGrid = svm_grid,
  preProcess = c("center", "scale"), cv.cores = 1
)

################################################
# STEP 22B: SVM RESULTS EXTRACTION
################################################

svm_n_folds <- length(svm_fit$outer_result)
svm_fold_aucs <- sapply(seq_len(svm_n_folds), function(f) {
  fold <- svm_fit$outer_result[[f]]
  as.numeric(pROC::auc(fold$preds$testy, fold$preds$predyp,
                       levels = c("Normal", "Psoriasis"), direction = "<"))
})
write.csv(data.frame(Fold = seq_len(svm_n_folds), AUC = round(svm_fold_aucs, 4)),
          out("svm_per_fold_auc.csv"), row.names = FALSE)

svm_true_class <- factor(svm_fit$output$testy, levels = c("Normal", "Psoriasis"))
svm_pred_class <- factor(svm_fit$output$predy, levels = c("Normal", "Psoriasis"))
svm_prob_pos   <- svm_fit$output$predyp
svm_cm_cv      <- caret::confusionMatrix(svm_pred_class, svm_true_class, positive = "Psoriasis")
svm_pooled_auc <- as.numeric(pROC::auc(svm_true_class, svm_prob_pos,
                                       levels = c("Normal", "Psoriasis"), direction = "<"))
svm_pooled_metrics <- data.frame(
  Metric = c("AUC (pooled)", "Accuracy", "Sensitivity", "Specificity", "F1", "MCC"),
  Value  = c(round(svm_pooled_auc, 4),
             round(as.numeric(svm_cm_cv$overall["Accuracy"]), 4),
             round(as.numeric(svm_cm_cv$byClass["Sensitivity"]), 4),
             round(as.numeric(svm_cm_cv$byClass["Specificity"]), 4),
             round(as.numeric(svm_cm_cv$byClass["F1"]), 4),
             round(mltools::mcc(preds = svm_pred_class, actuals = svm_true_class), 4))
)

write.csv(svm_fit$final_fit$bestTune, out("svm_final_hyperparameters.csv"), row.names = FALSE)

###########################################################################
# STEP 23: TRAINING-SET EVALUATION (NESTED CROSS-VALIDATED)
# Uses the pooled outer-fold held-out predictions (nested CV). 
# Produces, for all three models: confusion matrices,
# a combined pooled-metric table with mean +/- SD outer-fold AUC, and an
# overlaid ROC plot with bootstrap CIs.
############################################################################

cv_models <- list(
  XGBoost         = list(true = xgb_true_class, prob = xgb_prob_pos, cm = xgb_cm_cv,
                         fold_aucs = xgb_fold_aucs, metrics = xgb_pooled_metrics),
  `Random Forest` = list(true = rf_true_class,  prob = rf_prob_pos,  cm = rf_cm_cv,
                         fold_aucs = rf_fold_aucs,  metrics = rf_pooled_metrics),
  SVM             = list(true = svm_true_class, prob = svm_prob_pos, cm = svm_cm_cv,
                         fold_aucs = svm_fold_aucs, metrics = svm_pooled_metrics)
)

# Confusion matrices (long format)
cv_cm_long <- do.call(rbind, lapply(names(cv_models), function(m)
  transform(as.data.frame(cv_models[[m]]$cm$table), Model = m)))
write.csv(cv_cm_long, out("cv_confusion_matrices.csv"), row.names = FALSE)

# Combined pooled metrics + mean +/- SD outer-fold AUC (one row per model)
cv_metrics_all <- do.call(rbind, lapply(names(cv_models), function(m) {
  mt <- cv_models[[m]]$metrics
  get_v <- function(name) mt$Value[mt$Metric == name]
  data.frame(
    Model         = m,
    Mean_fold_AUC = round(mean(cv_models[[m]]$fold_aucs), 4),
    SD_fold_AUC   = round(sd(cv_models[[m]]$fold_aucs), 4),
    Pooled_AUC    = get_v("AUC (pooled)"),
    Accuracy      = get_v("Accuracy"),
    Sensitivity   = get_v("Sensitivity"),
    Specificity   = get_v("Specificity"),
    F1            = get_v("F1"),
    MCC           = get_v("MCC")
  )
}))

write.csv(cv_metrics_all, out("cv_pooled_metrics_all_models.csv"), row.names = FALSE)

# ROC objects + bootstrap CIs (used for the plot legend)
cv_roc <- lapply(names(cv_models), function(m) {
  d <- cv_models[[m]]
  roc_obj <- pROC::roc(d$true, d$prob, levels = c("Normal", "Psoriasis"), direction = "<")
  set.seed(42)
  ci_val <- suppressMessages(as.numeric(pROC::ci.auc(roc_obj, method = "bootstrap", B = 2000)))
  list(model = m, roc = roc_obj, auc = as.numeric(pROC::auc(roc_obj)),
       ci_low = ci_val[1], ci_high = ci_val[3])
})
names(cv_roc) <- names(cv_models)

# Overlaid CV ROC curves
png(out("cv_roc_overlay.png"), width = 1800, height = 1600, res = 300)
plot(cv_roc[["XGBoost"]]$roc, col = "#1b9e77", lwd = 2.5, lty = 1, legacy.axes = TRUE,
     main = "Nested cross-validated ROC curves (training set)",
     xlab = "False positive rate (1 - Specificity)", ylab = "True positive rate (Sensitivity)")
plot(cv_roc[["Random Forest"]]$roc, col = "#d95f02", lwd = 2.5, lty = 2, add = TRUE)
plot(cv_roc[["SVM"]]$roc,           col = "#7570b3", lwd = 2.5, lty = 3, add = TRUE)
legend("bottomright", bty = "n", lwd = c(2.5, 2.5, 2.5, 1), lty = c(1, 2, 3, 1),
       col = c("#1b9e77", "#d95f02", "#7570b3", "grey"),
       legend = c(
         sprintf("XGBoost       (AUC %.3f, 95%% CI %.3f\u2013%.3f)", cv_roc[["XGBoost"]]$auc, cv_roc[["XGBoost"]]$ci_low, cv_roc[["XGBoost"]]$ci_high),
         sprintf("Random Forest (AUC %.3f, 95%% CI %.3f\u2013%.3f)", cv_roc[["Random Forest"]]$auc, cv_roc[["Random Forest"]]$ci_low, cv_roc[["Random Forest"]]$ci_high),
         sprintf("SVM           (AUC %.3f, 95%% CI %.3f\u2013%.3f)", cv_roc[["SVM"]]$auc, cv_roc[["SVM"]]$ci_low, cv_roc[["SVM"]]$ci_high),
         "Chance (AUC = 0.5)"))
dev.off()
cat("Training (nested CV) evaluation saved\n")

###########################################################################
# STEP 24: TEST-SET EVALUATION (ALL 3 MODELS)
# Final models scored on the held-out test set, untouched by selection,
# tuning, or fitting. Confusion matrices, pooled metrics, mean AUC,
# overlaid ROC with bootstrap CIs.
###########################################################################

test_set_metrics <- function(fit, x_test, y_test, model_name) {
  pred_class <- factor(predict(fit, newdata = x_test), levels = c("Normal", "Psoriasis"))
  prob_pos   <- predict(fit, newdata = x_test, type = "prob")[, "Psoriasis"]
  cm  <- caret::confusionMatrix(pred_class, y_test, positive = "Psoriasis")
  roc_obj <- pROC::roc(y_test, prob_pos, levels = c("Normal", "Psoriasis"), direction = "<")
  set.seed(42)
  ci_val <- suppressMessages(as.numeric(pROC::ci.auc(roc_obj, method = "bootstrap", B = 2000)))
  metrics <- data.frame(
    Model = model_name, AUC = round(as.numeric(pROC::auc(roc_obj)), 4),
    AUC_CI_low = round(ci_val[1], 4), AUC_CI_high = round(ci_val[3], 4),
    Accuracy = round(as.numeric(cm$overall["Accuracy"]), 4),
    Sensitivity = round(as.numeric(cm$byClass["Sensitivity"]), 4),
    Specificity = round(as.numeric(cm$byClass["Specificity"]), 4),
    F1 = round(as.numeric(cm$byClass["F1"]), 4),
    MCC = round(mltools::mcc(preds = pred_class, actuals = y_test), 4))
  list(metrics = metrics, cm = cm, roc = roc_obj, model_name = model_name)
}

y_test <- factor(y_test, levels = c("Normal", "Psoriasis"))
xgb_test <- test_set_metrics(xgb_fit, x_test_untransformed, y_test, "XGBoost")
rf_test  <- test_set_metrics(rf_fit,  x_test_untransformed, y_test, "Random Forest")
svm_test <- test_set_metrics(svm_fit, x_test_log2,          y_test, "SVM")

test_metrics_all <- rbind(xgb_test$metrics, rf_test$metrics, svm_test$metrics)
write.csv(test_metrics_all, out("test_set_metrics_all_models.csv"), row.names = FALSE)

test_cm_long <- do.call(rbind, lapply(list(xgb_test, rf_test, svm_test), function(t)
  transform(as.data.frame(t$cm$table), Model = t$model_name)))
write.csv(test_cm_long, out("test_set_confusion_matrices.csv"), row.names = FALSE)

png(out("test_set_roc_overlay.png"), width = 1800, height = 1600, res = 300)
plot(xgb_test$roc, col = "#1b9e77", lwd = 2.5, lty = 1, legacy.axes = TRUE,
     main = "Test-set ROC curves", xlab = "False positive rate (1 - Specificity)",
     ylab = "True positive rate (Sensitivity)")
plot(rf_test$roc,  col = "#d95f02", lwd = 2.5, lty = 2, add = TRUE)
plot(svm_test$roc, col = "#7570b3", lwd = 2.5, lty = 3, add = TRUE)
legend("bottomright", bty = "n", lwd = c(2.5, 2.5, 2.5, 1), lty = c(1, 2, 3, 1),
       col = c("#1b9e77", "#d95f02", "#7570b3", "grey"),
       legend = c(
         sprintf("XGBoost       (AUC %.3f, 95%% CI %.3f\u2013%.3f)", xgb_test$metrics$AUC, xgb_test$metrics$AUC_CI_low, xgb_test$metrics$AUC_CI_high),
         sprintf("Random Forest (AUC %.3f, 95%% CI %.3f\u2013%.3f)", rf_test$metrics$AUC, rf_test$metrics$AUC_CI_low, rf_test$metrics$AUC_CI_high),
         sprintf("SVM           (AUC %.3f, 95%% CI %.3f\u2013%.3f)", svm_test$metrics$AUC, svm_test$metrics$AUC_CI_low, svm_test$metrics$AUC_CI_high),
         "Chance (AUC = 0.5)"))
dev.off()

############################################################################
# STEP 25: CV vs TEST COMPARISON (OVERFITTING CHECK)
# A small or zero gap between nested-CV and test performance confirms no
# meaningful overfitting.
############################################################################

get_metric <- function(df, m) as.numeric(df$Value[df$Metric == m])
cv_vs_test <- data.frame(
  Model = c("XGBoost", "Random Forest", "SVM"),
  CV_AUC_mean = round(c(mean(xgb_fold_aucs), mean(rf_fold_aucs), mean(svm_fold_aucs)), 4),
  CV_AUC_SD   = round(c(sd(xgb_fold_aucs),   sd(rf_fold_aucs),   sd(svm_fold_aucs)),   4),
  CV_Accuracy    = c(get_metric(xgb_pooled_metrics,"Accuracy"),    get_metric(rf_pooled_metrics,"Accuracy"),    get_metric(svm_pooled_metrics,"Accuracy")),
  CV_Sensitivity = c(get_metric(xgb_pooled_metrics,"Sensitivity"), get_metric(rf_pooled_metrics,"Sensitivity"), get_metric(svm_pooled_metrics,"Sensitivity")),
  CV_Specificity = c(get_metric(xgb_pooled_metrics,"Specificity"), get_metric(rf_pooled_metrics,"Specificity"), get_metric(svm_pooled_metrics,"Specificity")),
  CV_F1  = c(get_metric(xgb_pooled_metrics,"F1"),  get_metric(rf_pooled_metrics,"F1"),  get_metric(svm_pooled_metrics,"F1")),
  CV_MCC = c(get_metric(xgb_pooled_metrics,"MCC"), get_metric(rf_pooled_metrics,"MCC"), get_metric(svm_pooled_metrics,"MCC")),
  Test_AUC         = c(xgb_test$metrics$AUC,         rf_test$metrics$AUC,         svm_test$metrics$AUC),
  Test_Accuracy    = c(xgb_test$metrics$Accuracy,    rf_test$metrics$Accuracy,    svm_test$metrics$Accuracy),
  Test_Sensitivity = c(xgb_test$metrics$Sensitivity, rf_test$metrics$Sensitivity, svm_test$metrics$Sensitivity),
  Test_Specificity = c(xgb_test$metrics$Specificity, rf_test$metrics$Specificity, svm_test$metrics$Specificity),
  Test_F1          = c(xgb_test$metrics$F1,          rf_test$metrics$F1,          svm_test$metrics$F1),
  Test_MCC         = c(xgb_test$metrics$MCC,         rf_test$metrics$MCC,         svm_test$metrics$MCC)
)
write.csv(cv_vs_test, out("cv_vs_test_comparison.csv"), row.names = FALSE)

#################################################################
# STEP 26: SHAP ANALYSIS — XGBOOST
# Exact Tree SHAP (fast for shallow trees).
# Outputs: beeswarm (top 10), bar chart (top 10), top-10 CSV.
#################################################################

dissertation_theme <- function() {
  theme_classic(base_size = 12, base_family = "Arial") +
    theme(plot.title = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(size = 12, colour = "grey40"),
          plot.caption = element_text(size = 12, colour = "grey50"),
          legend.title = element_text(face = "bold"),
          panel.grid.major.y = element_line(colour = "grey90"),
          axis.title = element_text(size = 12))
}
replace_colour_scale <- function(p, ...) {
  p$scales$scales <- Filter(function(s) {
    !inherits(s, "ScaleContinuous") ||
      (!("colour" %in% s$aesthetics) && !("color" %in% s$aesthetics))
  }, p$scales$scales)
  p + scale_color_gradient2(...) + scale_colour_gradient2(...)
}

xgb_model_raw   <- xgb_fit$final_fit$finalModel
xgb_test_matrix <- x_test_untransformed[, xgb_fit$final_vars]
xgb_dmatrix     <- xgboost::xgb.DMatrix(data = xgb_test_matrix)
n_final_genes   <- ncol(xgb_test_matrix)
shp_xgb <- shapviz(xgb_model_raw, X_pred = xgb_dmatrix, X = xgb_test_matrix)

# Direction check: SHAP row-sums should correlate positively with Psoriasis probability
xgb_test_probs <- predict(xgb_fit$final_fit, newdata = xgb_test_matrix, type = "prob")[, "Psoriasis"]
shap_prob_cor <- cor(rowSums(shp_xgb$S), xgb_test_probs)
if (shap_prob_cor < 0) shp_xgb$S <- -shp_xgb$S
cat("XGBoost SHAP direction (cor):", round(shap_prob_cor, 3), "\n")

mean_abs_shap    <- colMeans(abs(shp_xgb$S))
shap_ranked      <- sort(mean_abs_shap, decreasing = TRUE)
xgb_nonzero_shap <- sum(mean_abs_shap > 0)
write.csv(data.frame(Rank = 1:10, Gene = names(shap_ranked)[1:10],
                     Mean_Abs_SHAP = round(shap_ranked[1:10], 6)),
          out("xgb_shap_top10_genes.csv"), row.names = FALSE)

p_bee <- replace_colour_scale(
  sv_importance(shp_xgb, kind = "beeswarm", max_display = 10) +
    labs(title = "XGBoost SHAP Feature Importance \u2014 Top 10 Genes",
         subtitle = paste0("Test set (n = ", nrow(xgb_test_matrix),
                           "). Positive SHAP = towards Psoriasis; negative = towards Normal."),
         x = "SHAP value (impact on model output)", y = "Gene",
         caption = paste0("Orange = high expression, blue = low. Only ", xgb_nonzero_shap,
                          " of ", n_final_genes, " genes had non-zero SHAP. Log-odds scale.")) +
    dissertation_theme(),
  low = "#0072B2", mid = "white", high = "#E69F00", midpoint = 0.5, name = "Gene\nexpression\nlevel")
ggsave(out("xgb_shap_beeswarm_top10.png"), p_bee, width = 10, height = 7, dpi = 300, bg = "white")

shap_bar_df <- data.frame(Gene = factor(names(shap_ranked)[1:10], levels = rev(names(shap_ranked)[1:10])),
                          Mean_Abs_SHAP = as.numeric(shap_ranked[1:10]))
p_bar <- ggplot(shap_bar_df, aes(x = Mean_Abs_SHAP, y = Gene)) +
  geom_col(fill = "#0072B2") +
  geom_text(aes(label = sprintf("%.3f", Mean_Abs_SHAP), x = Mean_Abs_SHAP + max(Mean_Abs_SHAP) * 0.02),
            fontface = "bold", hjust = 0, size = 3.5, colour = "grey30", family = "Arial") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "XGBoost SHAP Mean Absolute Importance \u2014 Top 10 Genes",
       subtitle = paste0("Test set (n = ", nrow(xgb_test_matrix), "). Mean absolute SHAP across all test samples."),
       x = "Mean |SHAP value|", y = "Gene",
       caption = paste0("Only ", xgb_nonzero_shap, " of ", n_final_genes, " genes had non-zero weight. Log-odds scale.")) +
  dissertation_theme() +
  theme(panel.grid.major.y = element_blank(),
        plot.caption = element_text(hjust = 0.5, size = 12, colour = "grey50", family = "Arial"),
        plot.caption.position = "plot")
ggsave(out("xgb_shap_bar_top10.png"), p_bar, width = 9, height = 6, dpi = 300, bg = "white")

######################################################################################
# STEP 27: SHAP ANALYSIS — RANDOM FOREST
# Exact Tree SHAP via treeshap. Compatibility patch: ranger 0.18.0 names
# columns pred.Psoriasis/pred.Normal in treeInfo() but treeshap 0.4.0
# expects pred.1 — the patch adds pred.1, runs unify(), then restores the
# original. Raw branch. Outputs: beeswarm, bar chart, top-10 + all-nonzero CSV.
#######################################################################################
rf_model_raw    <- rf_fit$final_fit$finalModel
rf_test_matrix  <- as.data.frame(x_test_untransformed[, rf_fit$final_vars])
rf_train_matrix <- as.data.frame(x_train_untransformed[, rf_fit$final_vars])

original_treeInfo <- ranger::treeInfo
patched_treeInfo <- function(object, tree = 1) {
  td <- original_treeInfo(object, tree)
  if ("pred.Psoriasis" %in% names(td) && !"pred.1" %in% names(td)) td$pred.1 <- td$pred.Psoriasis
  td
}
assignInNamespace("treeInfo", patched_treeInfo, ns = "ranger")
rf_unified <- unify(rf_model_raw, rf_train_matrix)
assignInNamespace("treeInfo", original_treeInfo, ns = "ranger")

rf_treeshap <- treeshap(rf_unified, rf_test_matrix, verbose = FALSE)
shp_rf <- shapviz(rf_treeshap, X = rf_test_matrix)

# Per-gene direction check against known-upregulated psoriasis genes
check_genes <- intersect(c("S100A12", "IL1F6", "IL1F9", "SERPINB4"), colnames(shp_rf$S))
for (g in check_genes) {
  r <- cor(rf_test_matrix[[g]], shp_rf$S[, g])
  cat(sprintf("RF %-10s cor(expr, SHAP) = %+.3f %s\n", g, r, ifelse(r > 0, "OK", "*** INVERTED ***")))
}

rf_mean_abs_shap <- colMeans(abs(shp_rf$S))
rf_shap_ranked   <- sort(rf_mean_abs_shap, decreasing = TRUE)
rf_nonzero_shap  <- sum(rf_mean_abs_shap > 0)
write.csv(data.frame(Rank = 1:10, Gene = names(rf_shap_ranked)[1:10],
                     Mean_Abs_SHAP = round(rf_shap_ranked[1:10], 6)),
          out("rf_shap_top10_genes.csv"), row.names = FALSE)
write.csv(data.frame(Rank = seq_len(rf_nonzero_shap), Gene = names(rf_shap_ranked)[1:rf_nonzero_shap],
                     Mean_Abs_SHAP = round(rf_shap_ranked[1:rf_nonzero_shap], 6)),
          out("rf_shap_all_nonzero_genes.csv"), row.names = FALSE)

p_bee_rf <- replace_colour_scale(
  sv_importance(shp_rf, kind = "beeswarm", max_display = 10) +
    labs(title = "Random Forest SHAP Feature Importance \u2014 Top 10 Genes",
         subtitle = paste0("Test set (n = ", nrow(rf_test_matrix),
                           "). Positive SHAP = towards Psoriasis; negative = towards Normal."),
         x = "SHAP value (impact on model output)", y = "Gene",
         caption = paste0("Orange = high expression, blue = low. Only ", rf_nonzero_shap,
                          " of ", n_final_genes, " genes had non-zero SHAP. Probability scale.")) +
    dissertation_theme() +
    theme(plot.caption = element_text(hjust = 0.5, size = 12, colour = "grey50", family = "Arial"),
          plot.caption.position = "plot"),
  low = "#0072B2", mid = "white", high = "#E69F00", midpoint = 0.5, name = "Gene\nexpression\nlevel")
ggsave(out("rf_shap_beeswarm_top10.png"), p_bee_rf, width = 10, height = 7, dpi = 300, bg = "white")

rf_shap_bar_df <- data.frame(Gene = factor(names(rf_shap_ranked)[1:10], levels = rev(names(rf_shap_ranked)[1:10])),
                             Mean_Abs_SHAP = as.numeric(rf_shap_ranked[1:10]))
p_bar_rf <- ggplot(rf_shap_bar_df, aes(x = Mean_Abs_SHAP, y = Gene)) +
  geom_col(fill = "#0072B2") +
  geom_text(aes(label = sprintf("%.3f", Mean_Abs_SHAP), x = Mean_Abs_SHAP + max(Mean_Abs_SHAP) * 0.02),
            fontface = "bold", hjust = 0, size = 3.5, colour = "grey30", family = "Arial") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Random Forest SHAP Mean Absolute Importance \u2014 Top 10 Genes",
       subtitle = paste0("Test set (n = ", nrow(rf_test_matrix), "). Mean absolute SHAP across all test samples."),
       x = "Mean |SHAP value|", y = "Gene",
       caption = paste0(rf_nonzero_shap, " of ", n_final_genes, " genes had non-zero weight. Probability scale.")) +
  dissertation_theme() +
  theme(panel.grid.major.y = element_blank(),
        plot.caption = element_text(hjust = 0.5, size = 12, colour = "grey50", family = "Arial"),
        plot.caption.position = "plot")
ggsave(out("rf_shap_bar_top10.png"), p_bar_rf, width = 9, height = 6, dpi = 300, bg = "white")

################################################
# STEP 28: SHAP ANALYSIS — SVM (LINEAR KERNEL)
# Kernel SHAP (TreeSHAP does not apply to SVMs). It did not
# converge within tolerance (combinatorial cost over 64 features), so the
# convergence extent and per-gene direction are checked below; directional
# attribution is the basis for interpretation and rank ordering is
# indicative only. Probability scale. Log2 branch.
# Outputs: beeswarm, bar chart, top-10 + all-nonzero CSV, diagnostics.
################################################
svm_test_matrix  <- x_test_log2[,  svm_fit$final_vars]
svm_train_matrix <- x_train_log2[, svm_fit$final_vars]
pred_fun_svm <- function(object, newdata) predict(object, newdata = newdata, type = "prob")[, "Psoriasis"]

set.seed(42)
ks_svm <- kernelshap(svm_fit$final_fit, X = svm_test_matrix, bg_X = svm_train_matrix,
                     pred_fun = pred_fun_svm, hybrid_degree = 2)

# Convergence diagnostic: how many test samples (rows) converged
print(ks_svm)
svm_converged <- tryCatch(ks_svm$converged, error = function(e) NA)
if (!all(is.na(svm_converged))) {
  cat("KernelSHAP converged for", sum(svm_converged), "of",
      length(svm_converged), "test samples\n")
} else {
  cat("Convergence flag not exposed by this kernelshap version; see print(ks_svm) above\n")
}

shp_svm <- shapviz(ks_svm)

# Rank genes by mean |SHAP|
svm_mean_abs_shap <- colMeans(abs(shp_svm$S))
svm_shap_ranked   <- sort(svm_mean_abs_shap, decreasing = TRUE)
svm_nonzero_shap  <- sum(svm_mean_abs_shap > 0)
svm_top10 <- data.frame(Rank = 1:10, Gene = names(svm_shap_ranked)[1:10],
                        Mean_Abs_SHAP = round(svm_shap_ranked[1:10], 6))

# Direction check: expression-to-SHAP correlation per top-10 gene.
# Given non-convergence, a strong, consistent relationship is the evidence
# that the direction of each gene's effect is still reliable.
svm_dir_cor <- sapply(svm_top10$Gene, function(g) cor(svm_test_matrix[, g], shp_svm$S[, g]))
for (g in svm_top10$Gene)
  cat(sprintf("SVM %-12s cor = %+.3f (%s)\n", g, svm_dir_cor[g],
              ifelse(svm_dir_cor[g] > 0, "towards Psoriasis", "towards Normal")))
cat("Minimum |correlation| across top-10:", round(min(abs(svm_dir_cor)), 3), "\n")

# Save top-10 and all non-zero attributions
write.csv(svm_top10, out("svm_shap_top10_genes.csv"), row.names = FALSE)
write.csv(data.frame(Rank = seq_len(svm_nonzero_shap), Gene = names(svm_shap_ranked)[1:svm_nonzero_shap],
                     Mean_Abs_SHAP = round(svm_shap_ranked[1:svm_nonzero_shap], 6)),
          out("svm_shap_all_nonzero_genes.csv"), row.names = FALSE)

# (1) Beeswarm
p_bee_svm <- replace_colour_scale(
  sv_importance(shp_svm, kind = "beeswarm", max_display = 10) +
    labs(title = "Support Vector Machine SHAP Feature Importance \u2014 Top 10 Genes",
         subtitle = paste0("Test set (n = ", nrow(svm_test_matrix),
                           "). Positive SHAP = towards Psoriasis; negative = towards Normal."),
         x = "SHAP value (impact on model output)", y = "Gene",
         caption = paste0("Orange = high expression, blue = low. All ", n_final_genes,
                          " genes received non-zero attribution (inherent to KernelSHAP).\nApproximate method, non-converged; probability scale.")) +
    dissertation_theme() +
    theme(plot.caption = element_text(hjust = 0.5, size = 12, colour = "grey50", family = "Arial"),
          plot.caption.position = "plot"),
  low = "#0072B2", mid = "white", high = "#E69F00", midpoint = 0.5, name = "Gene\nexpression\nlevel")
ggsave(out("svm_shap_beeswarm_top10.png"), p_bee_svm, width = 10, height = 7, dpi = 300, bg = "white")

# (2) Bar chart (ranking indicative only, given non-convergence)
svm_shap_bar_df <- data.frame(Gene = factor(names(svm_shap_ranked)[1:10], levels = rev(names(svm_shap_ranked)[1:10])),
                              Mean_Abs_SHAP = as.numeric(svm_shap_ranked[1:10]))
p_bar_svm <- ggplot(svm_shap_bar_df, aes(x = Mean_Abs_SHAP, y = Gene)) +
  geom_col(fill = "#0072B2") +
  geom_text(aes(label = sprintf("%.3f", Mean_Abs_SHAP), x = Mean_Abs_SHAP + max(Mean_Abs_SHAP) * 0.02),
            fontface = "bold", hjust = 0, size = 3.5, colour = "grey30", family = "Arial") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Support Vector Machine SHAP Mean Absolute Importance \u2014 Top 10 Genes",
       subtitle = paste0("Test set (n = ", nrow(svm_test_matrix), "). Mean absolute SHAP across all test samples."),
       x = "Mean |SHAP value|", y = "Gene",
       caption = paste0("All ", n_final_genes, " genes received non-zero attribution (inherent to KernelSHAP).\nApproximate method, non-converged; probability scale.")) +
  dissertation_theme() +
  theme(panel.grid.major.y = element_blank(),
        plot.caption = element_text(hjust = 0.5, size = 12, colour = "grey50", family = "Arial"),
        plot.caption.position = "plot")
ggsave(out("svm_shap_bar_top10.png"), p_bar_svm, width = 9, height = 6, dpi = 300, bg = "white")

################################################
# STEP 29: CROSS-MODEL GENE CONVERGENCE TABLE (RQ3)
# Integrates four complementary lines of evidence per gene: fold
# stability, differential expression, SHAP rank per model, and top-10
# membership per model.
################################################

recur_tab <- table(unlist(per_fold_genes))

shap_rank <- function(shap_vec, genes) {
  nz <- shap_vec[shap_vec > 0]
  if (length(nz) == 0) return(rep(NA_integer_, length(genes)))
  rk  <- rank(-nz, ties.method = "min")
  idx <- match(genes, names(rk))
  as.integer(ifelse(is.na(idx), NA_integer_, rk[idx]))
}
top_n_nonzero <- function(shap_vec, n = 10) names(utils::head(sort(shap_vec[shap_vec > 0], decreasing = TRUE), n))

rf_top10_genes  <- top_n_nonzero(rf_mean_abs_shap)
xgb_top10_genes <- top_n_nonzero(mean_abs_shap)
svm_top10_genes <- top_n_nonzero(svm_mean_abs_shap)

all_genes <- Reduce(union, list(
  names(rf_mean_abs_shap)[rf_mean_abs_shap   > 0],
  names(mean_abs_shap)[mean_abs_shap         > 0],
  names(svm_mean_abs_shap)[svm_mean_abs_shap > 0]))

convergence <- data.frame(
  Gene = all_genes,
  N_folds = as.integer(recur_tab[all_genes]),
  RF_rank  = shap_rank(rf_mean_abs_shap,  all_genes),
  XGB_rank = shap_rank(mean_abs_shap,     all_genes),
  SVM_rank = shap_rank(svm_mean_abs_shap, all_genes),
  RF_top10  = all_genes %in% rf_top10_genes,
  XGB_top10 = all_genes %in% xgb_top10_genes,
  SVM_top10 = all_genes %in% svm_top10_genes,
  stringsAsFactors = FALSE)
convergence$N_folds[is.na(convergence$N_folds)] <- 0
convergence$N_top10 <- rowSums(convergence[, c("RF_top10", "XGB_top10", "SVM_top10")])

de_stats <- limma_results[, c("gene", "logFC", "adj.P.Val")]
convergence <- merge(convergence, de_stats, by.x = "Gene", by.y = "gene", all.x = TRUE)
convergence$logFC     <- round(convergence$logFC, 3)
convergence$adj.P.Val <- signif(convergence$adj.P.Val, 3)
convergence$DE_sig    <- !is.na(convergence$adj.P.Val) & convergence$adj.P.Val < 0.05 & abs(convergence$logFC) > 1
convergence$DE_direction <- ifelse(is.na(convergence$logFC), NA, ifelse(convergence$logFC > 0, "Up", "Down"))

convergence$N_criteria <- as.integer(convergence$DE_sig) +
  as.integer(convergence$N_folds == 5) + as.integer(convergence$N_top10 >= 1)
convergence$RQ3_gate <- convergence$N_criteria == 3

convergence <- convergence[order(-convergence$RQ3_gate, -convergence$N_top10,
                                 -convergence$N_folds, -abs(convergence$logFC)), ]
convergence <- convergence[, c("Gene", "RQ3_gate", "N_criteria", "N_top10", "N_folds",
                               "DE_sig", "DE_direction", "logFC", "adj.P.Val",
                               "RF_rank", "XGB_rank", "SVM_rank", "RF_top10", "XGB_top10", "SVM_top10")]
write.csv(convergence, out("gene_convergence_table.csv"), row.names = FALSE)

rq3_genes <- convergence$Gene[convergence$RQ3_gate]
cat("\nRQ3 candidates (n =", length(rq3_genes), "):", paste(rq3_genes, collapse = ", "), "\n")
write.csv(convergence[convergence$RQ3_gate, ], out("rq3_candidate_biomarkers.csv"), row.names = FALSE)

relaxed <- convergence$Gene[convergence$DE_sig & convergence$N_folds >= 4 & convergence$N_top10 >= 1]
cat("Relaxed gate (>= 4/5 folds): n =", length(relaxed), "|",
    paste(setdiff(relaxed, rq3_genes), collapse = ", "), "additional\n")

cat("Top-10 in all three models:",
    paste(Reduce(intersect, list(rf_top10_genes, xgb_top10_genes, svm_top10_genes)), collapse = ", "), "\n")