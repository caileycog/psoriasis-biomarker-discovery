# File Guide

This page explains every file in the repository: what it is, where it comes from, and which step of `R/01_analysis.R` produces it. For the project overview, see the [README](../README.md).

---

## 1. How to run the code

1. **Install R (≥ 4.3) and RStudio.**
2. **Download the repository.** Either use the green **Code → Download ZIP** button on GitHub, or run:
   ```bash
   git clone https://github.com/<caileycog>/psoriasis-biomarker-discovery.git
   ```
3. **Open `psoriasis-biomarker-discovery.Rproj`** in RStudio. This sets the working directory to the repository root, so every file path in the code resolves correctly.
4. **Install the packages (once):**
   ```r
   source("R/00_install_packages.R")
   ```
5. **Run the analysis:**
   ```r
   source("R/01_analysis.R")
   ```
   Or open the file and run it section by section with `Ctrl/Cmd + Enter` to follow each step's printed output.
6. **Find the results** in `figures/` (PNG) and `tables/` (CSV).

**Runtime:** roughly 30–60 minutes on a standard laptop. Step 28 (SVM Kernel SHAP) takes the longest.

**Troubleshooting**

| Problem | Fix |
|---|---|
| `cannot open file 'data/raw/RPKM_Values.txt'` | Your working directory isn't the repo root. Open the `.Rproj` file, or run `setwd("path/to/psoriasis-biomarker-discovery")`. |
| Error in `xgbTree` / caret | Check `packageVersion("xgboost")` returns `1.7.8.1`, then re-run `00_install_packages.R`. |
| Font warnings about "Arial" | Harmless. The plots fall back to the default font. |

---

## 2. Code (`R/`)

| File | Purpose |
|---|---|
| `00_install_packages.R` | Installs all CRAN, Bioconductor (`limma`) and GitHub (`nestedcv`) dependencies, then pins `xgboost` to 1.7.8.1. |
| `01_analysis.R` | The full pipeline, organised into numbered steps (summarised below). |

### Pipeline steps in `01_analysis.R`

| Step | What it does | Dissertation section |
|---|---|---|
| 2–3 | Load libraries, set seed | 3.7 |
| 4–7 | Load data, match sample IDs, merge, check NAs, zeros and duplicates, drop all-zero genes | 3.2 |
| 8–11 | Descriptive statistics, raw density plot, log2(x+1) transform | 3.3.1–3.3.2 |
| 12A–12C | Per-sample boxplots, sample-correlation outlier check, PCA | 3.3.3–3.3.5 |
| 13–14 | limma differential expression + volcano plot (**RQ1**) | 3.4, 4.1 |
| 15 | Stratified 80/20 train/test split | 3.5.1 |
| 16 | Low-variance gene filter (fit on training data only) | 3.5.2 |
| 17–18 | Build model matrices; create shared outer CV folds | 3.5.4 |
| 19A–19C | Elastic-net α diagnostic, filter function, per-fold gene selection | 3.5.4.1 |
| 20–22 | Nested CV for XGBoost, Random Forest and SVM | 3.5.3–3.5.4.2 |
| 23–25 | CV metrics, test-set metrics, CV vs test overfitting check (**RQ2**) | 3.5.5, 4.2 |
| 26–28 | SHAP for XGBoost (Tree), Random Forest (Tree), SVM (Kernel) | 3.5.6, 4.2 |
| 29 | Cross-model convergence table → candidate biomarkers (**RQ3**) | 3.6, 4.3 |

---

## 3. Input data (`data/`)

| File | Description |
|---|---|
| `raw/RPKM_Values.txt` | Gene expression matrix from GEO **GSE54456**. Tab-separated; rows = 21,510 genes, columns = 174 samples; values are RPKM (Reads Per Kilobase per Million mapped reads). First column = gene symbol. |
| `raw/sample_metadata.csv` | One row per sample: `Sample ID` (matches the RPKM column names), `GEO Accession Number`, and `Tissue Type` (`Normal` or `Psoriasis`). |
| `disgenet/<GENE>_disgenet.xlsx` | DISGENET gene–disease association exports for each of the 14 candidate biomarkers, used for Table 4.6. Key columns: `Disease`, `score` (GDA score, 0–1), `N PMIDs` (supporting publications), `DSI g` (Disease Specificity Index). `IL1F9_IL36G` and `VNN3_VNN3P` are filed under their current HGNC symbols in DISGENET. |

---

## 4. Outputs (`figures/` and `tables/`)

The analysis script regenerates everything in these two folders. The full differential expression results (`limma_results_all_genes.csv`) and the two significant-gene lists (`limma_sig_genes_upregulated.csv` and `limma_sig_genes_downregulated.csv`) are included so the results can be browsed without running the code. Run `01_analysis.R` to produce the rest.

**Columns in the limma tables:** `log2FC` (log2 fold change, psoriasis vs normal), `AveExpr` (mean log2 expression), `t` (moderated t-statistic), `P.Value` (raw p-value), `adj.P.Val` (Benjamini–Hochberg FDR-adjusted p-value), `B` (log-odds that the gene is differentially expressed), `gene` (gene symbol).

### Figures (`figures/`)

| File | Step | What it shows |
|---|---|---|
| `density_rpkm_untransformed.png` | 9 | Extreme right skew of raw RPKM, which motivates the log transform |
| `boxplot_log2_transformed.png` | 12A | Per-sample distributions after log2 (no outlier samples) |
| `density_log2_transformed.png` | 12A | Distribution after log2(x+1) |
| `pca_plot.png` | 12C | PCA: psoriatic and normal samples separate along PC1 |
| `volcano_plot.png` | 14 | Differential expression: 1,131 significant genes |
| `variance_histogram_training.png` | 16 | Justification for the variance-filter threshold |
| `cv_roc_overlay.png` | 23 | Nested-CV ROC curves, all three models |
| `test_set_roc_overlay.png` | 24 | Held-out test ROC curves, all three models |
| `xgb_shap_bar_top10.png` / `xgb_shap_beeswarm_top10.png` | 26 | XGBoost SHAP importance |
| `rf_shap_bar_top10.png` / `rf_shap_beeswarm_top10.png` | 27 | Random Forest SHAP importance |
| `svm_shap_bar_top10.png` / `svm_shap_beeswarm_top10.png` | 28 | SVM SHAP importance |
| `nested_cv_diagram.png` | n/a | Diagram of the nested cross-validation procedure (dissertation Figure 3.6; not generated by the code) |

### Tables (`tables/`)

| File | Step | Contents |
|---|---|---|
| `sample_correlation_summary.csv` | 12B | Mean correlation of each sample with all others (outlier check) |
| `limma_results_all_genes.csv` | 13 | Full limma results for all 21,099 genes: log2FC, p-value, adjusted p-value |
| `limma_sig_genes_upregulated.csv` | 13 | 489 significantly upregulated genes |
| `limma_sig_genes_downregulated.csv` | 13 | 642 significantly downregulated genes |
| `elastic_net_alpha_diagnostic.csv` | 19A | CV error and number of genes selected at each α (0.1–1.0) |
| `genes_per_fold.csv` | 19C | Genes selected by the elastic net in each of the 5 outer folds |
| `final_gene_panel.csv` | 19C | The final 64-gene panel (selected on the full training set) |
| `genes_in_all_folds.csv` | 19C | The 30 genes selected in all 5 folds |
| `xgb_` / `rf_` / `svm_per_fold_auc.csv` | 20–22 | AUC on each outer fold |
| `xgb_` / `rf_` / `svm_final_hyperparameters.csv` | 20–22 | Best hyperparameters of each final model |
| `cv_pooled_metrics_all_models.csv` | 23 | AUC, accuracy, sensitivity, specificity, F1, MCC (nested CV) |
| `cv_confusion_matrices.csv` | 23 | Confusion matrices (nested CV) |
| `test_set_metrics_all_models.csv` | 24 | Same metrics on the 34-sample test set, with bootstrap 95% CI for AUC |
| `test_set_confusion_matrices.csv` | 24 | Confusion matrices (test set) |
| `cv_vs_test_comparison.csv` | 25 | Side-by-side CV vs test metrics (overfitting check) |
| `xgb_shap_top10_genes.csv` | 26 | Top 10 genes by mean \|SHAP\| (XGBoost) |
| `rf_shap_top10_genes.csv`, `rf_shap_all_nonzero_genes.csv` | 27 | Random Forest SHAP rankings |
| `svm_shap_top10_genes.csv`, `svm_shap_all_nonzero_genes.csv` | 28 | SVM SHAP rankings |
| `gene_convergence_table.csv` | 29 | Every panel gene scored against all three evidence lines |
| `rq3_candidate_biomarkers.csv` | 29 | **The 14 final candidate biomarkers** |

---

## 5. Documentation (`docs/`)

| File | Description |
|---|---|
| `dissertation.pdf` | The full MSc dissertation (methods, results, discussion, references) |
| `FILE_GUIDE.md` | This file |
