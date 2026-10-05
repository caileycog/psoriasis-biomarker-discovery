#################################################################################
# 00_install_packages.R
# Run this ONCE before R/01_analysis.R.
#
# Installs every package the analysis needs, then pins xgboost to 1.7.8.1.
# The pin matters: shapviz / treeshap can pull in a newer xgboost, which breaks
# caret's "xgbTree" method. Pinning last guarantees 1.7.8.1 is what remains.
#################################################################################

cran_packages <- c(
  "glmnet", "ranger", "caret", "e1071", "shapviz", "pROC", "ggplot2",
  "ggrepel", "matrixStats", "remotes", "mltools", "kernlab", "treeshap",
  "kernelshap", "scales"
)

missing <- cran_packages[!vapply(cran_packages, requireNamespace,
                                 logical(1), quietly = TRUE)]
if (length(missing) > 0) install.packages(missing)

# nestedcv (GitHub version, used for nested cross-validation)
if (!requireNamespace("nestedcv", quietly = TRUE)) {
  remotes::install_github("myles-lewis/nestedcv")
}

# limma (Bioconductor, used for differential expression)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("limma", quietly = TRUE)) {
  BiocManager::install("limma", update = FALSE, ask = FALSE)
}

# Pin xgboost to the caret-compatible version (must run last)
if (!requireNamespace("xgboost", quietly = TRUE) ||
    packageVersion("xgboost") != "1.7.8.1") {
  remotes::install_version("xgboost", version = "1.7.8.1", upgrade = "never")
}

cat("All packages installed. xgboost version:",
    as.character(packageVersion("xgboost")), "\n")
