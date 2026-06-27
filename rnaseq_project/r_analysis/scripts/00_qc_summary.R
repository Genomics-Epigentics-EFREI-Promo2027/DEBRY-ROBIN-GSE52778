# =============================================================================
# DEPENDENCY INSTALLATION SCRIPT
# Course: Genomics, Epigenetics and Applications — EFREI Paris
# INSTRUCTION: Run this script once to ensure all Bioconductor and CRAN
#              packages are perfectly configured before launching pipelines.
# =============================================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

required_pkgs <- c(
  "DESeq2", "airway", "org.Hs.eg.db", "AnnotationDbi",
  "clusterProfiler", "enrichplot", "ggplot2", "ggrepel", 
  "pheatmap", "patchwork", "RColorBrewer", "viridis", 
  "dplyr", "tibble", "scales", "tidyr"
)

cat("Checking and installing missing dependencies...\n")

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing package:", pkg, "\n")
    BiocManager::install(pkg, ask = FALSE)
  } else {
    cat("✓ Package already present:", pkg, "\n")
  }
}

cat("\n═══════════════════════════════════════════\n")
cat(" ALL DEPENDENCIES SUCCESSFULLY CONFIGURED!\n")
cat("═══════════════════════════════════════════\n")