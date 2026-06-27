# =============================================================================
# SETUP A: CUSTOM 4-SAMPLE COUNT MATRIX ANALYSIS
# Course: Genomics, Epigenetics and Applications — EFREI Paris
# =============================================================================
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(EnhancedVolcano)
library(AnnotationDbi)
library(org.Hs.eg.db)

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables",  recursive = TRUE, showWarnings = FALSE)


# ── 1. CLEAN LOADING, NUMERIC CONVERSION & NA FILTERING ──────────────────────
# Clear old cached variables from active memory cache
rm(list = ls(all.names = TRUE)) 

# Load custom matrix file directly as a data frame
counts_raw <- read.table(
  "data/counts_matrix.tsv", 
  header      = TRUE, 
  sep         = "\t", 
  check.names = FALSE
)

# 1. Force the first column (the Entrez IDs) to become the row names
rownames(counts_raw) <- counts_raw[, 1]

# 2. Slice out columns 2 to 5 (the actual numeric expression data)
counts_clean <- counts_raw[, 2:5]

# 3. Convert data frame to matrix and force numeric mode
counts <- as.matrix(counts_clean)
mode(counts) <- "numeric"

# 4. CRITICAL FIX: Strip away any rows containing NA values
counts <- counts[rowSums(is.na(counts)) == 0, ]

# 5. Re-assign clean column headers
colnames(counts) <- c("ctrl_1", "dex_1", "ctrl_2", "dex_2")

# VERIFICATION CHECK: This must print "numeric" and show 0 NAs
cat("Count matrix mode:", mode(counts), "\n")
cat("Number of NA cells remaining:", sum(is.na(counts)), "\n")


# ── 2. METADATA SETUP ────────────────────────────────────────────────────────
col_data <- data.frame(
  condition = factor(c("untreated", "dex", "untreated", "dex"), levels = c("untreated", "dex")),
  row.names = colnames(counts)
)
stopifnot(all(colnames(counts) == rownames(col_data)))


# ── 3. FILTERING & ANALYSIS ──────────────────────────────────────────────────
counts_filtered <- counts[rowSums(counts >= 10) >= 2, ]

# Setup the DESeq2 Object
dds <- DESeqDataSetFromMatrix(countData = counts_filtered, colData = col_data, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "dex", "untreated"), alpha = 0.05)


# ── 4. ANNOTATION WITH EXPLORATORY FALLBACK ──────────────────────────────────
res_df <- as.data.frame(res) %>%
  rownames_to_column("entrez_id") %>%
  mutate(
    symbol = mapIds(org.Hs.eg.db, keys = entrez_id, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first"),
    label  = ifelse(!is.na(symbol), symbol, entrez_id),
    direction = case_when(
      sum(padj < 0.05, na.rm = TRUE) > 0 & padj < 0.05 & log2FoldChange >= log2(1.5) ~ "Up in Dex",
      sum(padj < 0.05, na.rm = TRUE) > 0 & padj < 0.05 & log2FoldChange <= -log2(1.5) ~ "Down in Dex",
      sum(padj < 0.05, na.rm = TRUE) == 0 & pvalue < 0.05 & log2FoldChange >= log2(1.5) ~ "Up (p < 0.05 exploratory)",
      sum(padj < 0.05, na.rm = TRUE) == 0 & pvalue < 0.05 & log2FoldChange <= -log2(1.5) ~ "Down (p < 0.05 exploratory)",
      TRUE ~ "Not significant"
    )
  ) %>% arrange(pvalue)

write.csv(res_df, "results/tables/deseq2_results_custom_4samples.csv", row.names = FALSE)


# ── 5. FIGURES GENERATION ────────────────────────────────────────────────────
vst_data <- vst(dds, blind = FALSE)

# Figure 1: PCA Plot
pca_data <- plotPCA(vst_data, intgroup = "condition", returnData = TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"), 1)
p1 <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) + geom_point(size = 5) +
  labs(title = "PCA Plot (Custom 4 Samples)", x = paste0("PC1 (", pct_var[1], "%)"), y = paste0("PC2 (", pct_var[2], "%)")) + theme_bw()
ggsave("results/figures/Custom_Figure1_PCA.png", p1, width = 7, height = 5)

# Figure 2: Volcano Plot
p2 <- EnhancedVolcano(res_df, lab = res_df$label, x = "log2FoldChange", y = "pvalue",
                      title = "Volcano Plot (Custom 4 Samples)", pCutoff = 0.05, FCcutoff = log2(1.5), legendPosition = "right")
ggsave("results/figures/Custom_Figure2_Volcano.png", p2, width = 9, height = 7)

# Figure 3: Heatmap (Zero-variance protected)
top50_genes <- res_df %>% filter(!is.na(pvalue)) %>% arrange(pvalue) %>% slice_head(n = 50) %>% pull(entrez_id)
vst_mat     <- assay(vst_data)[top50_genes, ]
vst_mat     <- vst_mat[apply(vst_mat, 1, var) > 0, ]
rownames(vst_mat) <- res_df$label[match(rownames(vst_mat), res_df$entrez_id)]
vst_scaled  <- t(scale(t(vst_mat)))
vst_scaled[vst_scaled > 3] <- 3; vst_scaled[vst_scaled < -3] <- -3

png("results/figures/Custom_Figure3_Heatmap.png", width = 2000, height = 2500, res = 300)
pheatmap(vst_scaled, color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100), clustering_method = "ward.D2")
dev.off()

cat("✓ Custom analysis complete! Check your results/ folder.\n")