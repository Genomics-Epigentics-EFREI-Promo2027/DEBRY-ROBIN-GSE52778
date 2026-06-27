# =============================================================================
# SETUP B: AIRWAY COMPLETE 9-FIGURE GRAPHICS PIPELINE
# =============================================================================
suppressPackageStartupMessages({
  library(DESeq2)
  library(airway)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
  library(tibble)
  library(scales)
})

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables",  recursive = TRUE, showWarnings = FALSE)

# 1. PIPELINE INITIALIZATION & DESIGN FORMULA
data(airway)
dds <- DESeqDataSet(airway, design = ~ cell + dex)
dds$dex <- relevel(dds$dex, ref = "untrt")
dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]

# 2. RUN MODEL
dds      <- DESeq(dds)
res      <- results(dds, contrast = c("dex", "trt", "untrt"), alpha = 0.05, pAdjustMethod = "BH")
vst_data <- vst(dds, blind = FALSE)

# 3. ANNOTATE DATA
res_df <- as.data.frame(res) %>% tibble::rownames_to_column("ensembl_id") %>% dplyr::filter(!is.na(padj))
gene_anno <- suppressMessages(
  AnnotationDbi::select(x = org.Hs.eg.db, keys = res_df$ensembl_id, columns = c("ENSEMBL", "SYMBOL", "ENTREZID"), keytype = "ENSEMBL")
) %>% dplyr::distinct(ENSEMBL, .keep_all = TRUE) %>% dplyr::rename(ensembl_id = ENSEMBL)

res_ann <- res_df %>% dplyr::left_join(gene_anno, by = "ensembl_id") %>%
  dplyr::mutate(
    label = dplyr::if_else(!is.na(SYMBOL), SYMBOL, ensembl_id),
    direction = dplyr::case_when(padj < 0.05 & log2FoldChange > 1 ~ "UP", padj < 0.05 & log2FoldChange < -1 ~ "DOWN", TRUE ~ "NS"),
    neg_log10_padj = -log10(padj)
  )
write.csv(res_ann, "results/tables/DEG_annotated_results.csv", row.names = FALSE)

# ── FIGURES GENERATION ──

# Figure 1: Volcano Plot
top_labels <- res_ann %>% dplyr::filter(direction != "NS") %>% dplyr::arrange(padj) %>% dplyr::slice_head(n = 15) %>% dplyr::pull(ensembl_id)
res_ann <- res_ann %>% dplyr::mutate(plot_label = dplyr::if_else(ensembl_id %in% top_labels, label, ""))
fig1 <- ggplot(res_ann %>% dplyr::arrange(direction), aes(x = log2FoldChange, y = neg_log10_padj, color = direction, label = plot_label)) +
  geom_point(aes(size = direction, alpha = direction)) + geom_hline(yintercept = -log10(0.05), linetype = "dashed") + geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  ggrepel::geom_text_repel(data = dplyr::filter(res_ann, plot_label != ""), color = "black", fontface = "bold", max.overlaps = 30) +
  scale_color_manual(values = c("UP"="#E41A1C", "DOWN"="#377EB8", "NS"="grey70")) + scale_size_manual(values=c(UP=1.8, DOWN=1.8, NS=0.8)) + scale_alpha_manual(values=c(UP=0.9, DOWN=0.9, NS=0.4)) + theme_bw()
ggsave("results/figures/Figure1_Volcano.png", fig1, width = 10, height = 8)

# Figure 2: MA Plot
fig2 <- ggplot(res_ann, aes(x = log2(baseMean + 1), y = log2FoldChange, color = direction)) + geom_point(alpha = 0.4) + geom_hline(yintercept = 0, color = "red") + theme_bw()
ggsave("results/figures/Figure2_MA_Plot.png", fig2, width = 10, height = 7)

# Figure 3: Heatmap Top 50 DEGs
top50_ids <- res_ann %>% dplyr::filter(direction != "NS") %>% dplyr::arrange(padj) %>% dplyr::slice_head(n = 50) %>% dplyr::pull(ensembl_id)
vst_mat   <- assay(vst_data)[top50_ids, ]
rownames(vst_mat) <- res_ann$label[match(rownames(vst_mat), res_ann$ensembl_id)]
vst_scaled <- t(scale(t(vst_mat)))
col_anno <- data.frame(Treatment = dplyr::if_else(dds$dex == "trt", "Treated", "Untreated"), Cell_Line = as.character(dds$cell), row.names = colnames(vst_scaled))
png("results/figures/Figure3_Heatmap_Top50DEG.png", width = 4000, height = 5000, res = 300)
pheatmap(vst_scaled, annotation_col = col_anno, clustering_method = "ward.D2")
dev.off()

# Figure 4: PCA Plot
pca_data <- plotPCA(vst_data, intgroup = c("dex", "cell"), returnData = TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"), 1)
fig4 <- ggplot(pca_data, aes(x = PC1, y = PC2, color = dex, shape = cell)) + geom_point(size = 5) + labs(x = paste0("PC1 (", pct_var[1], "%)"), y = paste0("PC2 (", pct_var[2], "%)")) + theme_bw()
ggsave("results/figures/Figure4_PCA.png", fig4, width = 9, height = 7)

# Figure 5: Sample Distance Matrix
samp_mat <- as.matrix(dist(t(assay(vst_data))))
png("results/figures/Figure5_Sample_Distance.png", width = 2800, height = 2400, res = 300)
pheatmap(samp_mat, color = colorRampPalette(rev(brewer.pal(9, "Blues")))(100), clustering_method = "ward.D2")
dev.off()

# Figure 6 & 7: GO and KEGG Functional Enrichment Pathways
up_genes <- res_ann %>% dplyr::filter(direction == "UP", !is.na(ENTREZID)) %>% dplyr::pull(ENTREZID)
universe <- res_ann %>% dplyr::filter(!is.na(ENTREZID)) %>% dplyr::pull(ENTREZID)

go_up <- clusterProfiler::enrichGO(gene = up_genes, universe = universe, OrgDb = org.Hs.eg.db, ont = "BP", pvalueCutoff = 0.05, readable = TRUE)
if(!is.null(go_up) && nrow(go_up@result) > 0) {
  ggsave("results/figures/Figure6a_GO_UP.png", enrichplot::dotplot(go_up, showCategory = 20), width = 11, height = 9)
}
kegg_up <- clusterProfiler::enrichKEGG(gene = up_genes, universe = universe, organism = "hsa", pvalueCutoff = 0.05)
if(!is.null(kegg_up) && nrow(kegg_up@result) > 0) {
  ggsave("results/figures/Figure7a_KEGG_UP.png", enrichplot::dotplot(kegg_up, showCategory = 20), width = 11, height = 9)
}

# Figure 8: Top 12 Expression Boxplots
top12_ids <- res_ann %>% dplyr::filter(direction != "NS") %>% dplyr::arrange(padj) %>% dplyr::slice_head(n = 12) %>% dplyr::pull(ensembl_id)
box_df <- as.data.frame(t(counts(dds, normalized=TRUE)[top12_ids,])) %>% rownames_to_column("sample") %>%
  tidyr::pivot_longer(cols = -sample, names_to = "ensembl_id", values_to = "count") %>%
  dplyr::mutate(Treatment = dplyr::if_else(sample %in% colnames(dds)[dds$dex == "trt"], "Treated", "Untreated"), log_count = log2(count + 1))
fig8 <- ggplot(box_df, aes(x = Treatment, y = log_count, fill = Treatment)) + geom_boxplot() + facet_wrap(~ensembl_id, scales = "free_y") + theme_bw()
ggsave("results/figures/Figure8_Top12_Boxplots.png", fig8, width = 14, height = 10)

# Figure 9: Summary Statistics Bar Chart
summary_df <- data.frame(Category = c("Total tested", "Significant", "Upregulated", "Downregulated"), Count = c(nrow(res_ann), sum(res_ann$padj < 0.05), sum(res_ann$direction == "UP"), sum(res_ann$direction == "DOWN")))
fig9 <- ggplot(summary_df, aes(x = Category, y = Count, fill = Category)) + geom_col() + geom_text(aes(label = Count), vjust = -0.5) + theme_bw()
ggsave("results/figures/Figure9_Summary_Bar.png", fig9, width = 9, height = 7)

cat("\n ╔══════════════════════════════════════════════════════╗ \n ║ PIPELINE RUN SUCCESSFULLY WITHOUT ERRORS! ║ \n ╚══════════════════════════════════════════════════════╝ \n")