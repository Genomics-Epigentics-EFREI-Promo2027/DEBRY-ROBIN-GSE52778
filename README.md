# DEBRY-ROBIN-GSE52778

Final Project of Genomics and Epigenetics — RNA-seq analysis of Dexamethasone Response in Airway Smooth Muscle Cells

This project reproduces the differential gene expression analysis published by Himes et al. (2014) using the RNA-seq dataset GSE52778. The objective is to identify genes and biological pathways regulated by dexamethasone in human airway smooth muscle (HASM) cells, with particular attention to the novel glucocorticoid-responsive gene CRISPLD2 and its implications for asthma pharmacogenetics.

---

## Project structure

.

├── RNA_seq_GSE52778_Report_FV.Rmd   # Main analysis notebook

├── RNA_seq_GSE52778_Report_FV.html  # Rendered report

├── README.md

├── data/                            # Input count matrix (from Galaxy/featureCounts)

├── docs/                            # Supporting documentation

├── galaxy_workflow/                 # Galaxy Europe workflow and MultiQC reports

├── r_analysis/                      # R scripts and intermediate objects

└── results/

└── figures/                     # All generated plots (PNG)

├── Figure1_Volcano.png

├── Figure2_MA_Plot.png

├── Figure3_Heatmap_Top50DEG.png

├── Figure4_PCA.png

├── Figure5_Sample_Distance.png

├── Figure6a_GO_UP.png

├── Figure7a_KEGG_UP.png

├── Figure8_Top12_Boxplots.png

├── Figure9_Summary_Bar.png

├── Custom_Figure1_PCA.png

└── Custom_Figure2_Volcano.png

---

## Dataset

- **GEO accession:** GSE52778
- **Biological system:** Human airway smooth muscle (HASM) cells
- **Conditions:** Dexamethasone-treated (1 µM, 18h) vs Untreated
- **Donors:** 4 independent HASM cell lines (N61311, N052611, N080611, N061011)
- **Samples:** 8 (paired design — 4 treated, 4 untreated)
- **Genes tested:** ~14,224 (after low-count filtering)
- **Sequencing platform:** Illumina HiSeq 2000 (paired-end)

---

## Analysis pipeline

Two complementary analyses were performed:

- **Full analysis (8 samples, 4 donors):** complete dataset with donor as a blocking factor
- **Custom 4-sample analysis (2 donors):** end-to-end Galaxy Europe pipeline reproducibility check

### Phase A — Galaxy Europe

1. Raw FASTQ download from NCBI SRA
2. Quality control (FastQC)
3. Adapter trimming (Trim Galore)
4. Splice-aware alignment to hg19/GRCh37 (HISAT2)
5. Gene-level count matrix generation (featureCounts)
6. Aggregate QC report (MultiQC)

### Phase B — R / DESeq2

1. Data import and preprocessing
2. Low-count filtering (rowSums ≥ 10)
3. DESeq2 differential expression analysis (`design = ~ cell + condition`)
4. Quality control and visualisation
   - Sample-to-sample distance heatmap
   - Principal Component Analysis (PCA)
5. Differential expression visualisation
   - Volcano plot
   - MA plot
   - Heatmap of top 50 DEGs
   - Gene expression boxplots (top 12 DEGs)
6. Functional enrichment analysis
   - Gene Ontology (Biological Process)
   - KEGG pathway enrichment
7. Deep dive into the landmark finding: CRISPLD2

---

## Key parameters

- **FDR threshold:** 0.05
- **|log2FC| threshold:** 1 (minimum two-fold change)
- **Multiple testing correction:** Benjamini–Hochberg (BH)
- **Normalisation:** Variance Stabilising Transformation (VST)
- **LFC shrinkage:** apeglm estimator

---

## R packages required

- `DESeq2`
- `org.Hs.eg.db`
- `AnnotationDbi`
- `clusterProfiler`
- `EnhancedVolcano`
- `ggplot2`
- `pheatmap`
- `RColorBrewer`
- `dplyr`
- `tidyr`
- `tibble`

---

## How to run

1. Open `rnaseq_project.Rproj` in RStudio — this sets the working directory automatically.
2. Ensure the count matrix is present in `data/`.
3. Open `RNA_seq_GSE52778_Report_FV.Rmd` and click **Knit**.

All figures will be automatically generated in:
- `results/figures/`

---

## Reference

Himes BE, Jiang X, Wagner P, Hu R, Wang Q, Klanderman B, Whitaker RM, Duan Q, Lasky-Su J, Nikolos C, Jester W, Johnson M, Panettieri RA Jr, Tantisira KG, Weiss ST, Lu Q. (2014). RNA-Seq Transcriptome Profiling Identifies CRISPLD2 as a Glucocorticoid Responsive Gene that Modulates Cytokine Function in Airway Smooth Muscle Cells. *PLoS ONE*, 9(6):e99625. doi: 10.1371/journal.pone.0099625. PMID: 24926665.

---

## Author

Robin DEBRY — ING2-BIO  
Final Project — ST2GEA-2526PSP01  
Supervisor: MATHEW Mano Joseph
