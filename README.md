# DEBRY-ROBIN-GSE52778

## RNA-seq Analysis of Dexamethasone Response in Human Airway Smooth Muscle Cells

**Reproduction of Himes et al. (2014) - GSE52778**

> This project reproduces the differential gene expression analysis published by
> **Himes et al. (2014)** using the RNA-seq dataset **GSE52778**. The objective is to
> identify genes and biological pathways regulated by dexamethasone in human
> airway smooth muscle (**HASM**) cells, with particular attention to the novel
> glucocorticoid-responsive gene **CRISPLD2** and its implications for asthma
> pharmacogenetics.

---

## Course Information

| Field | Details |
|---|---|
| **Course** | ST2GEA – Genomics, Epigenetics and Applications (ST2GEA-2526PSP01) |
| **Supervisor** | MATHEW Mano Joseph |
| **Academic year** | 2025–2026 |
| **Author** | Robin DEBRY - 20240681 - ING2-BIO |

---

## Biological Question

**Which genes are differentially expressed in human airway smooth muscle cells
treated with dexamethasone, and what biological pathways do they implicate?**

Glucocorticoids such as dexamethasone are among the most widely prescribed
anti-inflammatory drugs for asthma and COPD. The landmark finding of Himes et
al. was the identification of **CRISPLD2** as a novel glucocorticoid-responsive
gene acting as a negative modulator of the inflammatory response - a result with
direct pharmacogenetic implications. This project tests whether that finding is
reproducible with a modern bioinformatics pipeline (**Galaxy Europe + DESeq2**),
using both the full 8-sample dataset and a 4-sample reproducibility subset.

---

## Dataset

| Property | Value |
|---|---|
| **GEO accession** | [GSE52778](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778) |
| **Biological system** | Human airway smooth muscle (HASM) cells |
| **Conditions** | Dexamethasone (1 µM, 18h) vs Untreated |
| **Donors** | 4 independent HASM cell lines (N61311, N052611, N080611, N061011) |
| **Design** | Paired - 4 treated / 4 untreated (8 samples total) |
| **Genes tested** | ~14,224 (after low-count filtering) |
| **Sequencing platform** | Illumina HiSeq 2000 (paired-end) |

### Samples

| SRA Accession | Cell Line | Condition | Analysis |
|---|---|---|---|
| SRR1039508 | N61311 | Untreated | Full |
| SRR1039509 | N61311 | Dexamethasone | Full |
| SRR1039512 | N052611 | Untreated | Full + Custom |
| SRR1039513 | N052611 | Dexamethasone | Full + Custom |
| SRR1039516 | N080611 | Untreated | Full |
| SRR1039517 | N080611 | Dexamethasone | Full |
| SRR1039520 | N061011 | Untreated | Full |
| SRR1039521 | N061011 | Dexamethasone | Full |

---

## Repository Structure

```text
.
├── RNA_seq_GSE52778_Report_FV.Rmd
├── RNA_seq_GSE52778_Report_FV.html
├── rnaseq_project.Rproj
├── README.md
│
├── data/
├── docs/
├── galaxy_workflow/
├── r_analysis/
│
└── results/
    └── figures/
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
```

---

## Analysis Pipeline

Two complementary analyses were performed:

- **Full analysis (8 samples, 4 donors):** complete dataset with donor as a
  blocking factor in the DESeq2 model (design = ~ cell + condition)
- **Custom 4-sample analysis (2 donors):** independent Galaxy Europe run on
  SRR1039508/09 and SRR1039512/13, demonstrating end-to-end pipeline
  reproducibility without a blocking factor

### Phase A - Galaxy Europe

| Step | Tool | Version | Purpose |
|---|---|---|---|
| 1 | SRA Download | - | Fetch raw FASTQ reads |
| 2 | FastQC | 0.74 | Per-sample quality metrics |
| 3 | Trim Galore | 0.6.7 | Adapter removal, quality trimming |
| 4 | HISAT2 | 2.2.1 | Splice-aware alignment to hg19 |
| 5 | featureCounts | 2.0.3 | Gene-level count matrix |
| 6 | MultiQC | 1.11 | Aggregate QC report |

Alignment target: **hg19 / GRCh37**, Ensembl release 75 annotation.

### Phase B - R / DESeq2

1. Data import and low-count filtering (rowSums ≥ 10)
2. DESeq2 model fitting (design = ~ cell + condition)
3. LFC shrinkage (apeglm estimator)
4. Quality control
   - Sample-to-sample distance heatmap (VST-normalised)
   - Principal Component Analysis (PCA)
5. Differential expression visualisation
   - Volcano plot
   - MA plot
   - Heatmap of top 50 DEGs
   - Gene expression boxplots (top 12 DEGs)
6. Functional enrichment analysis
   - Gene Ontology - Biological Process (clusterProfiler)
   - KEGG pathway enrichment (clusterProfiler)
7. Deep dive: CRISPLD2 as landmark reproducibility target

---

## Key Parameters

| Parameter | Value |
|---|---|
| **FDR threshold** | 0.05 |
| **Log2FC threshold** | 1.0 (minimum two-fold change) |
| **Multiple testing correction** | Benjamini–Hochberg (BH) |
| **Normalisation** | Variance Stabilising Transformation (VST) |
| **LFC shrinkage** | apeglm estimator |

---

## Main Results

- Dexamethasone induced a strong and reproducible transcriptional response,
  with **PC1 (47.5%)** cleanly separating treated from untreated samples
- The upregulation asymmetry (more genes up than down) is consistent with
  glucocorticoid receptor biology: direct GRE-driven activation dominates over
  indirect NF-κB/AP-1 repression
- **Well-established GC targets** recovered as top hits: FKBP5, TSC22D3, DUSP1,
  KLF15 - validating the pipeline
- **CRISPLD2** identified among the most significantly upregulated genes,
  reproducing the landmark finding of Himes et al.
- Top GO term: *cellular response to hormone stimulus* (highest GeneRatio and
  most significant adjusted p-value)
- Top KEGG pathway: *Cytoskeleton in muscle cells* (largest gene count and
  lowest adjusted p-value among enriched pathways)

---

## R Packages Required

```r
# Bioconductor
BiocManager::install(c(
  "DESeq2",
  "org.Hs.eg.db",
  "AnnotationDbi",
  "clusterProfiler",
  "enrichplot"
))

# CRAN
install.packages(c(
  "ggplot2",
  "EnhancedVolcano",
  "pheatmap",
  "RColorBrewer",
  "dplyr",
  "tidyr",
  "tibble"
))
```

---

## How to Run

1. **Clone the repository**

```bash
git clone https://github.com/Genomics-Epigenetics-EFREI-Promo2027/DEBRY-ROBIN-GSE52778.git
cd DEBRY-ROBIN-GSE52778
```

2. **Download the raw count matrix** from
   [GEO - GSE52778](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778)
   and place it in **data/**.

3. **Open the RStudio project** by double-clicking **rnaseq_project.Rproj** -
   this sets the working directory automatically.

4. **Open RNA_seq_GSE52778_Report_FV.Rmd** and click **Knit**, or run:

```r
rmarkdown::render("RNA_seq_GSE52778_Report_FV.Rmd")
```

All figures will be generated in **results/figures/**.

---

## Note on Pipeline Divergence

The original Himes et al. study used **TopHat + Cufflinks/Cuffdiff** with RefSeq
annotation. This re-analysis uses the more modern **HISAT2 / featureCounts /
DESeq2** pipeline with Ensembl release 75 annotation, offering improved accuracy
and more robust statistical modelling. Minor differences in the total DEG count
relative to the original 316 are expected and reflect methodological improvements
rather than inconsistencies in the underlying data.

---

## Reference

Himes BE, Jiang X, Wagner P, Hu R, Wang Q, Klanderman B, Whitaker RM, Duan Q,
Lasky-Su J, Nikolos C, Jester W, Johnson M, Panettieri RA Jr, Tantisira KG,
Weiss ST, Lu Q. (2014). RNA-Seq Transcriptome Profiling Identifies CRISPLD2 as
a Glucocorticoid Responsive Gene that Modulates Cytokine Function in Airway
Smooth Muscle Cells. *PLoS ONE*, 9(6):e99625.
**doi**: [10.1371/journal.pone.0099625](https://doi.org/10.1371/journal.pone.0099625)
**PMID**: 24926665.

---

## Author

**Robin DEBRY**  
*Student ID*: 20240681  
ING2-BIO - EFREI Paris  
Academic Year 2025–2026  
Final Project - ST2GEA-2526PSP01  
*Supervisor*: **MATHEW Mano Joseph**
