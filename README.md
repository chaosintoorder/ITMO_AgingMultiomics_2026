# ITMO_AgingMultiomics_2026

Cross-sectional multi-omics meta-analysis of **non-monotonic human blood aging**.
The pipeline pools public microarray **gene-expression** and Illumina 450K **DNA-methylation**
datasets, recovers age trends per sex, locates the ages of fastest change ("critical ages"),
integrates the two omics layers, and validates the result against an external aging-gene database.

Concept inherited from the **SODA** project (https://github.com/CaptnClementine/SODA),
re-implemented with statistical control (FDR, cell-composition adjustment, batch correction),
a second omics layer, a sex split, and external validation.

---

## Requirements

- **R ≥ 4.3** (developed on R 4.5)
- ~8 GB RAM (the pipeline is memory-aware; methylation steps process one dataset at a time)
- Internet access for the first run (GEO download, KEGG/Reactome/STRING annotation)
- RStudio recommended — open `ITMO_AgingMultiomics_2026.Rproj` so the working directory is set automatically

## Dependencies

Installed in one step by `R/install_packages.R`.

- **CRAN:** data.table, R.utils, ggplot2, ggridges, mgcv, segmented, matrixStats, future.apply,
  e1071, RColorBrewer, UpSetR, readxl, svglite, uwot, igraph
- **Bioconductor:** GEOquery, limma, sva, impute, org.Hs.eg.db,
  hgu133plus2.db, hugene10sttranscriptcluster.db, hugene11sttranscriptcluster.db,
  huex10sttranscriptcluster.db, illuminaHumanv4.db,
  minfi, IlluminaHumanMethylation450kanno.ilmn12.hg19, EpiDISH, missMethyl,
  clusterProfiler, ReactomePA, WGCNA, STRINGdb, Mfuzz

```r
source("R/install_packages.R")   # run once
```

---

## Script description

| File | Role |
|------|------|
| `config.R` | Central configuration: paths, dataset lists (`EXPR_GSE`, `METH_GSE`), all parameters (FDR, age grid, variable-feature counts), colour palettes, platform→annotation maps, methylation file manifest. Sourced by **every** script and calls `setup_session()`. Single source of truth. |
| `R/install_packages.R` | One-time installation of all CRAN and Bioconductor dependencies. |
| `R/utils.R` | Generic helpers: logging, plotting theme, `save_plot()` (png + svg), p-value formatting, age/sex parsing, RDS I/O. |
| `R/utils_io.R` | Heavy I/O and data functions: GEO reading, probe→gene annotation, sex prediction, beta↔M conversion, EpiDISH and marker-based deconvolution, `build_covars()`, ComBat wrapper, PCA/UMAP and distribution plots. |
| `R/utils_analysis.R` | Analysis functions: `age_limma()`, `critical_ages()` (GAM + segmented), `choose_k()` + `trend_clusters()` (Mfuzz), enrichment (GO/Reactome) and hallmark mapping, GenAge concordance. |
| `01_metadata_expression.R` | Downloads expression series from GEO, parses sample metadata, extracts **age** and **sex** → expression metadata table. |
| `02_metadata_methylation.R` | Same for the methylation datasets → methylation metadata table. |
| `03_preprocess_expression.R` | Annotates probes to genes, collapses multi-probe genes, log2-transforms, kNN-imputes, intersects to common genes, merges all datasets → one gene × sample matrix. Predicts missing sex (XIST + Y-genes). |
| `04_preprocess_methylation.R` | Reads 450K betas, removes sex-chromosome / cross-reactive / high-NA probes, runs **EpiDISH** cell-type deconvolution, converts beta→M, keeps top-variable CpGs per dataset. |
| `04b_harmonize_methylation.R` | Re-reads betas once and rewrites every dataset to a **common 30 000-CpG set** (CpGs most variable across all datasets), so the matrices overlap for ComBat. Produces the sex-composition QC plot. |
| `05_deconv_expression.R` | Marker-gene z-score deconvolution (7 leukocyte types) → cell scores, composition heatmap, cell-score-vs-age plots. |
| `05_deconv_methylation.R` | Gathers the EpiDISH fractions from `04`, draws composition heatmap and cell-fraction-vs-age plots. |
| `06_combat_expression.R` | **ComBat** batch correction of the expression matrix (age, sex, cell scores preserved); PCA/UMAP before vs after. |
| `06_combat_methylation.R` | ComBat batch correction of the methylation M-values (same covariates). |
| `07_prepare_genage_expression.R` | Loads and formats the GenAge meta-analytic aging signature for use as an external validation reference. |
| `07_trends_expression.R` | Per-sex age-association testing (**limma + BH-FDR**), critical ages (GAM derivative + segmented breakpoints), trajectory clustering (Mfuzz), GenAge concordance, GO/Reactome enrichment. |
| `07_trends_methylation.R` | Same for methylation (limma on M-values; `missMethyl::gometh` enrichment with CpG-per-gene correction). |
| `08_integration.R` | The only step where the omics meet: gene-overlap **hypergeometric** test, direction concordance (**silenced** / **activated**), candidate drivers, GO/Reactome enrichment + hallmark mapping, critical-age overlay, transition-gene export. |
| `09_networks.R` | **STRING** protein–protein interaction networks (candidate and transition genes) with degree-centrality hubs, and **WGCNA** co-expression modules correlated with age. |

---

## Setup instructions

### 1. Input data (download from NCBI GEO)

Place series-matrix files in `data/raw/expression/` and the supplementary beta-value /
series-matrix files in `data/raw/methylation/`.

**Expression (9):** GSE14642, GSE30453, GSE47353, GSE49058, GSE56033, GSE56045, GSE56580, GSE65907, GSE68759
**Methylation (5):** GSE36054, GSE40279, GSE51032, GSE64495, GSE87571

Each dataset page is at `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSEXXXXX`, e.g.
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40279
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87571

(Most expression series are fetched automatically by `GEOquery` in script 01; large methylation
beta matrices are downloaded as supplementary files.)

### 2. Install dependencies

```r
source("R/install_packages.R")
```

### 3. Run the pipeline (from the project root)

```r
# metadata
source("scripts/01_metadata_expression.R")
source("scripts/02_metadata_methylation.R")
# preprocessing
source("scripts/03_preprocess_expression.R")
source("scripts/04_preprocess_methylation.R")
source("scripts/04b_harmonize_methylation.R")
# deconvolution
source("scripts/05_deconv_expression.R")
source("scripts/05_deconv_methylation.R")
# batch correction
source("scripts/06_combat_expression.R")
source("scripts/06_combat_methylation.R")
# trends, integration, networks
source("scripts/07_prepare_genage_expression.R")
source("scripts/07_trends_expression.R")
source("scripts/07_trends_methylation.R")
source("scripts/08_integration.R")
source("scripts/09_networks.R")
```

Outputs (tables and figures, png + svg) are written to `data/output/tables/` and `data/output/pics/`.

---

## Notes

- The design is **cross-sectional** (many cohorts of unrelated individuals): trends are population
  averages, not individual trajectories, and are subject to survivorship bias at older ages.
- Integration is performed at the **trend level**, not the sample level, because the expression and
  methylation cohorts are different people.

## Author

Done by Victor Florin, ITMO University, 2026. Supervisor: Aleksey Alekseev.
