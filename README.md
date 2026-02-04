# <img src="https://raw.githubusercontent.com/kuijjerlab/SCORPION/main/inst/logoSCORPION.png" width="30" title="SCORPION"> SCORPION

**SCORPION** (**S**ingle-**C**ell **O**riented **R**econstruction of **P**ANDA **I**ndividually **O**ptimized Gene Regulatory **N**etworks) is an R package that uses coarse-graining of single-cell/nuclei RNA-seq data to reduce sparsity and improve detection of the gene regulatory network's underlying correlation structure. The coarse-grained data is then used to reconstruct gene regulatory networks through the **PANDA** (**P**assing **A**ttributes between **N**etworks for **D**ata **A**ssimilation) message passing algorithm, which integrates protein-protein interactions, gene expression, and sequence motif data to predict regulatory relationships. Thanks to the use of the same baseline priors in each instance, this approach can reconstruct comparable, fully-connected, weighted, and directed transcriptome-wide single-cell gene regulatory networks suitable for use in population-level studies.

![method](https://raw.githubusercontent.com/kuijjerlab/SCORPION/main/inst/methodSCORPION.png)

---

## Table of Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [Example Data](#example-data)
- [Core Functions](#core-functions)
  - [`scorpion()` – Build a Single Network](#scorpion-function)
  - [`runSCORPION()` – Build Networks by Group](#runscopion-function)
- [Statistical Analysis](#statistical-analysis)
  - [`testEdges()` – Compare Network Edges](#testedges-function)
  - [`regressEdges()` – Trend Analysis](#regressedges-function)
- [Parameter Reference](#parameter-reference)
- [Citation](#citation)

---

## Installation

SCORPION is available through CRAN:

```r
install.packages("SCORPION")
library(SCORPION)
```

---

## Quick Start

```r
# Load example data
data(scorpionTest)

# Build a single network
network <- scorpion(
  tfMotifs = scorpionTest$tf,
  gexMatrix = scorpionTest$gex,
  ppiNet = scorpionTest$ppi
)

# Build networks for each group (e.g., tissue regions)
networks <- runSCORPION(
  gexMatrix = scorpionTest$gex,
  tfMotifs = scorpionTest$tf,
  ppiNet = scorpionTest$ppi,
  cellsMetadata = scorpionTest$metadata,
  groupBy = "region"
)
```

---

## Example Data

The package includes `scorpionTest`, a list with four objects:

| Object | Type | Description |
|--------|------|-------------|
| `gex` | dgCMatrix | Gene expression matrix (genes × cells) |
| `tf` | data.frame | TF-target motif pairs from dorothea (Human) |
| `ppi` | data.frame | Protein-protein interactions |
| `metadata` | data.frame | Cell annotations: `donor`, `region` (T/B/N), `cell_type` |

> **Region codes:** T = Tumor, B = Border, N = Normal

```r
data(scorpionTest)
str(scorpionTest)
```

<details>
<summary>View output structure</summary>

```
List of 4
 $ gex     :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
   ..@ Dim     : int [1:2] 300 1954
   ..@ Dimnames:List of 2
   .. ..$ : chr [1:300] "IGHM" "IGHG2" "IGLC3" ...
   .. ..$ : chr [1:1954] "P31-T_AAACGGGTCGGTTAAC" ...
 $ tf      :'data.frame': 371738 obs. of  3 variables:
   ..$ source_genesymbol: chr [1:371738] "MYC" "SPI1" ...
   ..$ target_genesymbol: chr [1:371738] "TERT" "BGLAP" ...
   ..$ weight           : num [1:371738] 1 1 1 ...
 $ ppi     :'data.frame': 4076 obs. of  3 variables:
   ..$ source_genesymbol: chr [1:4076] "ZIC1" "HES5" ...
   ..$ target_genesymbol: chr [1:4076] "ATOH1" "ATOH1" ...
   ..$ weight           : num [1:4076] 1 1 1 ...
 $ metadata:'data.frame': 1954 obs. of  4 variables:
   ..$ cell_id  : chr [1:1954] "P31-T_AAACGGGTCGGTTAAC" ...
   ..$ donor    : chr [1:1954] "P31" "P31" ...
   ..$ region   : chr [1:1954] "T" "T" ...
   ..$ cell_type: Factor w/ 1 level "Epithelial": 1 1 ...
```

</details>

---

## Core Functions

### `scorpion()` Function

Builds a single gene regulatory network using coarse-graining and the PANDA algorithm.

```r
result <- scorpion(
  tfMotifs = scorpionTest$tf,
  gexMatrix = scorpionTest$gex,
  ppiNet = scorpionTest$ppi,
  alphaValue = 0.1,      # Weight of prior networks (0-1)
  gammaValue = 10,       # Coarse-graining level
  nPC = 25,              # Principal components for kNN
  assocMethod = "pearson"
)
```

**Returns** a list containing:

| Component | Description |
|-----------|-------------|
| `regNet` | Regulatory network matrix (TF → target) |
| `coregNet` | Co-regulation network (gene-gene) |
| `coopNet` | Cooperation network (TF-TF) |
| `numGenes` | Number of genes |
| `numTFs` | Number of transcription factors |
| `numEdges` | Total edges in regulatory network |

---

### `runSCORPION()` Function

Builds separate networks for cell groups defined by metadata columns, returning a wide-format data frame suitable for population-level comparisons.

#### Basic Usage

```r
# Group by single column
nets_by_region <- runSCORPION(
  gexMatrix = scorpionTest$gex,
  tfMotifs = scorpionTest$tf,
  ppiNet = scorpionTest$ppi,
  cellsMetadata = scorpionTest$metadata,
  groupBy = "region"
)

# Group by multiple columns (creates donor-region combinations)
nets_by_donor_region <- runSCORPION(
  gexMatrix = scorpionTest$gex,
  tfMotifs = scorpionTest$tf,
  ppiNet = scorpionTest$ppi,
  cellsMetadata = scorpionTest$metadata,
  groupBy = c("donor", "region")
)
```

**Output format:**

| tf | target | P31--T | P31--B | P31--N | P32--T | ... |
|----|--------|--------|--------|--------|--------|-----|
| AATF | ACKR1 | -0.326 | -0.337 | -0.344 | -0.298 | ... |
| ABL1 | ACKR1 | -0.340 | -0.339 | -0.351 | -0.312 | ... |

#### With Batch Correction

```r
nets_corrected <- runSCORPION(
  gexMatrix = scorpionTest$gex,
  tfMotifs = scorpionTest$tf,
  ppiNet = scorpionTest$ppi,
  cellsMetadata = scorpionTest$metadata,
  groupBy = "region",
  removeBatchEffect = TRUE,
  batch = scorpionTest$metadata$donor
)
```

#### With GPU Acceleration

```r
nets_gpu <- runSCORPION(
  gexMatrix = scorpionTest$gex,
  tfMotifs = scorpionTest$tf,
  ppiNet = scorpionTest$ppi,
  cellsMetadata = scorpionTest$metadata,
  groupBy = "region",
  computingEngine = "gpu"
)
```

#### Custom PANDA Parameters

```r
nets_custom <- runSCORPION(
  gexMatrix = scorpionTest$gex,
  tfMotifs = scorpionTest$tf,
  ppiNet = scorpionTest$ppi,
  cellsMetadata = scorpionTest$metadata,
  groupBy = "region",
  alphaValue = 0.5,         # Stronger prior influence
  nPC = 30,                 # More principal components
  assocMethod = "spearman"  # Rank-based correlation
)
```

---

## Statistical Analysis

### `testEdges()` Function

Performs statistical testing on network edges to identify differential regulatory relationships between groups. Uses fully vectorized computations for efficiency with large-scale datasets (millions of edges).

**Supports three test types:**

| Test Type | Use Case | Method |
|-----------|----------|--------|
| Single-sample | Test if edges differ from zero | One-sample t-test |
| Two-sample | Compare edges between groups | Welch's t-test |
| Paired | Compare matched samples | Paired t-test |

#### Two-Sample Comparison

```r
# Build networks
nets <- runSCORPION(
  gexMatrix = scorpionTest$gex,
  tfMotifs = scorpionTest$tf,
  ppiNet = scorpionTest$ppi,
  cellsMetadata = scorpionTest$metadata,
  groupBy = c("donor", "region")
)

# Select groups
tumor_nets <- grep("--T$", colnames(nets), value = TRUE)
normal_nets <- grep("--N$", colnames(nets), value = TRUE)

# Compare Tumor vs Normal
results <- testEdges(
  networksDF = nets,
  testType = "two.sample",
  group1 = tumor_nets,
  group2 = normal_nets
)

# Top differential edges
head(results[order(results$pAdj), ])
```

#### Paired T-Test (Matched Samples)

When samples are matched (e.g., tumor and normal from the same patient), paired tests are more powerful:

```r
# Ensure columns are ordered by patient
tumor_ordered <- c("P31--T", "P32--T", "P33--T")
normal_ordered <- c("P31--N", "P32--N", "P33--N")

results_paired <- testEdges(
  networksDF = nets,
  testType = "two.sample",
  group1 = tumor_ordered,
  group2 = normal_ordered,
  paired = TRUE
)
```

#### Single-Sample Test

```r
# Test if tumor edges differ from zero
results_single <- testEdges(
  networksDF = nets,
  testType = "single",
  group1 = tumor_nets
)
```

**Output columns:**

| Column | Description |
|--------|-------------|
| `tf`, `target` | TF-target pair identifiers |
| `meanGroup1`, `meanGroup2` | Group means (two-sample only) |
| `diffMean` | Difference in means (two-sample only) |
| `log2FoldChange` | Log2 fold change (two-sample only) |
| `tStatistic` | Test statistic |
| `pValue` | Raw p-value |
| `pAdj` | Adjusted p-value (FDR) |

---

### `regressEdges()` Function

Performs linear regression to identify edges with significant trends across ordered conditions (e.g., disease progression: Normal → Border → Tumor).

```r
# Define progression: Normal → Border → Tumor
ordered_conditions <- list(
  Normal = grep("--N$", colnames(nets), value = TRUE),
  Border = grep("--B$", colnames(nets), value = TRUE),
  Tumor = grep("--T$", colnames(nets), value = TRUE)
)

# Analyze trends
results_reg <- regressEdges(
  networksDF = nets,
  orderedGroups = ordered_conditions
)

# Find edges increasing along progression
increasing <- results_reg[results_reg$pAdj < 0.05 & results_reg$slope > 0, ]

# Find edges decreasing along progression
decreasing <- results_reg[results_reg$pAdj < 0.05 & results_reg$slope < 0, ]
```

**Output columns:**

| Column | Description |
|--------|-------------|
| `slope` | Change in edge weight per condition step |
| `intercept` | Regression intercept |
| `rSquared` | Model fit quality (0-1) |
| `fStatistic` | F-statistic for regression |
| `pValue`, `pAdj` | Raw and adjusted p-values |
| `meanNormal`, `meanBorder`, `meanTumor` | Condition-specific means |

---

## Parameter Reference

### Data Preprocessing

| Parameter | Description | Default |
|-----------|-------------|---------|
| `normalizeData` | Apply log normalization | `TRUE` |
| `removeBatchEffect` | Correct for batch effects | `FALSE` |
| `batch` | Batch assignment vector | — |
| `filterExpr` | Remove zero-expression genes | `TRUE` |

### Network Construction

| Parameter | Description | Default |
|-----------|-------------|---------|
| `gammaValue` | Coarse-graining level (higher = more cells/super-cell) | `10` |
| `nPC` | Principal components for kNN | `25` |
| `assocMethod` | Correlation method: `pearson`, `spearman`, `pcNet` | `pearson` |

### PANDA Algorithm

| Parameter | Description | Default |
|-----------|-------------|---------|
| `alphaValue` | Prior network weight (0-1) | `0.1` |
| `hammingValue` | Convergence threshold | `0.001` |
| `nIter` | Maximum iterations | `Inf` |
| `zScaling` | Output as Z-scores | `TRUE` |

### Computing

| Parameter | Description | Default |
|-----------|-------------|---------|
| `computingEngine` | `cpu` or `gpu` | `cpu` |
| `nCores` | Processors for BLAS/MPI | `1` |

### Output Options

| Parameter | Description | Default |
|-----------|-------------|---------|
| `outNet` | Network type: `regNet`, `coregNet`, `coopNet` | `regNet` |
| `scaleByPresent` | Scale by % positive samples | `FALSE` |

---

## Citation

Osorio, D., Capasso, A., Eckhardt, S.G. et al. Population-level comparisons of gene regulatory networks modeled on high-throughput single-cell transcriptomics data. *Nat Comput Sci* **4**, 237–250 (2024). https://doi.org/10.1038/s43588-024-00597-5

---

## Supplementary Information

For supplementary materials, visit: [SCORPION Supplementary Information](https://github.com/dosorio/SCORPION/)
