# <img src="https://raw.githubusercontent.com/kuijjerlab/SCORPION/main/inst/logoSCORPION.png" width="30" title="SCORPION"> SCORPION

**SCORPION** (<ins>S</ins>ingle-<ins>C</ins>ell <ins>O</ins>riented <ins>R</ins>econstruction of <ins>P</ins>ANDA <ins>I</ins>ndividually <ins>O</ins>ptimized Gene Regulatory <ins>N</ins>etworks), is an R package that uses coarse-graining of single-cell/nuclei RNA-seq data to reduce sparsity and improve the ability to detect the gene regulatory network's underlying correlation structure. The coarse-grained data generated is then used to reconstruct the gene regulatory network using a network refinement strategy through the **PANDA** (<ins>P</ins>assing <ins>A</ins>ttributes between <ins>N</ins>etworks for <ins>D</ins>ata <ins>A</ins>ssimilation) message passing algorithm. This algorithm is designed to integrate multiple sources of information such as protein-protein interaction, gene expression, and sequence motif data to predict accurate regulatory relationships. Thanks to the use of the same baseline priors in each instance, this approach can reconstruct comparable, fully-connected, weighted, and directed transcriptome-wide single-cell gene regulatory networks suitable for use in population-level studies.
## Method
![method](https://raw.githubusercontent.com/kuijjerlab/SCORPION/main/inst/methodSCORPION.png)

## Usage
**SCORPION** is available through the CRAN repositories, you can install and load it, using the following command:
```{r}
install.packages('SCORPION')
library(SCORPION)
```
## Example Data
We provide an example dataset (formally a ```list```) containing four objects. The motif ```data.frame``` describes a set of pairwise connections where a specific known sequence motif of a transcription factor was found upstream of the corresponding gene. For this particular example, the data is a subset of the transcription-factor and target gene pairs provided by the ```dorothea``` package for *Homo sapiens*. The expression ```dgCMatrix``` contains gene expression levels measured across single cells. The ppi ```data.frame``` describes a set of known pairwise protein-protein interactions. Finally, the metadata ```data.frame``` contains cell-level information including donor, region, and cell_type annotations.
```{R}
data(scorpionTest)
```
The structure of the data can be accessed as follows:
```{R}
str(scorpionTest)

# List of 4
# $ gex     :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
# .. ..@ i       : int [1:46171] 29 32 41 43 61 170 208 245 251 269 ...
# .. ..@ p       : int [1:1955] 0 11 62 97 112 163 184 215 257 274 ...
# .. ..@ Dim     : int [1:2] 300 1954
# .. ..@ Dimnames:List of 2
# .. .. ..$ : chr [1:300] "IGHM" "IGHG2" "IGLC3" "IGLL5" ...
# .. .. ..$ : chr [1:1954] "P31-T_AAACGGGTCGGTTAAC" "P31-T_AAAGATGGTGGCCCTA" ...
# .. ..@ x       : num [1:46171] 1 1 1 1 2 2 1 1 2 1 ...
# .. ..@ factors : list()
# $ tf      :'data.frame':	371738 obs. of  3 variables:
# ..$ source_genesymbol: chr [1:371738] "MYC" "SPI1" "JUN_JUND" "FOS_JUND" ...
# ..$ target_genesymbol: chr [1:371738] "TERT" "BGLAP" "JUN" "JUN" ...
# ..$ weight           : num [1:371738] 1 1 1 1 1 1 1 1 1 1 ...
# $ ppi     :'data.frame':	4076 obs. of  3 variables:
# ..$ source_genesymbol: chr [1:4076] "ZIC1" "HES5" "ATOH1" "DLL1" ...
# ..$ target_genesymbol: chr [1:4076] "ATOH1" "ATOH1" "HES5" "NOTCH1" ...
# ..$ weight           : num [1:4076] 1 1 1 1 1 1 1 1 1 1 ...
# $ metadata:'data.frame':	1954 obs. of  4 variables:
# ..$ cell_id  : chr [1:1954] "P31-T_AAACGGGTCGGTTAAC" "P31-T_AAAGATGGTGGCCCTA"...
# ..$ donor    : chr [1:1954] "P31" "P31" "P31" "P31" ...
# ..$ region   : chr [1:1954] "T" "T" "T" "T" ...
# ..$ cell_type: Factor w/ 1 level "Epithelial": 1 1 1 1 1 1 1 1 1 1 ...
```

### Using the `scorpion()` function
The `scorpion()` function builds a single gene regulatory network by applying coarse-graining to reduce sparsity and then running the PANDA (Passing Attributes between Networks for Data Assimilation) algorithm to integrate transcription factor motifs, protein-protein interactions, and gene expression data.

**Key parameters:**
- `gexMatrix`: Gene expression matrix with genes in rows and cells in columns
- `tfMotifs`: TF-target motifs with columns [TF, target gene, motif score]
- `ppiNet`: Protein-protein interactions with columns [protein 1, protein 2, interaction score]
- `alphaValue`: Update parameter (0-1) controlling relative contribution of prior networks. Default 0.1
- `gammaValue`: Coarse-graining level (default 10). Higher values = more cells per super-cell
- `nPC`: Principal components for kNN network construction (default 25)
- `assocMethod`: Association method - 'pearson', 'spearman', or 'pcNet' (default 'pearson')
- `computingEngine`: 'cpu' or 'gpu' for computation

Here, we are running SCORPION with `alphaValue = 0.8` for testing purposes. The default value is `0.1`.
```{R}
scorpionOutput <- scorpion(tfMotifs = scorpionTest$tf,
                           gexMatrix = scorpionTest$gex,
                           ppiNet = scorpionTest$ppi,
                           alphaValue = 0.8)
── SCORPION ──────────────────────────────────────────────
✔ Initializing and validating
✔ Verified sufficient samples
ℹ Normalizing networks
ℹ Learning Network
ℹ Using tanimoto similarity
✔ Successfully ran SCORPION on 281 Genes and 963 TFs
ℹ Time elapsed: 2.12 seconds            
```

The structure of the output can be accessed as follows:
```{R}
str(scorpionOutput)

# List of 6
# $ regNet  : num [1:963, 1:281] -0.1556 -0.0455 -0.1461 1.6881 0.8746 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:963] "AATF" "ABL1" "ACSS2" "ADNP" ...
# .. ..$ : chr [1:281] "ACKR1" "ACTA2" "ACTG2" "ADAMDEC1" ...
# $ coregNet: num [1:281, 1:281] 2.02e+06 3.84 4.10 -1.26 8.81e-01 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:281] "ACKR1" "ACTA2" "ACTG2" "ADAMDEC1" ...
# .. ..$ : chr [1:281] "ACKR1" "ACTA2" "ACTG2" "ADAMDEC1" ...
# $ coopNet : num [1:963, 1:963] 1.17e+07 -2.66 8.13 -1.31 4.95 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:963] "AATF" "ABL1" "ACSS2" "ADNP" ...
# .. ..$ : chr [1:963] "AATF" "ABL1" "ACSS2" "ADNP" ...
# $ numGenes: int 281
# $ numTFs  : int 963
# $ numEdges: int 270603
```

**Output components:**
- `regNet`: Regulatory network (TF → target gene interactions)
- `coregNet`: Co-regulation network (gene-gene co-regulation)
- `coopNet`: Cooperation network (TF-TF interactions)
- `numGenes`: Number of genes in the network
- `numTFs`: Number of transcription factors
- `numEdges`: Total edges in the regulatory network

### Using the `runSCORPION()` function for grouped analysis
The `runSCORPION()` function extends SCORPION by allowing you to build separate gene regulatory networks for cell groups defined by metadata columns. This wrapper automatically handles data grouping, filtering (keeping groups with ≥30 cells by default), and combining results into a wide-format data frame suitable for population-level comparisons.

**Key parameters:**
- `countMatrix`: Gene expression matrix (genes × cells)
- `tfTargets`: TF-target motifs with columns [TF, target gene, score]
- `ppiNet`: Protein-protein interactions with columns [protein 1, protein 2, score]
- `cellsMetadata`: Cell metadata data.frame containing grouping variables
- `groupBy`: Column name(s) in `cellsMetadata` for grouping cells (creates separate networks for each group)
- `minCells`: Minimum cells required per group (default 30)
- `normalizeData`: Boolean for log normalization (default TRUE)
- `removeBatchEffect`: Boolean for batch correction (default FALSE)
- `batch`: Batch assignment vector (required if `removeBatchEffect = TRUE`)
- `computingEngine`: 'cpu' or 'gpu' for computation (default 'cpu')
- Additional PANDA parameters: `alphaValue`, `nPC`, `assocMethod`, `gammaValue`, `hammingValue`, etc.

**Output:** Wide-format data.frame where rows are TF-target pairs and columns are network identifiers. Cell values contain edge weights from each network.

#### Example 1: Group by single column
```{R}
# Build networks grouped by tissue region
nets_by_region <- runSCORPION(
  countMatrix = scorpionTest$gex,
  tfTargets = scorpionTest$tf,
  ppiNet = scorpionTest$ppi,
  cellsMetadata = scorpionTest$metadata,
  groupBy = "region"
)
# ── SCORPION ────────────────────────────────────────────────────────────────
# ✔ Normalizing data (log scale)
# ℹ 3 networks requested
# ✔ 3 networks meet the minimum cell requirement (30)
# ℹ Computing networks
# ✔ Networks successfully constructed
# ✔ Networks successfully combined
# 
# head(nets_by_region)
#                           tf target           T          B           N
# 1                       AATF  ACKR1 -0.31433856 -0.3569918 -0.33734920
# 2                       ABL1  ACKR1 -0.32915008 -0.3648895 -0.34437341
# 3                      ACSS2  ACKR1 -0.31418599 -0.3557854 -0.33663144
# 4                       ADNP  ACKR1  0.04105895  0.1109288  0.09910822
# 5                      AEBP2  ACKR1 -0.18964574 -0.2202269 -0.17558140
# 6 AEBP2_EED_EZH2_RBBP4_SUZ12  ACKR1 -0.31024700 -0.3508320 -0.33054519
```

#### Example 2: Group by multiple columns
```{R}
# Build networks grouped by both donor and region
nets_by_donor_region <- runSCORPION(
  countMatrix = scorpionTest$gex,
  tfTargets = scorpionTest$tf,
  ppiNet = scorpionTest$ppi,
  cellsMetadata = scorpionTest$metadata,
  groupBy = c("donor", "region")
)
# ── SCORPION ────────────────────────────────────────────────────────────────
# ✔ Normalizing data (log scale)
# ℹ 9 networks requested
# ✔ 9 networks meet the minimum cell requirement (30)
# ℹ Computing networks
# ✔ Networks successfully constructed
# ✔ Networks successfully combined
#
# head(nets_by_donor_region)
#                           tf target      P31--T      P31--B     P31--N
# 1                       AATF  ACKR1 -0.32634975 -0.33717677 -0.3442886
# 2                       ABL1  ACKR1 -0.34048759 -0.33890429 -0.3509986
# 3                      ACSS2  ACKR1 -0.32570697 -0.33600811 -0.3436603
# 4                       ADNP  ACKR1  0.07975735  0.05354279  0.1048301
# 5                      AEBP2  ACKR1 -0.21472437 -0.20545660 -0.1815737
# 6 AEBP2_EED_EZH2_RBBP4_SUZ12  ACKR1 -0.31861592 -0.32809314 -0.3375652
```

#### Example 3: Enable batch effect correction
```{R}
# Build networks with batch effect correction using donor as batch
nets_batch_corrected <- runSCORPION(
  countMatrix = scorpionTest$gex,
  tfTargets = scorpionTest$tf,
  ppiNet = scorpionTest$ppi,
  cellsMetadata = scorpionTest$metadata,
  groupBy = "region",
  removeBatchEffect = TRUE,
  batch = scorpionTest$metadata$donor
)
# ── SCORPION ────────────────────────────────────────────────────────────────
# ✔ Normalizing data (log scale)
# ✔ Correcting for batch effects
# ℹ 3 networks requested
# ✔ 3 networks meet the minimum cell requirement (30)
# ℹ Computing networks
# ✔ Networks successfully constructed
# ✔ Networks successfully combined
#
# head(nets_batch_corrected)
#                           tf target          T           B           N
# 1                       AATF  ACKR1 -0.3337298 -0.34885471 -0.13011777
# 2                       ABL1  ACKR1 -0.3408020 -0.35409813 -0.17694266
# 3                      ACSS2  ACKR1 -0.3325270 -0.35115311 -0.12661518
# 4                       ADNP  ACKR1  0.1117504  0.08691481  0.01608898
# 5                      AEBP2  ACKR1 -0.2334648 -0.22113011  0.12519312
# 6 AEBP2_EED_EZH2_RBBP4_SUZ12  ACKR1 -0.3274770 -0.34475499 -0.12449908
```

#### Example 4: Use GPU computing engine
```{R}
# Build networks using GPU acceleration (if available)
nets_gpu <- runSCORPION(
  countMatrix = scorpionTest$gex,
  tfTargets = scorpionTest$tf,
  ppiNet = scorpionTest$ppi,
  cellsMetadata = scorpionTest$metadata,
  groupBy = "region",
  computingEngine = "gpu"
)
# ── SCORPION ────────────────────────────────────────────────────────────────
# ✔ Normalizing data (log scale)
# ℹ 3 networks requested
# ✔ 3 networks meet the minimum cell requirement (30)
# ℹ Computing networks
# ✔ Networks successfully constructed
# ✔ Networks successfully combined
#
# head(nets_gpu)
#                           tf target           T          B           N
# 1                       AATF  ACKR1 -0.31433821 -0.3569913 -0.33734894
# 2                       ABL1  ACKR1 -0.32915005 -0.3648892 -0.34437302
# 3                      ACSS2  ACKR1 -0.31418574 -0.3557851 -0.33663106
# 4                       ADNP  ACKR1  0.04105883  0.1109285  0.09910798
# 5                      AEBP2  ACKR1 -0.18964562 -0.2202267 -0.17558131
# 6 AEBP2_EED_EZH2_RBBP4_SUZ12  ACKR1 -0.31024694 -0.3508317 -0.33054504
```

#### Example 5: Customize PANDA parameters
```{R}
# Build networks with custom PANDA parameters
nets_custom <- runSCORPION(
  countMatrix = scorpionTest$gex,
  tfTargets = scorpionTest$tf,
  ppiNet = scorpionTest$ppi,
  cellsMetadata = scorpionTest$metadata,
  groupBy = "region",
  alphaValue = 0.5,        # Increase weight of prior networks
  nPC = 30,                # More PCs for kNN network
  assocMethod = "spearman" # Use rank correlation
)
# ── SCORPION ────────────────────────────────────────────────────────────────
# ✔ Normalizing data (log scale)
# ℹ 3 networks requested
# ✔ 3 networks meet the minimum cell requirement (30)
# ℹ Computing networks
# ✔ Networks successfully constructed
# ✔ Networks successfully combined
#
# head(nets_custom)
#                           tf target           T          B           N
# 1                       AATF  ACKR1  0.14665281 -0.3813583 -0.13632369
# 2                       ABL1  ACKR1 -0.04234458 -0.1687202 -0.13003033
# 3                      ACSS2  ACKR1  0.16172807 -0.3757015 -0.13184588
# 4                       ADNP  ACKR1  0.37028678  0.9410779  0.76758912
# 5                      AEBP2  ACKR1 -0.08311113 -0.4529925 -0.07731317
# 6 AEBP2_EED_EZH2_RBBP4_SUZ12  ACKR1  0.15790771 -0.3597740 -0.11774763
```

## Advanced Usage

### Parameter Guide

**Data Normalization and Preprocessing:**
- `normalizeData = TRUE`: Applies log normalization to count matrix
- `removeBatchEffect = TRUE` + `batch`: Corrects for batch effects before network inference
- `filterExpr = TRUE`: Removes genes with zero expression across all cells

**Network Construction:**
- `gammaValue`: Controls coarse-graining level (default 10). Lower values = finer graining
- `nPC`: Principal components for kNN construction (default 25). Higher values capture more variance
- `assocMethod`: Choose 'pearson' (default), 'spearman' (rank-based), or 'pcNet' (principal component)

**PANDA Algorithm:**
- `alphaValue`: (0-1) Relative weight of prior networks (default 0.1). Higher values = stronger prior influence
- `hammingValue`: Convergence threshold based on Hamming distance (default 0.001)
- `nIter`: Maximum iterations before stopping (default Inf)
- `zScaling = TRUE`: Output as Z-scores; `FALSE` for [0,1] scale

**Computing:**
- `computingEngine`: 'cpu' (default) or 'gpu' for GPU acceleration
- `nCores`: Number of processors for BLAS/MPI
- `randomizationMethod`: 'None' (default), 'within.gene', or 'by.gene' for shuffled controls

**Output:**
- `outNet`: Which network(s) to extract: 'regNet' (regulatory, default), 'coregNet' (co-regulation), 'coopNet' (cooperation)
- `scaleByPresent`: Boolean to scale correlations by percentage of positive samples

### Understanding the Output

For `scorpion()`: Returns a list with three network matrices (regNet, coregNet, coopNet) and summary statistics.

For `runSCORPION()`: Returns wide-format data.frame where:
- **Rows**: Unique TF-target pairs (union across all groups)
- **Columns**: Network identifiers (combinations of `groupBy` variables)
- **Values**: Edge weights from corresponding networks, allowing direct comparison across groups

## Citation
Please cite: Osorio, D., Capasso, A., Eckhardt, S.G. et al. Population-level comparisons of gene regulatory networks modeled on high-throughput single-cell transcriptomics data. *Nat Comput Sci* **4**, 237–250 (2024). https://doi.org/10.1038/s43588-024-00597-5

## Supplementary Information
For the supplementary information of *Population-level comparisons of gene regulatory networks modeled on high-throughput single-cell transcriptomics data*. Please visit: [SCORPION Supplementary Information](https://github.com/dosorio/SCORPION/)