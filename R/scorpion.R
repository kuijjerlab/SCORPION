#' @importFrom Matrix Matrix
#' @export scorpion
#' @title Build gene regulatory networks from single-cell RNA-seq data using PANDA
#' @description
#' Constructs gene regulatory networks from single-cell/nuclei RNA-seq data by first
#' applying coarse-graining to reduce sparsity, then running the PANDA (Passing
#' Attributes between Networks for Data Assimilation) message-passing algorithm to
#' integrate transcription factor motifs, protein-protein interactions, and gene
#' expression data into unified regulatory networks.
#' @author Daniel Osorio <daniecos@uio.no>
#' @param tfMotifs A motif dataset (data.frame or matrix) with 3 columns: TF, target gene, and motif score. Pass NULL for co-expression analysis only.
#' @param gexMatrix An expression dataset, with genes in the rows and barcodes (cells) in the columns.
#' @param ppiNet A Protein-Protein-Interaction dataset (data.frame or matrix) with 3 columns: protein 1, protein 2, and interaction score. Pass NULL to disable protein interaction integration.
#' @param computingEngine Character specifying computing device: 'cpu' or 'gpu' (if available). Default 'cpu'.
#' @param nCores Number of processors to be used if BLAS or MPI is active.
#' @param gammaValue Graining level of data (proportion of number of single cells in the initial dataset to the number of super-cells in the final dataset)
#' @param nPC Number of principal components to use for construction of single-cell kNN network.
#' @param alphaValue Numeric update parameter (0 to 1) controlling relative contribution of prior networks. Default 0.1.
#' @param hammingValue Numeric convergence threshold based on Hamming distance. Algorithm stops when updates fall below this. Default 0.001.
#' @param assocMethod Association method. Must be one of 'pearson', 'spearman' or 'pcNet'.
#' @param nIter Sets the maximum number of iterations PANDA can run before exiting.
#' @param outNet A vector containing which networks to return. Options include "regNet", "coregNet", "coopNet".
#' @param zScaling Boolean to indicate use of Z-Scores in output. FALSE will use [0,1] scale.
#' @param showProgress Boolean to indicate printing of output for algorithm progress.
#' @param randomizationMethod Method by which to randomize gene expression matrix. Default "None". Must be one of "None", "within.gene", "by.genes". "within.gene" randomization scrambles each row of the gene expression matrix, "by.gene" scrambles gene labels.
#' @param scaleByPresent Boolean to indicate scaling of correlations by percentage of positive samples.
#' @param filterExpr Boolean to remove genes with zero expression across all cells
#'   before network inference. Default FALSE.
#' @return A list of 6 elements describing the inferred networks at convergence:
#'   \itemize{
#'     \item{regNet: Regulatory network matrix (TFs × genes)}
#'     \item{coregNet: Co-regulation network matrix (genes × genes)}
#'     \item{coopNet: Cooperation network matrix (TFs × TFs)}
#'     \item{numGenes: Number of genes in the network}
#'     \item{numTFs: Number of transcription factors}
#'     \item{numEdges: Total number of edges in regulatory network}
#'   }
#' @seealso \code{\link{runSCORPION}} for building networks across cell groups.
#' @examples
#' # Loading example data
#' data(scorpionTest)
#'
#' # The structure of the data
#' str(scorpionTest)
#'
#' # List of 4
#' # $ gex     :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#' # .. ..@ i       : int [1:46171] 29 32 41 43 61 170 208 245 251 269 ...
#' # .. ..@ p       : int [1:1955] 0 11 62 97 112 163 184 215 257 274 ...
#' # .. ..@ Dim     : int [1:2] 300 1954
#' # .. ..@ Dimnames:List of 2
#' # .. .. ..$ : chr [1:300] "IGHM" "IGHG2" "IGLC3" "IGLL5" ...
#' # .. .. ..$ : chr [1:1954] "P31-T_AAACGGGTCGGTTAAC" "P31-T_AAAGATGGTGGCCCTA" ...
#' # .. ..@ x       : num [1:46171] 1 1 1 1 2 2 1 1 2 1 ...
#' # .. ..@ factors : list()
#' # $ tf      :'data.frame':	371738 obs. of  3 variables:
#' # ..$ source_genesymbol: chr [1:371738] "MYC" "SPI1" "JUN_JUND" "FOS_JUND" ...
#' # ..$ target_genesymbol: chr [1:371738] "TERT" "BGLAP" "JUN" "JUN" ...
#' # ..$ weight           : num [1:371738] 1 1 1 1 1 1 1 1 1 1 ...
#' # ..- attr(*, "origin")= chr "cache"
#' # ..- attr(*, "url")= chr "https://omnipathdb.org/interactions? __truncated__
#' # $ ppi     :'data.frame':	4076 obs. of  3 variables:
#' # ..$ source_genesymbol: chr [1:4076] "ZIC1" "HES5" "ATOH1" "DLL1" ...
#' # ..$ target_genesymbol: chr [1:4076] "ATOH1" "ATOH1" "HES5" "NOTCH1" ...
#' # ..$ weight           : num [1:4076] 1 1 1 1 1 1 1 1 1 1 ...
#' # ..- attr(*, "origin")= chr "cache"
#' # ..- attr(*, "url")= chr "https://omnipathdb.org/interactions?__truncated__
#' # $ metadata:'data.frame':	1954 obs. of  4 variables:
#' # ..$ cell_id  : chr [1:1954] "P31-T_AAACGGGTCGGTTAAC" "P31-T_AAAGATGGTGGCCCTA"...
#' # ..$ donor    : chr [1:1954] "P31" "P31" "P31" "P31" ...
#' # ..$ region   : chr [1:1954] "T" "T" "T" "T" ...
#' # ..$ cell_type: Factor w/ 1 level "Epithelial": 1 1 1 1 1 1 1 1 1 1 ...
#'
#' # Running SCORPION for epithelial cells from the normal tissue
#' # We are using alphaValue = 0.8 for testing purposes (Default = 0.1).
#' scorpionOutput <- scorpion(
#'   tfMotifs = scorpionTest$tf,
#'   gexMatrix = scorpionTest$gex[, scorpionTest$metadata$region == "N"],
#'   ppiNet = scorpionTest$ppi,
#'   alphaValue = 0.8
#' )
#'
#' # -- SCORPION --------------------------------------------------------------------------------------
#' # + Initializing and validating
#' # + Verified sufficient samples
#' # i Normalizing networks
#' # i Learning Network
#' # i Using tanimoto similarity
#' # + Successfully ran SCORPION on 281 Genes and 963 TFs
#'
#' # Structure of the output.
#' str(scorpionOutput)
#'
#' # List of 6
#' # $ regNet  : num [1:963, 1:281] -0.1556 -0.0455 -0.1461 1.6881 0.8746 ...
#' # ..- attr(*, "dimnames")=List of 2
#' # .. ..$ : chr [1:963] "AATF" "ABL1" "ACSS2" "ADNP" ...
#' # .. ..$ : chr [1:281] "ACKR1" "ACTA2" "ACTG2" "ADAMDEC1" ...
#' # $ coregNet: num [1:281, 1:281] 2.02e+06 3.84 4.10 -1.26 8.81e-01 ...
#' # ..- attr(*, "dimnames")=List of 2
#' # .. ..$ : chr [1:281] "ACKR1" "ACTA2" "ACTG2" "ADAMDEC1" ...
#' # .. ..$ : chr [1:281] "ACKR1" "ACTA2" "ACTG2" "ADAMDEC1" ...
#' # $ coopNet : num [1:963, 1:963] 1.17e+07 -2.66 8.13 -1.31 4.95 ...
#' # ..- attr(*, "dimnames")=List of 2
#' # .. ..$ : chr [1:963] "AATF" "ABL1" "ACSS2" "ADNP" ...
#' # .. ..$ : chr [1:963] "AATF" "ABL1" "ACSS2" "ADNP" ...
#' # $ numGenes: int 281
#' # $ numTFs  : int 963
#' # $ numEdges: int 270603
scorpion <- function(tfMotifs = NULL,
                     gexMatrix,
                     ppiNet = NULL,
                     computingEngine = "cpu",
                     nCores = 1,
                     gammaValue = 10,
                     nPC = 25,
                     assocMethod = "pearson",
                     alphaValue = 0.1,
                     hammingValue = 0.001,
                     nIter = Inf,
                     outNet = c("regNet", "coregNet", "coopNet"),
                     zScaling = TRUE,
                     showProgress = TRUE,
                     randomizationMethod = "None",
                     scaleByPresent = FALSE,
                     filterExpr = FALSE) {
  # Input validation
  if (!is.numeric(alphaValue) || alphaValue < 0 || alphaValue > 1) {
    cli::cli_abort("alphaValue must be a numeric value between 0 and 1")
  }
  if (!is.numeric(hammingValue) || hammingValue < 0) {
    cli::cli_abort("hammingValue must be a non-negative numeric value")
  }
  if (!assocMethod %in% c("pearson", "spearman", "pcNet")) {
    cli::cli_abort("assocMethod must be one of: 'pearson', 'spearman', 'pcNet'")
  }
  if (!computingEngine %in% c("cpu", "gpu")) {
    cli::cli_abort("computingEngine must be either 'cpu' or 'gpu'")
  }
  if (ncol(gexMatrix) < 30) {
    cli::cli_warn("gexMatrix has fewer than 30 cells. Network inference may be unreliable.")
  }
  if (nrow(gexMatrix) < 10) {
    cli::cli_abort("gexMatrix must have at least 10 genes (rows)")
  }

  if (showProgress) {
    cli::cli_h1("SCORPION")
  }

  if (isTRUE(filterExpr)) {
    gexMatrix <- gexMatrix[rowSums(gexMatrix) > 0, ]
  }
  gexMatrix <- makeSuperCells(X = gexMatrix, gamma = gammaValue, n.pc = nPC, fast.pca = FALSE)

  if (is.null(ppiNet) & is.null(tfMotifs)) {
    if (assocMethod == "spearman") {
      gexMatrix <- Matrix(t(apply(gexMatrix, 1, rank)))
    }
    geneCoExpr <- gexMatrix - rowMeans(gexMatrix)
    geneCoExpr <- geneCoExpr / sqrt(rowSums(geneCoExpr^2))
    geneCoExpr <- tcrossprod(geneCoExpr)
    return(geneCoExpr)
  }
  outNetworks <- runPANDA(
    motif = tfMotifs,
    ppi = ppiNet,
    expr = gexMatrix,
    computing.engine = computingEngine,
    n.cores = nCores,
    alpha = alphaValue,
    hamming = hammingValue,
    iter = nIter,
    output = outNet,
    zScale = zScaling,
    progress = showProgress,
    randomize = randomizationMethod,
    assoc.method = assocMethod,
    scale.by.present = scaleByPresent
  )
  return(outNetworks)
}
