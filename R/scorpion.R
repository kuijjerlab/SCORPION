#' @importFrom Matrix Matrix
#' @export scorpion
#' @title Build gene regulatory networks from single-cell RNA-seq data using PANDA
#' @description
#' Constructs gene regulatory networks from single-cell/nuclei RNA-seq data by first
#' applying coarse-graining to reduce sparsity, then running the PANDA (Passing
#' Attributes between Networks for Data Assimilation) message-passing algorithm to
#' integrate transcription factor motifs, protein-protein interactions, and gene
#' expression data into unified regulatory networks. When \code{"regSign"} is requested
#' (the default), the regulatory network is also returned with a direction: each edge is
#' signed as activating or repressing using \code{\link{inferEdgeSign}}, computed from the
#' same cells after PANDA, so \code{regNet} is identical with or without it.
#' @author Daniel Osorio <daniecos@uio.no>
#' @param tfMotifs A motif dataset (data.frame or matrix) with 3 columns: TF, target gene, and motif score, and optionally a fourth column \code{mor} (mode of regulation: 1 activation, -1 repression, 0 unknown) used only for \code{regSign}. Pass NULL for co-expression analysis only.
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
#' @param outNet A vector containing which networks to return. Options include "regNet", "coregNet", "coopNet" and "regSign".
#' @param zScaling Boolean to indicate use of Z-Scores in output. FALSE will use [0,1] scale.
#' @param showProgress Boolean to indicate printing of output for algorithm progress.
#' @param randomizationMethod Method by which to randomize gene expression matrix. Default "None". Must be one of "None", "within.gene", "by.genes". "within.gene" randomization scrambles each row of the gene expression matrix, "by.gene" scrambles gene labels.
#' @param scaleByPresent Boolean to indicate scaling of correlations by percentage of positive samples.
#' @param filterExpr Boolean to remove genes with zero expression across all cells
#'   before network inference. Default FALSE.
#' @param samples Optional vector with the sample (donor) of each cell, used only for
#'   \code{regSign} to remove sample effects.
#' @param signParams Optional named list of further arguments passed to
#'   \code{\link{inferEdgeSign}} for \code{regSign} (e.g. \code{list(nBoot = 50)}). The
#'   direction step uses single cells (\code{gammaValue = 1}) by default, independently of
#'   the \code{gammaValue} used to build the network; set \code{gammaValue} here to change
#'   it. Prior edges are weighted by this network's \code{regNet} strength.
#' @return A list describing the inferred networks at convergence:
#'   \itemize{
#'     \item{regNet: Regulatory network matrix (TFs × genes)}
#'     \item{coregNet: Co-regulation network matrix (genes × genes)}
#'     \item{coopNet: Cooperation network matrix (TFs × TFs)}
#'     \item{regSign: Signed regulatory network (TFs × genes, same dimensions as
#'       \code{regNet}): \eqn{d \cdot \Phi^{-1}((1 + \Phi(z)) / 2)}, where \eqn{z} is the
#'       \code{regNet} edge score and \eqn{d} is the edge direction (+1 activating, -1
#'       repressing). For a prior edge, \eqn{d} is the direction called from this network's
#'       cells by \code{\link{inferEdgeSign}} when the data support one (including when it
#'       conflicts with the prior), and otherwise the prior's expected sign: its \code{mor}
#'       when signed, and for unsigned edges the prevailing sign of the prior's signed edges
#'       (+1 when the prior has no signs). TF activity is estimated from the TF's targets,
#'       never from its own mRNA, and everything is computed from this network's cells only.
#'       TF-gene pairs absent from the prior are NA. The magnitude grows with the edge score;
#'       it carries direction but not certainty, which is reported in \code{edgeSign}. With
#'       \code{zScaling = FALSE} it is \eqn{d \cdot \Phi(z)}, in [-1, 1].}
#'     \item{edgeSign: The per-edge table from \code{\link{inferEdgeSign}} (returned with
#'       \code{regSign}).}
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
#' # $ tf      :'data.frame':	356709 obs. of  4 variables:
#' # ..$ source_genesymbol: chr [1:356709] "MYC" "SPI1" "JUN" "JUND" ...
#' # ..$ target_genesymbol: chr [1:356709] "TERT" "BGLAP" "JUN" "JUN" ...
#' # ..$ weight           : num [1:356709] 1 1 1 1 1 1 1 1 1 1 ...
#' # ..$ mor              : int [1:356709] 1 1 1 1 1 1 1 1 1 1 ...
#' # ..- attr(*, "origin")= chr "cache"
#' # ..- attr(*, "url")= chr "https://omnipathdb.org/interactions? __truncated__
#' # $ ppi     :'data.frame':	4013 obs. of  3 variables:
#' # ..$ source_genesymbol: chr [1:4013] "ZIC1" "HES5" "DLL1" "SIX1" ...
#' # ..$ target_genesymbol: chr [1:4013] "ATOH1" "ATOH1" "NOTCH1" "OTX2" ...
#' # ..$ weight           : num [1:4013] 1 1 1 1 1 1 1 1 1 1 ...
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
#' # List of 8
#' # $ regNet  : num [1:966, 1:281] 0.0583 -0.1429 0.0367 -0.0133 1.868 ...
#' # ..- attr(*, "dimnames")=List of 2
#' # .. ..$ : chr [1:966] "AATF" "ABL1" "ACSS2" "ACTL6A" ...
#' # .. ..$ : chr [1:281] "ACKR1" "ACTA2" "ACTG2" "ADAMDEC1" ...
#' # $ coregNet: num [1:281, 1:281] 1.25e+06 0.672 -0.0285 1.25 0.614 ...
#' # ..- attr(*, "dimnames")=List of 2
#' # .. ..$ : chr [1:281] "ACKR1" "ACTA2" "ACTG2" "ADAMDEC1" ...
#' # .. ..$ : chr [1:281] "ACKR1" "ACTA2" "ACTG2" "ADAMDEC1" ...
#' # $ coopNet : num [1:966, 1:966] 9.25e+06 -4.57 6.24 6.65 -1.63 ...
#' # ..- attr(*, "dimnames")=List of 2
#' # .. ..$ : chr [1:966] "AATF" "ABL1" "ACSS2" "ACTL6A" ...
#' # .. ..$ : chr [1:966] "AATF" "ABL1" "ACSS2" "ACTL6A" ...
#' # $ regSign : num [1:966, 1:281] NA NA NA NA NA ...
#' # ..- attr(*, "dimnames")=List of 2
#' # .. ..$ : chr [1:966] "AATF" "ABL1" "ACSS2" "ACTL6A" ...
#' # .. ..$ : chr [1:281] "ACKR1" "ACTA2" "ACTG2" "ADAMDEC1" ...
#' # $ edgeSign:'data.frame':	5942 obs. of  15 variables: (see inferEdgeSign)
#' # $ numGenes: int 281
#' # $ numTFs  : int 966
#' # $ numEdges: int 271446
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
                     outNet = c("regNet", "coregNet", "coopNet", "regSign"),
                     zScaling = TRUE,
                     showProgress = TRUE,
                     randomizationMethod = "None",
                     scaleByPresent = FALSE,
                     filterExpr = FALSE,
                     samples = NULL,
                     signParams = list()) {
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
  if (!all(outNet %in% c("regNet", "coregNet", "coopNet", "regSign"))) {
    cli::cli_abort("outNet must be one or more of: regNet, coregNet, coopNet, regSign")
  }
  if (!is.null(samples) && length(samples) != ncol(gexMatrix)) {
    cli::cli_abort("samples must have one value per column of gexMatrix")
  }
  if (!is.list(signParams)) {
    cli::cli_abort("signParams must be a named list")
  }

  if (showProgress) {
    cli::cli_h1("SCORPION")
  }

  # Control BLAS threading to respect nCores
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    old_blas <- RhpcBLASctl::blas_get_num_procs()
    RhpcBLASctl::blas_set_num_threads(nCores)
    on.exit(RhpcBLASctl::blas_set_num_threads(old_blas), add = TRUE)
  }

  if (isTRUE(filterExpr)) {
    gexMatrix <- gexMatrix[rowSums(gexMatrix) > 0, ]
  }
  wantSign <- "regSign" %in% outNet && !is.null(tfMotifs)
  cellGex <- if (wantSign) gexMatrix else NULL
  gexMatrix <- makeSuperCells(X = gexMatrix, gamma = gammaValue, n.pc = nPC, fast.pca = FALSE)

  if (is.null(ppiNet) & is.null(tfMotifs)) {
    # Densify once; centering already destroys sparsity, so keep a single dense
    # base copy and use BLAS tcrossprod (dsyrk) for the co-expression matrix.
    gexMatrix <- as.matrix(gexMatrix)
    if (assocMethod == "spearman") {
      gexMatrix <- t(apply(gexMatrix, 1, rank))
    }
    geneCoExpr <- gexMatrix - rowMeans(gexMatrix)
    geneCoExpr <- geneCoExpr / sqrt(rowSums(geneCoExpr^2))
    return(tcrossprod(geneCoExpr))
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
    output = if (wantSign) union(outNet, "regNet") else outNet,
    zScale = zScaling,
    progress = showProgress,
    randomize = randomizationMethod,
    assoc.method = assocMethod,
    scale.by.present = scaleByPresent
  )

  # Direction of regulatory edges, computed after PANDA so regNet is unaffected
  if (wantSign) {
    if (showProgress) {
      cli::cli_alert_info("Inferring edge directions (regSign)")
    }
    edgeSign <- .scorpion_edgeSign(cellGex, tfMotifs, outNetworks$regNet, zScaling, samples, nPC, signParams)
    regSign <- .es_signNetwork(outNetworks$regNet, edgeSign, .es_prior(tfMotifs), zScaling)
    if (!"regNet" %in% outNet) {
      outNetworks$regNet <- NULL
    }
    counts <- outNetworks[c("numGenes", "numTFs", "numEdges")]
    outNetworks <- c(outNetworks[setdiff(names(outNetworks), names(counts))],
                     list(regSign = regSign, edgeSign = edgeSign), counts)
  }
  return(outNetworks)
}
