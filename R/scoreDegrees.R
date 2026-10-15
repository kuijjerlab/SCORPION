#' @importFrom Matrix sparseMatrix
#' @export scoreDegrees
#' @title Per-cell out-degrees and in-degrees from a signed network (for visualization)
#' @description
#' Scores every cell with the signed, weighted network (e.g. \code{regSign} from
#' \code{\link{scorpion}} or \code{\link{runSCORPION}}) without building a network per cell.
#' With \eqn{w_{ij}} the signed weight of edge TF \eqn{i} to gene \eqn{j} and \eqn{z_{jc}} the
#' standardized expression of gene \eqn{j} in cell \eqn{c}:
#' \itemize{
#'   \item{Out-degree of TF \eqn{i} (its activity in the cell):
#'     \deqn{out_{ic} = \frac{\sum_j w_{ij} z_{jc}}{\sum_j |w_{ij}|}}
#'     Activated targets that are high raise it and repressed targets that are high lower it.
#'     The TF's own expression is never used.}
#'   \item{In-degree of gene \eqn{j} (the net regulatory input it receives):
#'     \deqn{in_{jc} = \frac{\sum_i w_{ij} \, out_{ic}^{(-j)}}{\sum_i |w_{ij}|}}
#'     positive for net activating input and negative for net repressive input, where
#'     \eqn{out_{ic}^{(-j)}} is the regulator's out-degree computed without gene \eqn{j}, so a
#'     gene never contributes to its own input.}
#' }
#' Each TF keeps its \code{topTargets} strongest edges (largest \eqn{|w_{ij}|}): with a
#' network that signs every edge, the many weak edges otherwise dilute the scores.
#'
#' The scores are intended for visualization (e.g. on an embedding or in a heatmap). The
#' edge directions were estimated from the same cells, so they should not be used for
#' formal tests between groups; use \code{\link{testEdges}} on the networks instead.
#' Score a group's cells with that group's own network.
#' @author Daniel Osorio <daniecos@uio.no>
#' @param gexMatrix An expression dataset with genes in the rows and cells in the columns.
#' @param network Signed network: either a TF x gene matrix (e.g. \code{scorpion()$regSign})
#'   or a data.frame with columns \code{tf}, \code{target} and the weights in the column
#'   named by \code{weightColumn} (e.g. the \code{"tf-target-sign"} rows of
#'   \code{runSCORPION()} output).
#' @param weightColumn For a data.frame \code{network}, the column holding the weights
#'   (e.g. a network identifier such as \code{"T"}). Ignored for a matrix.
#' @param topTargets Number of strongest edges (by absolute weight) kept per TF, or NULL to
#'   keep all. Default 50.
#' @param normalizeData Boolean; log-normalize \code{gexMatrix} (counts) before scoring.
#'   Default TRUE.
#' @param scaleGenes Boolean; standardize each gene across cells before scoring. Default TRUE.
#' @param minTargets Minimum number of targets with a non-zero weight in \code{gexMatrix}
#'   required to score a TF. Default 5.
#' @return A list with \code{outDegree} (TF x cell) and \code{inDegree} (gene x cell, for
#'   genes with at least one regulator).
#' @seealso \code{\link{scorpion}}, \code{\link{runSCORPION}}, \code{\link{inferEdgeSign}}
#' @examples
#' \donttest{
#' data(scorpionTest)
#' cellsN <- scorpionTest$metadata$region == "N"
#' net <- scorpion(
#'   tfMotifs = scorpionTest$tf,
#'   gexMatrix = scorpionTest$gex[, cellsN],
#'   ppiNet = scorpionTest$ppi,
#'   alphaValue = 0.8,
#'   outNet = c("regNet", "regSign"),
#'   signParams = list(nBoot = 20)
#' )
#' degrees <- scoreDegrees(scorpionTest$gex[, cellsN], net$regSign)
#' dim(degrees$outDegree)
#' dim(degrees$inDegree)
#' }
scoreDegrees <- function(gexMatrix,
                         network,
                         weightColumn = NULL,
                         topTargets = 50,
                         normalizeData = TRUE,
                         scaleGenes = TRUE,
                         minTargets = 5) {
  if (is.null(rownames(gexMatrix)) || is.null(colnames(gexMatrix))) {
    cli::cli_abort("gexMatrix must have gene names as rownames and cell barcodes as colnames")
  }
  if (!is.null(topTargets) && (!is.numeric(topTargets) || length(topTargets) != 1 || topTargets < 1)) {
    cli::cli_abort("topTargets must be a positive number or NULL")
  }
  # Signed weights as a TF x gene matrix
  if (is.data.frame(network)) {
    if (!all(c("tf", "target") %in% colnames(network))) {
      cli::cli_abort("A data.frame network must have columns 'tf' and 'target'")
    }
    if (is.null(weightColumn) || !weightColumn %in% colnames(network)) {
      cli::cli_abort("weightColumn must name the column of network holding the weights")
    }
    w <- as.numeric(network[[weightColumn]])
    keep <- is.finite(w) & w != 0
    tfs <- sort(unique(network$tf[keep]))
    targets <- sort(unique(network$target[keep]))
    W <- matrix(0, length(tfs), length(targets), dimnames = list(tfs, targets))
    W[cbind(match(network$tf[keep], tfs), match(network$target[keep], targets))] <- w[keep]
  } else if (is.matrix(network) || methods::is(network, "Matrix")) {
    if (is.null(rownames(network)) || is.null(colnames(network))) {
      cli::cli_abort("A matrix network must have TF rownames and gene colnames")
    }
    W <- as.matrix(network)
    W[!is.finite(W)] <- 0
  } else {
    cli::cli_abort("network must be a TF x gene matrix or a data.frame")
  }
  genes <- rownames(gexMatrix)[rownames(gexMatrix) %in% colnames(W)]   # gexMatrix order
  if (length(genes) == 0) {
    cli::cli_abort("No network targets are present in gexMatrix")
  }
  W <- W[, genes, drop = FALSE]
  # Keep each TF's strongest edges
  if (!is.null(topTargets)) {
    k <- min(as.integer(topTargets), ncol(W))
    for (i in seq_len(nrow(W))) {
      W[i, order(-abs(W[i, ]))[-seq_len(k)]] <- 0
    }
  }
  W <- W[rowSums(W != 0) >= minTargets, , drop = FALSE]
  if (nrow(W) == 0) {
    cli::cli_abort("No TF has at least {minTargets} targets with a non-zero weight in gexMatrix")
  }
  W <- W[, colSums(W != 0) > 0, drop = FALSE]

  X <- gexMatrix[colnames(W), , drop = FALSE]
  if (isTRUE(normalizeData)) {
    X <- log_normalize_data(gexMatrix)[colnames(W), , drop = FALSE]
  }
  Z <- as.matrix(X)
  if (isTRUE(scaleGenes)) {
    Z <- Z - rowMeans(Z)
    s <- sqrt(rowSums(Z^2) / max(ncol(Z) - 1, 1))
    Z <- Z / ifelse(s > 0, s, 1)
  }

  D <- rowSums(abs(W))                                     # per TF
  S <- W %*% Z                                             # TF x cell
  outDegree <- S / D
  # In-degree with leave-one-out regulator out-degrees:
  #   out_i^(-j) = (S_i - w_ij z_j) / (D_i - |w_ij|)
  inv <- ifelse(W != 0, 1 / pmax(D - abs(W), 1e-12), 0)    # TF x gene
  num <- t(W * inv) %*% S - Z * colSums(W^2 * inv)          # gene x cell
  inDegree <- num / colSums(abs(W))
  dimnames(outDegree) <- list(rownames(W), colnames(gexMatrix))
  dimnames(inDegree) <- list(colnames(W), colnames(gexMatrix))
  list(outDegree = outDegree, inDegree = inDegree)
}
