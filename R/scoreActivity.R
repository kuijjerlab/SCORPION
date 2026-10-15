#' @importFrom Matrix sparseMatrix
#' @export scoreActivity
#' @title Score TF activity in each cell from a signed network (for visualization)
#' @description
#' Computes a per-cell activity score for each transcription factor (TF) as the
#' weighted mean of its targets' standardized expression, using a signed network such as
#' \code{regSign} from \code{\link{scorpion}} or \code{\link{runSCORPION}}:
#' \deqn{activity_{ic} = \frac{\sum_j w_{ij} z_{jc}}{\sum_j |w_{ij}|}}
#' where \eqn{w_{ij}} is the signed edge weight and \eqn{z_{jc}} the standardized
#' expression of target \eqn{j} in cell \eqn{c}. Activated targets that are high in a cell
#' raise the score and repressed targets that are high lower it. Edges with NA or zero
#' weight are ignored, and the TF's own expression is never used.
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
#' @param normalizeData Boolean; log-normalize \code{gexMatrix} (counts) before scoring.
#'   Default TRUE.
#' @param scaleGenes Boolean; standardize each gene across cells before scoring. Default TRUE.
#' @param minTargets Minimum number of targets with a non-zero weight in \code{gexMatrix}
#'   required to score a TF. Default 5.
#' @return A TF x cell matrix of activity scores. The number of targets used per TF is in
#'   the \code{"nTargets"} attribute.
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
#' activity <- scoreActivity(scorpionTest$gex[, cellsN], net$regSign)
#' dim(activity)
#' }
scoreActivity <- function(gexMatrix,
                          network,
                          weightColumn = NULL,
                          normalizeData = TRUE,
                          scaleGenes = TRUE,
                          minTargets = 5) {
  if (is.null(rownames(gexMatrix)) || is.null(colnames(gexMatrix))) {
    cli::cli_abort("gexMatrix must have gene names as rownames and cell barcodes as colnames")
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
    W <- Matrix::sparseMatrix(i = match(network$tf[keep], tfs), j = match(network$target[keep], targets),
                              x = w[keep], dims = c(length(tfs), length(targets)), dimnames = list(tfs, targets))
  } else if (is.matrix(network) || methods::is(network, "Matrix")) {
    if (is.null(rownames(network)) || is.null(colnames(network))) {
      cli::cli_abort("A matrix network must have TF rownames and gene colnames")
    }
    W <- as.matrix(network)
    W[!is.finite(W)] <- 0
    W <- Matrix::Matrix(W, sparse = TRUE)
  } else {
    cli::cli_abort("network must be a TF x gene matrix or a data.frame")
  }
  genes <- intersect(colnames(W), rownames(gexMatrix))
  if (length(genes) == 0) {
    cli::cli_abort("No network targets are present in gexMatrix")
  }
  W <- W[, genes, drop = FALSE]
  nTargets <- Matrix::rowSums(W != 0)
  W <- W[nTargets >= minTargets, , drop = FALSE]
  if (nrow(W) == 0) {
    cli::cli_abort("No TF has at least {minTargets} targets with a non-zero weight in gexMatrix")
  }

  X <- gexMatrix[genes, , drop = FALSE]
  if (isTRUE(normalizeData)) {
    X <- log_normalize_data(X)
  }
  X <- as.matrix(X)
  if (isTRUE(scaleGenes)) {
    X <- X - rowMeans(X)
    s <- sqrt(rowSums(X^2) / max(ncol(X) - 1, 1))
    X <- X / ifelse(s > 0, s, 1)
  }
  activity <- as.matrix(W %*% X) / Matrix::rowSums(abs(W))
  dimnames(activity) <- list(rownames(W), colnames(gexMatrix))
  attr(activity, "nTargets") <- nTargets[rownames(W)]
  activity
}
