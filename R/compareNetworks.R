#' @importFrom stats mahalanobis cov
#' @export compareNetworks
#' @title Calculate differences between two gene regulatory networks
#' @description Using the mahalinobis distance, calculates the gene in-degree for two different PANDA gene regulatory networks.
#' @author Daniel Osorio <daniecos@uio.no>
#' @param X GRN1
#' @param Y GRN2

compareNetworks <- function(X, Y){
  tfList <- intersect(rownames(X), rownames(Y))
  geneList <- intersect(colnames(X), colnames(Y))
  nGenes <- length(geneList)
  X <- t(X[tfList, geneList])
  Y <- t(Y[tfList, geneList])
  U <- rbind(X, Y)
  D <- stats::mahalanobis(U, colMeans(U), stats::cov(U))
  D1 <- D[1:nGenes]
  D2 <- D[(nGenes + 1):(2 * nGenes)]
  D <- D1 - D2
  names(D) <- geneList
  return(D)
}
