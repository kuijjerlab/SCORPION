#' @export scorpion
#' @title Constructs PANDA gene regulatory networks from single-cell gene expression data
#' @description  Constructs gene regulatory networks from single-cell gene expression data using the PANDA (Passing Attributes between Networks for Data Assimilation) algorithm
#' @author Daniel Osorio <daniecos@uio.no>
#' @param motif TF Motif
#' @param expr Single-cell Expression Data
#' @param ppi PPI
#' @return A list
scorpion <- function(motif = NULL,
                     expr,
                     ppi = NULL,
                     nCores = 1) {

  expr <- makeSuperCells(expr)
  expr <- as.matrix(expr)
  O <- runPANDA(motif = motif, ppi = ppi, expr = expr, n.cores = nCores)
  return(O)
}
