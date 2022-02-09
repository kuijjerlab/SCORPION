#' @export netScoring
#' @importFrom GSVA gsva
#' @author Daniel Osorio <daniecos@uio.no>
#' @title Gene Set Network Scoring
#' @description Perform gene set scoring using a gene regulatory network as an input.
#' @param regNetwork A gene regulatory network with transcription factors in the rows and target genes in the columns
#' @param geneSets Gene sets provided as a \emph{list}
#' @param scoreType Either \emph{'gSet'} or \emph{'gene'}. If \emph{'gSet'} the single-sample Gene Set Enrichment Score (ssGSEA) for each gene set is returned. If \emph{'gene'} the mahalinobis distance is returned for each gene.
#' @return A named vector with the single-sample gene set enrichment score for each gene set or the mahalinobis distance for each target gene in the network.

netScoring <- function(regNetwork, geneSets, scoreType = 'gSet') {
  MD <-
    as.matrix(stats::mahalanobis(
      t(regNetwork),
      center = Matrix::rowMeans(regNetwork),
      cov = stats::cov(t(regNetwork))
    ))
  colnames(MD) <- 'MD'
  O <- MD[, 1]
  if (scoreType == 'gSet') {
    O <-
      base::suppressWarnings(GSVA::gsva(MD, geneSets, method = 'ssgsea', verbose = FALSE)[, 1])
  }
  return(O)
}
