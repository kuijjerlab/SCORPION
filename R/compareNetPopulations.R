#' @export compareNetworkPopulations
#' @importFrom GSVA gsva
#' @importFrom stats t.test wilcox.test kruskal.test p.adjust
#' @author Daniel Osorio <dcosorioh@utexas.edu>
#' @title Evaluates gene differential regulation between population of networks based on mahalinobis distances.
#' @description Perform gene set scoring using a gene regulatory network as an input.
#' @param grnList A list of gene regulatory networks with transcription factors in the rows and target genes in the columns
#' @param groupID A vector with the group ID for each gene regulatory network
#' @param geneSets Gene sets provided as a \emph{list}
#' @param scoringType Either \emph{'gSet'} or \emph{'gene'}. If \emph{'gSet'} the single-sample Gene Set Enrichment Score (ssGSEA) for each gene set is returned. If \emph{'gene'} the mahalinobis distance is returned for each gene.
#' @param test The statistical test to be performed. By default uses the Wilcoxon Sum Rank Test, if more than two groups are about to be compared the Kruskal-Wallis Rank Sum Test should be used.
#' @return A \emph{data.frame} single-sample gene set enrichment score for each gene set or the mahalinobis distance for each target gene in each network. P-values and FDR values are also provided.

compareNetworkPopulations <- function(grnList, groupID, geneSets, scoringType = 'gSet', test = wilcox.test()){
  # Getting the list of all the genes present in all networks
  geneList <- unique(unlist(lapply(grnList, colnames)))
  # Getting the score for each network
  netScores <- sapply(grnList, function(X){
    mahalanobis(t(X), center = rowMeans(X), cov = cov(t(X)))[geneList]
  })
  # First option of return
  if (scoringType == 'gene'){
    gDiff <- apply(netScores, 1, function(g){test(g~groupID)$p.value})
    gDiff <- data.frame(netScores, P = gDiff)
    gDiff$FDR <- p.adjust(gDiff$P, method = 'fdr')
    return(gDiff)
  }
  # Second option of return
  if (scoringType == 'gene'){
    eScores <- suppressWarnings(gsva(netScores, geneSets, method = 'ssgsea'))
    gSetDiff <- apply(eScores, 1, function(g){test(g~groupID)$p.value})
    gSetDiff <- data.frame(eScores, P = gSetDiff)
    gSetDiff$FDR <- p.adjust(gSetDiff$P, method = 'fdr')
    return(gSetDiff)
  }
}
