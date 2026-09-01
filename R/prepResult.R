prepResult <- function(zScale, output, regulatoryNetwork, geneCoreg, tfCoopNetwork, edgelist, motif) {
  # Coerce to base R dense matrices so the returned classes are consistent
  # regardless of assocMethod (e.g. pcNet yields an S4 Matrix, pearson a base matrix)
  regulatoryNetwork <- as.matrix(regulatoryNetwork)
  geneCoreg <- as.matrix(geneCoreg)
  tfCoopNetwork <- as.matrix(tfCoopNetwork)
  if (!zScale) {
    regulatoryNetwork <- stats::pnorm(regulatoryNetwork)
    geneCoreg <- stats::pnorm(geneCoreg)
    tfCoopNetwork <- stats::pnorm(tfCoopNetwork)
  }
  resList <- list()
  numGenes <- dim(geneCoreg)[1]
  numTFs <- dim(tfCoopNetwork)[1]
  numEdges <- sum(apply(regulatoryNetwork, 1, function(X) {
    X != 0
  }))
  if ('regNet' %in% output) {
    resList$regNet <- regulatoryNetwork
  }
  if ('coregNet' %in% output) {
    resList$coregNet <- geneCoreg
  }
  if ('coopNet' %in% output) {
    resList$coopNet <- tfCoopNetwork
  }
  resList$numGenes <- numGenes
  resList$numTFs <- numTFs
  resList$numEdges <- numEdges

  return(resList)
}
