prepResult <- function(zScale, output, regulatoryNetwork, geneCoreg, tfCoopNetwork, edgelist, motif) {
  # if (!zScale){
  #   regulatoryNetwork <- stats::pnorm(regulatoryNetwork)
  #   geneCoreg         <- stats::pnorm(geneCoreg)
  #   tfCoopNetwork     <- stats::pnorm(tfCoopNetwork)
  # }
  # if("regulatory"%in%output){
  #   # if(edgelist){
  #   #   regulatoryNetwork <- melt.array(regulatoryNetwork)
  #   #   colnames(regulatoryNetwork) <- c("TF", "Gene", "Score")
  #   #   regulatoryNetwork$Motif <- as.numeric(with(regulatoryNetwork, paste0(TF, Gene))%in%paste0(motif[,1],motif[,2]))
  #   # }
  #   resList$regNet <- regulatoryNetwork
  # }
  # if("coexpression"%in%output){
  #   # if(edgelist){
  #   #   geneCoreg <- melt.array(geneCoreg)
  #   #   colnames(geneCoreg) <- c("Gene.x", "Gene.y", "Score")
  #   # }
  #   resList$coregNet <- geneCoreg
  # }
  # if("cooperative"%in%output){
  #   # if(edgelist){
  #   #   tfCoopNetwork <- melt.array(tfCoopNetwork)
  #   #   colnames(tfCoopNetwork) <- c("TF.x", "TF.y", "Score")
  #   # }
  #   resList$coopNet <- tfCoopNetwork
  # }
  resList <- list()
  numGenes <- dim(geneCoreg)[1]
  numTFs <- dim(tfCoopNetwork)[1]
  numEdges <- sum(apply(regulatoryNetwork, 1, function(X) {
    X != 0
  }))
  list(regNet = regulatoryNetwork, coregNet = geneCoreg, coopNet = tfCoopNetwork, numGenes = numGenes, numTFs = numTFs, numEdges = numEdges)
}
