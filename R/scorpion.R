#' @importFrom Matrix Matrix
#' @export scorpion
#' @title Constructs PANDA gene regulatory networks from single-cell gene expression data
#' @description  Constructs gene regulatory networks from single-cell gene expression data using the PANDA (Passing Attributes between Networks for Data Assimilation) algorithm.
#' @author Daniel Osorio <daniecos@uio.no>
#' @param tfMotifs A motif dataset, a data.frame or a matrix containing 3 columns. Each row describes an motif associated with a transcription factor (column 1) a gene (column 2) and a score (column 3) for the motif.
#' @param gexMatrix An expression dataset, with genes in the rows and barcodes (cells) in the columns.
#' @param ppiNet A Protein-Protein-Interaction dataset, a data.frame or matrix containing 3 columns. Each row describes a protein-protein interaction between transcription factor 1(column 1), transcription factor 2 (column 2) and a score (column 3) for the interaction.
#' @param nCores Number of processors to be used if BLAS or MPI is active.
#' @param gammaValue Graining level of data (proportion of number of single cells in the initial dataset to the number of super-cells in the final dataset)
#' @param nPC Number of principal components to use for construction of single-cell kNN network.
#' @param alphaValue Value to be used for update variable.
#' @param hammingValue Value at which to terminate the process based on Hamming distance.
#' @param assocMethod Association method. Must be one of 'pearson', 'spearman' or 'pcNet'.
#' @param nIter Sets the maximum number of iterations PANDA can run before exiting.
#' @param outNet A vector containing which networks to return. Options include "regulatory", "coregulatory", "cooperative".
#' @param zScaling Boolean to indicate use of Z-Scores in output. False will use [0,1] scale.
#' @param showProgress Boolean to indicate printing of output for algorithm progress.
#' @param randomizationMethod Method by which to randomize gene expression matrix. Default "None". Must be one of "None", "within.gene", "by.genes". "within.gene" randomization scrambles each row of the gene expression matrix, "by.gene" scrambles gene labels.
#' @param scaleByPresent Boolean to indicate scaling of correlations by percentage of positive samples
#' @return A list of matrices describing networks achieved by convergence with PANDA algorithm.

scorpion <- function(tfMotifs = NULL,
                     gexMatrix,
                     ppiNet = NULL,
                     nCores = 1,
                     gammaValue = 10,
                     nPC = 25,
                     assocMethod = 'pearson',
                     alphaValue = 0.1,
                     hammingValue = 0.001,
                     nIter = Inf,
                     outNet = c('regulatory', 'coregulatory', 'cooperative'),
                     zScaling = TRUE,
                     showProgress = TRUE,
                     randomizationMethod = 'None',
                     scaleByPresent = FALSE) {
  gexMatrix <- gexMatrix[rowSums(gexMatrix) > 0,]
  gexMatrix <- makeSuperCells(X = gexMatrix, gamma = gammaValue, n.pc = nPC, fast.pca = TRUE)

  if(is.null(ppiNet) & is.null(tfMotifs)){
    if(assocMethod == 'spearman'){
      gexMatrix <- Matrix(t(apply(gexMatrix, 1, rank)))
    }
    geneCoExpr <- gexMatrix - rowMeans(gexMatrix)
    geneCoExpr <- geneCoExpr/sqrt(rowSums(geneCoExpr^2))
    geneCoExpr <- tcrossprod(geneCoExpr)
    return(geneCoExpr)
  }
  outNetworks <- runPANDA(motif = tfMotifs,
                          ppi = ppiNet,
                          expr = gexMatrix,
                          n.cores = nCores,
                          alpha = alphaValue,
                          hamming = hammingValue,
                          iter = nIter,
                          output = outNet,
                          zScale = zScaling,
                          progress = showProgress ,
                          randomize = randomizationMethod ,
                          assoc.method = assocMethod,
                          scale.by.present = scaleByPresent)
  return(outNetworks)
}
