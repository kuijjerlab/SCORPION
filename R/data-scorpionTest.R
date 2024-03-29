#' @name scorpionTest
#' @docType data
#' @title Example single-cell gene expression, motif, and ppi data
#' @description This data is a list containing three objects. The motif \code{data.frame} describes a set of pairwise connections where a specific known sequence motif of a transcription factor was found upstream of the corresponding gene. The expression \code{dgCMatrix} is a set of 230 gene expression levels measured across 80 cells. Finally, the ppi \code{data.frame} describes a set of known pairwise protein-protein interactions.
#' @usage data(scorpionTest)
#' @format A list containing three datasets.
#' \describe{
#' \item{\code{gex}}{A subsetted version of 10X Genomics' 3k PBMC dataset provided by the \code{Seurat} package.}
#' \item{\code{tf}}{Subset of the transcription-factor and target gene list provided by the \code{dorothea} package for Homo sapiens.}
#' \item{\code{ppi}}{The known protein-protein interactions and the combined score downloaded from the STRING database}
#' }
#' @examples
#' # Loading example data
#' data(scorpionTest)
#'
#' # The structure of the data
#' str(scorpionTest)
#'
#' # List of 3
#' # $ gex:Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#' # .. ..@ i       : int [1:4456] 1 5 8 11 22 30 33 34 36 38 ...
#' # .. ..@ p       : int [1:81] 0 47 99 149 205 258 306 342 387 423 ...
#' # .. ..@ Dim     : int [1:2] 230 80
#' # .. ..@ Dimnames:List of 2
#' # .. .. ..$ : chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
#' # .. .. ..$ : chr [1:80] "ATGCCAGAACGACT" "CATGGCCTGTGCAT" "GAACCTGATGAACC" "TGACTGGATTCTCA" ...
#' # .. ..@ x       : num [1:4456] 1 1 3 1 1 4 1 5 1 1 ...
#' # .. ..@ factors : list()
#' # $ tf :'data.frame':	4485 obs. of  3 variables:
#' #   ..$ tf    : chr [1:4485] "ADNP" "ADNP" "ADNP" "AEBP2" ...
#' # ..$ target: chr [1:4485] "PRF1" "TMEM40" "TNFRSF1B" "CFP" ...
#' # ..$ mor   : num [1:4485] 1 1 1 1 1 1 1 1 1 1 ...
#' # $ ppi:'data.frame':	12754 obs. of  3 variables:
#' #   ..$ X.node1       : chr [1:12754] "ADNP" "ADNP" "ADNP" "AEBP2" ...
#' # ..$ node2         : chr [1:12754] "ZBTB14" "NFIA" "CDC5L" "YY1" ...
#' # ..$ combined_score: num [1:12754] 0.769 0.64 0.581 0.597 0.54 0.753 0.659 0.548 0.59 0.654 ...

NULL
