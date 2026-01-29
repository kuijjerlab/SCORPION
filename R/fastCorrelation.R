fastCorrelation <- function(X, Y, method = 'pearson'){
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  
  stopifnot(nrow(X) == nrow(Y))
  
  if(method == 'spearman'){
    # Rank columns
    RX <- apply(X, 2, rank)
    RY <- apply(Y, 2, rank)  
  } else {
    RX <- X
    RY <- Y
  }
  
  # Center
  RX <- sweep(RX, 2, colMeans(RX), "-")
  RY <- sweep(RY, 2, colMeans(RY), "-")
  
  # Norms
  nx <- sqrt(colSums(RX * RX))
  ny <- sqrt(colSums(RY * RY))
  
  # Correlation
  t(RX) %*% RY / (nx %o% ny)
}
