log_normalize_data <- function(X){
  X <- log1p(t(t(X) / colSums(X))* 10000)
  return(X)
}