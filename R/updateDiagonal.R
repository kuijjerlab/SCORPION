update.diagonal <- function(diagMat, num, alpha, step) {
  diag(diagMat) <- NaN
  # Use column standard deviation to match Python implementation
  diagstd <- colSds(diagMat, na.rm = TRUE)
  diag(diagMat) <- diagstd * num * exp(2 * alpha * step)
  return(diagMat)
}
