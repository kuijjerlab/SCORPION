update.diagonal <- function(diagMat, num, alpha, step) {
  isGPU <- grepl('gpu', class(diagMat))
  diagMat <- as.matrix(diagMat)
  diag(diagMat) <- NaN
  # Use column standard deviation to match Python implementation
  diagstd <- colSds(diagMat, na.rm = TRUE)
  diag(diagMat) <- diagstd * num * exp(2 * alpha * step)
  if(isGPU){
    gpuMatrix <- getFromNamespace("gpuMatrix", "gpuR")
    diagMat <- gpuMatrix(diagMat)
  }
  return(diagMat)
}
