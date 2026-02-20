update.diagonal <- function(diagMat, num, alpha, step) {
  isGPU <- grepl('gpu', class(diagMat))
  diagMat <- as.matrix(diagMat)
  diag(diagMat) <- NaN
  diagstd <- rowSds(diagMat, na.rm = TRUE) * sqrt((num - 2) / (num - 1))
  diag(diagMat) <- diagstd * num * exp(2 * alpha * step)
  if(isGPU){
    gpuMatrix <- getFromNamespace("gpuMatrix", "gpuR")
    diagMat <- gpuMatrix(diagMat)
  }
  return(diagMat)
}
