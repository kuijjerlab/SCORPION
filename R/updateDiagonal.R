update.diagonal <- function(diagMat, num, alpha, step) {
  seqs = seq(1, num * num, num + 1)
  diagMat[seqs] = NaN
  diagstd = rowSds(diagMat, na.rm = TRUE) * sqrt((num - 2) / (num - 1))
  diagMat[seqs] = diagstd * num * exp(2 * alpha * step)
  return(diagMat)
}
