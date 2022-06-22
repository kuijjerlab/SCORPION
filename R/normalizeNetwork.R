normalizeNetwork <- function(X) {

  nr = nrow(X)
  nc = ncol(X)
  dm = c(nr, nc)

  # overall values
  mu0 = mean(X)
  std0 = sd(X) * sqrt((nr * nc - 1) / (nr * nc))

  # operations on rows
  mu1 = rowMeans(X) # operations on rows
  std1 = rowSds(X) * sqrt((nc - 1) / nc)

  mu1 = rep(mu1, nc)
  dim(mu1) = dm
  std1 = rep(std1, nc)
  dim(std1) = dm

  Z1 = (X - mu1) / std1

  # operations on columns
  mu2 = colMeans(X) # operations on columns
  std2 = colSds(X) * sqrt((nr - 1) / nr)

  mu2 = rep(mu2, each = nr)
  dim(mu2) = dm
  std2 = rep(std2, each = nr)
  dim(std2) = dm

  Z2 = (X - mu2) / std2

  # combine and return
  normMat = Z1 / sqrt(2) + Z2 / sqrt(2)

  return(normMat)
}
