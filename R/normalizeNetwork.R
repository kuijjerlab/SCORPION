normalizeNetwork <- function(X) {
  nr <- nrow(X)
  nc <- ncol(X)

  # overall values
  mu0 <- mean(X)
  std0 <- sd(X)

  # Row normalization (R's column-major recycling handles row-wise broadcast)
  mu1 <- rowMeans(X)
  std1 <- rowSds(X) * sqrt((nc - 1) / nc)
  Z1 <- (X - mu1) / std1

  # Column normalization (needs explicit broadcast)
  mu2 <- colMeans(X)
  std2 <- colSds(X) * sqrt((nr - 1) / nr)
  mu2 <- rep(mu2, each = nr)
  dim(mu2) <- c(nr, nc)
  std2 <- rep(std2, each = nr)
  dim(std2) <- c(nr, nc)
  Z2 <- (X - mu2) / std2

  # combine and return
  normMat <- (Z1 + Z2) / sqrt(2)

  Z0 <- (X - mu0) / std0
  f1 <- is.na(Z1)
  f2 <- is.na(Z2)

  normMat[f1] <- (Z2[f1] + Z0[f1]) / sqrt(2)
  normMat[f2] <- (Z1[f2] + Z0[f2]) / sqrt(2)
  normMat[f1 & f2] <- 2 * Z0[f1 & f2] / sqrt(2)

  normMat
}
