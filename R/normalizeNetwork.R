normalizeNetwork <- function(X) {
  nr <- nrow(X)
  nc <- ncol(X)

  # overall values
  mu0 <- mean(X)
  std0 <- sd(X)

  # Row stats
  mu1 <- rowMeans(X)
  std1 <- rowSds(X) * sqrt((nc - 1) / nc)

  # Column stats (vectors only — no full-matrix broadcast)
  mu2 <- colMeans(X)
  std2 <- colSds(X) * sqrt((nr - 1) / nr)

  has_zero_var <- any(std1 == 0) || any(std2 == 0)

  if (!has_zero_var) {
    # Fast path: accumulate Z1 + Z2 in a single matrix
    normMat <- (X - mu1) / std1
    for (j in seq_len(nc)) {
      normMat[, j] <- normMat[, j] + (X[, j] - mu2[j]) / std2[j]
    }
    normMat <- normMat / sqrt(2)
  } else {
    # Slow path: handle NaN from zero-variance rows/columns
    Z1 <- (X - mu1) / std1

    Z2 <- X
    for (j in seq_len(nc)) {
      Z2[, j] <- (X[, j] - mu2[j]) / std2[j]
    }

    normMat <- (Z1 + Z2) / sqrt(2)

    Z0 <- (X - mu0) / std0
    f1 <- is.na(Z1)
    f2 <- is.na(Z2)

    normMat[f1] <- (Z2[f1] + Z0[f1]) / sqrt(2)
    normMat[f2] <- (Z1[f2] + Z0[f2]) / sqrt(2)
    normMat[f1 & f2] <- 2 * Z0[f1 & f2] / sqrt(2)
  }

  normMat
}
