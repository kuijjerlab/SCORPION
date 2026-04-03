tanimoto <- function(X, Y = NULL,
                     x_norm_sq = NULL, y_norm_sq = NULL,
                     type = c("general", "tcrossprod", "crossprod")) {
  type <- match.arg(type)

  if (type == "tcrossprod") {
    Amat <- tcrossprod(X)
    if (is.null(x_norm_sq)) x_norm_sq <- rowSums(X * X)
    y_norm_sq <- x_norm_sq

  } else if (type == "crossprod") {
    Amat <- crossprod(X)
    if (is.null(x_norm_sq)) x_norm_sq <- colSums(X * X)
    y_norm_sq <- x_norm_sq

  } else {
    Amat <- X %*% Y
    if (is.null(y_norm_sq)) y_norm_sq <- colSums(Y * Y)
    if (is.null(x_norm_sq)) x_norm_sq <- rowSums(X * X)
  }

  # Column-by-column normalization to avoid full-size outer() temporary
  nc <- ncol(Amat)
  for (j in seq_len(nc)) {
    col_j <- Amat[, j]
    Amat[, j] <- col_j / sqrt(x_norm_sq + y_norm_sq[j] - abs(col_j))
  }
  Amat
}

#' @importFrom utils getFromNamespace
tanimoto_gpu <- function(X, Y = NULL,
                         x_norm_sq = NULL, y_norm_sq = NULL,
                         type = c("general", "tcrossprod", "crossprod"),
                         gpuMatrixFn = NULL) {
  type <- match.arg(type)
  if (is.null(gpuMatrixFn)) gpuMatrixFn <- getFromNamespace("gpuMatrix", "gpuR")

  if (type == "tcrossprod") {
    gX <- gpuMatrixFn(X)
    Amat <- as.matrix(gX %*% gpuMatrixFn(t(X)))
    if (is.null(x_norm_sq)) x_norm_sq <- rowSums(X * X)
    y_norm_sq <- x_norm_sq

  } else if (type == "crossprod") {
    gX <- gpuMatrixFn(X)
    Amat <- as.matrix(gpuMatrixFn(t(X)) %*% gX)
    if (is.null(x_norm_sq)) x_norm_sq <- colSums(X * X)
    y_norm_sq <- x_norm_sq

  } else {
    Amat <- as.matrix(gpuMatrixFn(X) %*% gpuMatrixFn(Y))
    if (is.null(y_norm_sq)) y_norm_sq <- colSums(Y * Y)
    if (is.null(x_norm_sq)) x_norm_sq <- rowSums(X * X)
  }

  nc <- ncol(Amat)
  for (j in seq_len(nc)) {
    col_j <- Amat[, j]
    Amat[, j] <- col_j / sqrt(x_norm_sq + y_norm_sq[j] - abs(col_j))
  }
  Amat
}
