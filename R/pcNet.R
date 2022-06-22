#' @importFrom pbapply pbsapply
pcNet <- function(X,
                  nComp = 3,
                  scaleScores = TRUE,
                  symmetric = FALSE,
                  q = 0, verbose = TRUE,
                  nCores = 1) {
  if (!all(Matrix::rowSums(X) > 0)) {
    stop('Quality control has not been applied over the matrix.')
  }
  xClass <- class(X)[[1]]
  validClass <- xClass %in% c('matrix', 'dgCMatrix')
  if (!validClass) {
    stop('Input should be a matrix with cells as columns and genes as rows')
  }
  if (nComp < 2 | nComp >= nrow(X)) {
    stop('nCom should be greater or equal than 2 and lower than the total number of genes')
  }
  gNames <- rownames(X)
  pcCoefficients <- function(K) {
    # Taking out the gene to be regressed out
    y <- X[, K]
    Xi <- X
    Xi <- Xi[, -K]
    # Step 1: Perform PCA on the observed covariates data matrix to obtain $n$ number of the principal components.
    coeff <- irlba::irlba(Xi, nComp)$v
    score <- Xi %*% coeff
    # Step 2: Regress the observed vector of outcomes on the selected principal components as covariates, using ordinary least squares regression to get a vector of estimated regression coefficients.
    score <-
      Matrix::t(Matrix::t(score) / (apply(score, 2, function(X) {
        sqrt(sum(X ^ 2))
      }) ^ 2))
    # Step 3: Transform this vector back to the scale of the actual covariates, using the eigenvectors corresponding to the selected principal components to get the final PCR estimator for estimating the regression coefficients characterizing the original model.
    Beta <- colSums(y * score)
    Beta <- coeff %*% (Beta)

    return(Beta)
  }

  # Standardizing the data
  X <- (scale(Matrix::t(X)))

  # Identify the number of rows in the input matrix
  n <- ncol(X)

  # Generate the output matrix
  A <- 1 - diag(n)

  # Apply the principal component regression for each gene
  # RhpcBLASctl::omp_set_num_threads(nCores)
  # RhpcBLASctl::blas_set_num_threads(nCores)

  if(verbose){
    B <- pbapply::pbsapply(seq_len(n), pcCoefficients)
  } else {
    B <- sapply(seq_len(n), pcCoefficients)
  }

  # Transposition of the Beta coefficient matrix
  B <- t(B)

  # Replacing the values in the output matrix
  for (K in seq_len(n)) {
    A[K, A[K, ] == 1] = B[K, ]
  }

  # Making the output matrix symmetric
  if (isTRUE(symmetric)) {
    A <- (A + t(A)) / 2
  }

  # Absolute values for scaling and filtering
  absA <- abs(A)

  # Scaling the output matrix
  if (isTRUE(scaleScores)) {
    A <- (A / max(absA))
  }

  # Filtering the output matrix
  A[absA < stats::quantile(absA, q)] <- 0

  # Setting the diagonal to be 0
  diag(A) <- 0

  # Adding names
  colnames(A) <- rownames(A) <- gNames

  # Making the output a sparse matrix
  A <- Matrix(A)

  # Return
  return(A)
}
