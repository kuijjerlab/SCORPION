colSds <-
  function(X, na.rm = FALSE) {
    n <- colSums(!is.na(X))
    col_means <- colMeans(X, na.rm = na.rm)
    col_mean_sq <- colMeans(X^2, na.rm = na.rm)
    sqrt(n / (n - 1) * (col_mean_sq - col_means^2))
  }
