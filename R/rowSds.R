rowSds <-
  function(X, na.rm = FALSE) {
    n <- rowSums(!is.na(X))
    row_means <- rowMeans(X, na.rm = na.rm)
    row_mean_sq <- rowMeans(X^2, na.rm = na.rm)
    sqrt(n / (n - 1) * (row_mean_sq - row_means^2))
  }
