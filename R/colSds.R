colSds <-
  function(X, na.rm = FALSE) {
    apply(X, 2, function(r) {
      stats::sd(r, na.rm = na.rm)
    })
  }
