colSds <-
  function(X, na.rm = FALSE) {
    apply(X, 2, function(r) {
      sd(r, na.rm = na.rm)
    })
  }
