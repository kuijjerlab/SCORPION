rowSds <-
  function(X, na.rm = FALSE) {
    apply(X, 1, function(r) {
      sd(r, na.rm = na.rm)
    })
  }
