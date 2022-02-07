dFunction <- function(X, Y) {
  A <- cor(X, Y)
  A[A < 0] <- 0
  A
}
