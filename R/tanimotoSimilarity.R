tanimoto <- function(X, Y) {
  nc = ncol(Y)
  nr = nrow(X)
  dm = c(nr, nc)

  Amat = (X %*% Y)
  Bmat = colSums(Y * Y)

  Bmat = rep(Bmat, each = nr)
  dim(Bmat) = dm
  #Bmat=matrix(rep(Bmat, each=nr), dm)

  Cmat = rowSums(X * X)
  Cmat = rep(Cmat, nc)
  dim(Cmat) = dm
  #Cmat=matrix(rep(Cmat, nc), dm)

  den = (Bmat + Cmat - abs(Amat))
  Amat = Amat / sqrt(den)

  return(Amat)
}
