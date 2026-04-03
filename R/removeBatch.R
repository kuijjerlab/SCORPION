remove_batch <- function(X, batch) {
    batch <- as.factor(batch)

    # Design matrix
    H <- model.matrix(~batch)

    # Solve coefficients: beta = (H'H)^(-1) H' X'
    HtH_inv <- solve(crossprod(H))
    beta <- HtH_inv %*% t(X %*% H)

    # Correction: subtract H %*% beta from X
    X <- X - tcrossprod(t(beta), H)

    gc()
    return(X)
}
