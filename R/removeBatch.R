remove_batch <- function(X, batch) {
    batch <- as.factor(batch)

    # Design matrix
    H <- model.matrix(~batch)

    # Solve coefficients: beta = (H'H)^(-1) H' X'
    HtH_inv <- solve(crossprod(H))
    HtX <- crossprod(H, t(X))
    beta <- HtH_inv %*% HtX

    # Correction: H * beta, then subtract
    X <- t(t(X) - H %*% beta)

    # Clear memory
    gc()
    return(X)
}
