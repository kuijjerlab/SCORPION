remove_batch <- function(X, batch) {
    batch <- droplevels(as.factor(batch))

    # With fewer than two batch levels there is nothing to correct
    if (nlevels(batch) < 2) {
        return(X)
    }

    # Design matrix
    H <- model.matrix(~batch)

    # Solve coefficients: beta = (H'H)^(-1) H' X'
    # Fall back to a Moore-Penrose pseudo-inverse when H'H is rank-deficient
    # (e.g. collinear or singleton batches) so the run is not aborted with a
    # Lapack 'system is exactly singular' error.
    HtH <- crossprod(H)
    HtH_inv <- tryCatch(
        solve(HtH),
        error = function(e) {
            s <- svd(HtH)
            tol <- max(dim(HtH)) * .Machine$double.eps * max(s$d)
            keep <- s$d > tol
            s$v[, keep, drop = FALSE] %*%
                ((1 / s$d[keep]) * t(s$u[, keep, drop = FALSE]))
        }
    )
    beta <- HtH_inv %*% t(X %*% H)

    # Correction: subtract H %*% beta from X
    X <- X - tcrossprod(t(beta), H)

    gc()
    return(X)
}
