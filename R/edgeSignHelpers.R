# Internal helpers for inferEdgeSign(). Notation follows the method description:
# Z (genes x observations) is expression residualized on donor/covariates and
# standardized per gene; M holds the target weights of each TF's activity and
# S = M Z is the un-normalized activity. How a TF's activity is built depends only
# on what the prior provides, never on user flags, and the output is the same:
#   - more than minTargets signed targets: m_ij = mor_ij * w_ij (fixed);
#   - otherwise: m_ij = w_ij * u_ij, with u_i the leading axis of the TF's targets,
#     cross-fitted over 5 folds (each fold's activity comes from axes learned on the
#     other folds), and oriented per edge by the expected signs of the TF's other
#     targets (mor when signed, otherwise the mean sign of the signed prior edges, or 1
#     if there are none).
# TF mRNA is never used: under SCORPION it reflects a future, not current, state.

# Collapse duplicated TF-target pairs: weights summed, sign kept only if unambiguous
.es_prior <- function(tfMotifs) {
  tfMotifs <- as.data.frame(tfMotifs)
  if (ncol(tfMotifs) < 3) {
    cli::cli_abort("tfMotifs must have at least 3 columns: TF, target and weight")
  }
  mor <- if ("mor" %in% colnames(tfMotifs)) {
    tfMotifs$mor
  } else if (ncol(tfMotifs) >= 4) {
    tfMotifs[[4]]
  } else {
    cli::cli_alert_warning("tfMotifs has no mor column; all prior signs are treated as unknown (0)")
    rep(0, nrow(tfMotifs))
  }
  mor <- suppressWarnings(as.numeric(mor))
  if (anyNA(mor) || any(!mor %in% c(-1, 0, 1))) {
    cli::cli_abort("mor must contain only -1, 0 or 1")
  }
  prior <- data.frame(
    tf = as.character(tfMotifs[[1]]),
    target = as.character(tfMotifs[[2]]),
    weight = as.numeric(tfMotifs[[3]]),
    mor = mor,
    stringsAsFactors = FALSE
  )
  if (any(!is.finite(prior$weight)) || any(prior$weight < 0)) {
    cli::cli_abort("tfMotifs weights must be finite and non-negative")
  }
  key <- paste(prior$tf, prior$target, sep = "\r")
  if (anyDuplicated(key)) {
    weight <- tapply(prior$weight, key, sum)
    mor <- tapply(prior$mor, key, function(v) {
      v <- unique(v[v != 0])
      if (length(v) == 1) v else 0
    })
    prior <- prior[!duplicated(key), ]
    key <- key[!duplicated(key)]
    prior$weight <- as.numeric(weight[key])
    prior$mor <- as.numeric(mor[key])
  }
  rownames(prior) <- NULL
  prior
}

# Aggregate cells into super-cells within each donor so that observations never mix donors.
# Covariates are appended as extra rows so they are averaged with the same membership.
.es_aggregate <- function(X, donor, covMat, gamma, nPC, seed) {
  if (gamma <= 1) {
    return(list(Y = as.matrix(X), donor = donor, cov = covMat))
  }
  geneNames <- rownames(X)
  out <- lapply(split(seq_len(ncol(X)), donor), function(idx) {
    nSC <- round(length(idx) / gamma)
    if (nSC < 2) {
      return(NULL)
    }
    Xd <- X[, idx, drop = FALSE]
    gv <- apply(as.matrix(Xd), 1, stats::var)
    genesUse <- names(sort(gv[gv > 0], decreasing = TRUE))[seq_len(min(1000, sum(gv > 0)))]
    if (!is.null(covMat)) {
      Xd <- rbind(Xd, t(covMat[idx, , drop = FALSE]))
    }
    SC <- suppressWarnings(makeSuperCells(
      X = Xd, genes.use = genesUse, gamma = gamma,
      n.pc = max(2, min(nPC, length(genesUse) - 1, length(idx) - 2)),
      fast.pca = FALSE, seed = seed
    ))
    SC <- as.matrix(SC)
    list(
      Y = SC[geneNames, , drop = FALSE],
      cov = if (is.null(covMat)) NULL else t(SC[colnames(covMat), , drop = FALSE]),
      donor = rep(donor[idx[1]], ncol(SC))
    )
  })
  out <- out[!vapply(out, is.null, logical(1))]
  if (length(out) == 0) {
    return(NULL)
  }
  list(
    Y = do.call(cbind, lapply(out, `[[`, "Y")),
    donor = unlist(lapply(out, `[[`, "donor"), use.names = FALSE),
    cov = if (is.null(covMat)) NULL else do.call(rbind, lapply(out, `[[`, "cov"))
  )
}

# Residualize on intercept + donor fixed effects + covariates (Frisch-Waugh-Lovell)
# and standardize each gene; returns Z and the residual degrees of freedom used.
.es_standardize <- function(Y, donor, cov) {
  C <- if (length(unique(donor)) > 1) {
    stats::model.matrix(~ factor(donor))
  } else {
    matrix(1, nrow = ncol(Y), ncol = 1)
  }
  if (!is.null(cov)) {
    C <- cbind(C, cov)
  }
  qrC <- qr(C)
  Q <- qr.Q(qrC)[, seq_len(qrC$rank), drop = FALSE]
  R <- Y - (Y %*% Q) %*% t(Q)
  ss <- rowSums(R^2)
  ok <- ss > 1e-10 * pmax(1, rowSums(Y^2))
  Z <- R
  Z[ok, ] <- R[ok, ] / sqrt(ss[ok] / (ncol(Y) - qrC$rank))
  Z[!ok, ] <- 0
  list(Z = Z, q = qrC$rank, ok = ok)
}

# Targets, weights and expected signs of each TF, and how its activity is built
.es_structure <- function(prior, genes, okGenes, minTargets) {
  piBar <- if (any(prior$mor != 0)) mean(prior$mor[prior$mor != 0]) else 1
  pr <- prior[prior$target %in% genes[okGenes] & prior$tf != prior$target, ]
  tfs <- sort(unique(prior$tf))
  f <- factor(pr$tf, levels = tfs)
  nTarg <- tabulate(f, nbins = length(tfs))
  nSigned <- vapply(split(pr$mor != 0, f), sum, numeric(1))
  mode <- ifelse(nSigned > minTargets, "signed", ifelse(nTarg > minTargets, "axis", "none"))
  names(mode) <- tfs
  testable <- tfs[mode != "none"]
  sg <- pr[pr$tf %in% testable[mode[testable] == "signed"] & pr$mor != 0, ]
  ax <- pr[pr$tf %in% testable[mode[testable] == "axis"], ]
  axTF <- testable[mode[testable] == "axis"]
  fa <- factor(ax$tf, levels = axTF)
  st <- list(
    tfs = testable, genes = genes, mode = mode,
    signedI = match(sg$tf, testable), signedJ = match(sg$target, genes), signedX = sg$mor * sg$weight,
    targets = split(match(ax$target, genes), fa),
    weights = split(ax$weight, fa),
    expSign = split(ifelse(ax$mor != 0, ax$mor, piBar), fa)
  )
  .es_setM(st, lapply(st$weights, function(w) rep(1, length(w))))
}

# Assemble M (and W = rowSums(M^2)) from signed weights and axis loadings
.es_setM <- function(st, loadings) {
  axTF <- names(st$targets)
  i <- c(st$signedI, rep(match(axTF, st$tfs), lengths(st$targets)))
  j <- c(st$signedJ, unlist(st$targets, use.names = FALSE))
  x <- c(st$signedX, unlist(Map(`*`, st$weights, loadings[axTF]), use.names = FALSE))
  st$M <- Matrix::sparseMatrix(i = i, j = j, x = x, dims = c(length(st$tfs), length(st$genes)),
                               dimnames = list(st$tfs, st$genes))
  st$W <- Matrix::rowSums(st$M^2)
  st$loadings <- loadings
  st
}

# Leading axis of the weighted targets of each axis-based TF in Z. The sign is
# arbitrary unless `orient` is TRUE, in which case it follows the expected signs.
.es_axes <- function(st, Z, orient = FALSE, ref = NULL) {
  pc1 <- function(X) {
    if (min(dim(X)) > 50) irlba::irlba(X, nv = 1, nu = 1)$u[, 1] else svd(X, nu = 1, nv = 0)$u[, 1]
  }
  loadings <- lapply(stats::setNames(names(st$targets), names(st$targets)), function(t) {
    X <- Z[st$targets[[t]], , drop = FALSE] * st$weights[[t]]
    u <- pc1(X)
    if (!is.null(ref)) {
      if (sum(u * ref[[t]]) < 0) -u else u
    } else if (orient && sum(st$expSign[[t]] * st$weights[[t]] * u) < 0) {
      -u
    } else {
      u
    }
  })
  .es_setM(st, loadings)
}

# Per-edge orientation: 1 for signed TFs; for axis-based TFs, from the expected
# signs of the TF's OTHER targets, o_ij = sign( sum_k e_ik m_ik - e_ij m_ij ), so the
# tested target (and its own prior sign) never orients its estimate.
.es_orientation <- function(st, tf, target) {
  o <- rep(1, length(tf))
  ax <- st$mode[tf] == "axis"
  if (any(ax)) {
    axTF <- names(st$targets)
    cTF <- vapply(axTF, function(t) sum(st$expSign[[t]] * st$weights[[t]] * st$loadings[[t]]), numeric(1))
    jj <- match(target[ax], st$genes)
    k <- mapply(function(t, j) match(j, st$targets[[t]]), tf[ax], jj)
    eij <- mapply(function(t, kk) st$expSign[[t]][kk], tf[ax], k)
    mij <- st$M[cbind(match(tf[ax], st$tfs), jj)]
    o[ax] <- sign(round(cTF[tf[ax]] - eij * mij, 12))
  }
  o
}

# TF activity S (un-normalized) for all testable TFs
.es_activity <- function(Z, st) {
  S <- as.matrix(st$M %*% Z)
  rownames(S) <- st$tfs
  S
}

# Joint per-target regression of z_j on leave-one-target-out activities of all of
# its testable regulators, computed from Gram matrices only. Observations may come in
# folds (cross-fitting), each with its own activity weights M_f; per fold
#   G_f = S_f S_f', H_f = S_f Z_f', nu_f = ||z_j,f||^2
# and the normal equations are the sums over folds. A single fold is the ordinary fit.
# Returns one row per (regulator, target) with estimate, SE, df and VIF.
.es_fitParts <- function(parts, n, q, tfs, regsByTarget, ridge = 1e-8) {
  pre <- lapply(parts, function(p) {
    list(G = tcrossprod(p$S), H = p$S %*% t(p$Z), nu = rowSums(p$Z^2), M = p$M, W = p$W)
  })
  genes <- rownames(parts[[1]]$Z)
  res <- lapply(names(regsByTarget), function(tg) {
    R <- regsByTarget[[tg]]
    j <- match(tg, genes)
    p <- length(R)
    df <- n - q - p
    if (df < 3) {
      return(data.frame(tf = tfs[R], target = tg, estimate = NA_real_, SE = NA_real_,
                        df = df, vif = NA_real_, stringsAsFactors = FALSE))
    }
    XtX <- matrix(0, p, p)
    Xty <- numeric(p)
    nuj <- 0
    for (f in pre) {
      m <- as.numeric(f$M[R, j])
      d <- sqrt(pmax(f$W[R] - m^2, 1e-12))
      h <- f$H[R, j]
      XtX <- XtX + (f$G[R, R, drop = FALSE] - outer(h, m) - outer(m, h) + outer(m, m) * f$nu[j]) / outer(d, d)
      Xty <- Xty + (h - m * f$nu[j]) / d
      nuj <- nuj + f$nu[j]
    }
    A <- XtX + diag(ridge * mean(diag(XtX)) + 1e-12, p)
    Ainv <- tryCatch(solve(A), error = function(e) NULL)
    if (is.null(Ainv)) {
      return(data.frame(tf = tfs[R], target = tg, estimate = NA_real_, SE = NA_real_,
                        df = df, vif = Inf, stringsAsFactors = FALSE))
    }
    beta <- as.numeric(Ainv %*% Xty)
    rss <- max(nuj - 2 * sum(beta * Xty) + as.numeric(t(beta) %*% XtX %*% beta), 0)
    covB <- (rss / df) * (Ainv %*% XtX %*% Ainv)
    vif <- if (p == 1) {
      1
    } else {
      dd <- sqrt(pmax(diag(XtX), 1e-12))
      Rm <- XtX / outer(dd, dd)
      tryCatch(diag(solve(Rm)), error = function(e) rep(Inf, p))
    }
    data.frame(tf = tfs[R], target = tg, estimate = beta, SE = sqrt(pmax(diag(covB), 0)),
               df = df, vif = vif, stringsAsFactors = FALSE)
  })
  do.call(rbind, res)
}

# Single-fold (full data) fit
.es_fit <- function(Z, q, st, regsByTarget, ridge = 1e-8) {
  .es_fitParts(list(list(Z = Z, S = .es_activity(Z, st), M = st$M, W = st$W)),
               ncol(Z), q, st$tfs, regsByTarget, ridge)
}

# Bootstrap resample of observations: donors with replacement when there are at
# least 5, otherwise observations with replacement within each donor.
.es_resample <- function(donor) {
  ud <- unique(donor)
  if (length(ud) >= 5) {
    pick <- sample(ud, length(ud), replace = TRUE)
    idx <- unlist(lapply(seq_along(pick), function(k) which(donor == pick[k])))
    lab <- unlist(lapply(seq_along(pick), function(k) rep(paste0(pick[k], "#", k), sum(donor == pick[k]))))
    return(list(idx = idx, donor = lab))
  }
  idx <- unlist(lapply(split(seq_along(donor), donor), function(i) i[sample.int(length(i), length(i), replace = TRUE)]))
  list(idx = idx, donor = donor[idx])
}

# Random split of observations into K folds, balanced within each donor
.es_folds <- function(donor, K) {
  f <- integer(length(donor))
  for (d in unique(donor)) {
    i <- which(donor == d)
    f[i] <- sample(rep(seq_len(K), length.out = length(i)))
  }
  f
}

# Edge effects. When every TF is signed, a single full-data fit is used. Otherwise
# axes are cross-fitted: for each of K folds, axes are learned on the other folds and
# give the activity of the fold's observations, so every observation enters the
# regression once with an activity it did not help estimate. Axis signs are aligned
# to the full-data axes, and estimates are oriented per edge from the expected signs
# of the TF's other targets.
.es_estimate <- function(Y, donor, cov, st, regsByTarget, edgesKey, K = 5) {
  key <- function(f) paste(f$tf, f$target)
  std <- .es_standardize(Y, donor, cov)
  if (length(st$targets) == 0) {
    fit <- .es_fit(std$Z, std$q, st, regsByTarget)
    return(fit[match(edgesKey, key(fit)), c("tf", "target", "estimate", "SE", "df", "vif")])
  }
  Z <- std$Z
  ref <- .es_axes(st, Z)
  folds <- .es_folds(donor, K)
  parts <- lapply(sort(unique(folds)), function(f) {
    stF <- .es_axes(st, Z[, folds != f, drop = FALSE], ref = ref$loadings)
    Zf <- Z[, folds == f, drop = FALSE]
    list(Z = Zf, S = .es_activity(Zf, stF), M = stF$M, W = stF$W)
  })
  fit <- .es_fitParts(parts, ncol(Z), std$q, st$tfs, regsByTarget)
  fit <- fit[match(edgesKey, key(fit)), c("tf", "target", "estimate", "SE", "df", "vif")]
  o <- .es_orientation(ref, fit$tf, fit$target)
  fit$estimate <- ifelse(o == 0, NA, o * fit$estimate)
  fit
}

# Run inferEdgeSign() on the cells of a single SCORPION network (no grouping).
# Prior edges are weighted by this network's own regNet strength, Phi^-1((1 + Phi(z)) / 2)
# cubed, so TF activities are built mainly from the targets PANDA supports in this
# network. Directions use single cells (gammaValue = 1) unless signParams says otherwise;
# the network itself is built with the gammaValue given to scorpion().
.scorpion_edgeSign <- function(gex, tfMotifs, regNet, zScaling, samples, nPC, signParams) {
  prior <- .es_prior(tfMotifs)
  prior <- prior[prior$tf %in% rownames(regNet) & prior$target %in% colnames(regNet), ]
  z <- regNet[cbind(prior$tf, prior$target)]
  if (!zScaling) {
    z <- stats::qnorm(pmin(pmax(z, 1e-12), 1 - 1e-12))
  }
  prior$weight <- (-stats::qnorm(stats::pnorm(-z, log.p = TRUE) - log(2), log.p = TRUE))^3
  prior <- prior[is.finite(prior$weight), ]
  meta <- data.frame(network = rep("network", ncol(gex)),
                     sample = if (is.null(samples)) "all" else as.character(samples),
                     stringsAsFactors = FALSE)
  args <- utils::modifyList(list(
    gexMatrix = gex, tfMotifs = prior, cellsMetadata = meta, groupBy = "network",
    sampleBy = if (is.null(samples)) NULL else "sample", normalizeData = FALSE,
    gammaValue = 1, nPC = nPC, minCells = min(30, ncol(gex)), showProgress = FALSE
  ), signParams)
  es <- tryCatch(do.call(inferEdgeSign, args)$edges, error = function(e) {
    cli::cli_alert_warning("Edge directions could not be estimated: {conditionMessage(e)}")
    NULL
  })
  if (!is.null(es)) {
    es$group <- NULL
  }
  es
}

# Signed regulatory network: d * Phi^{-1}((1 + Phi(z)) / 2) on the z scale, d * Phi(z)
# on the [0, 1] scale, computed from one network's cells only. The direction d of a prior
# edge is its label from inferEdgeSign() in this network: the data direction when called
# (including when it conflicts with the prior), otherwise the prior's expected sign, i.e.
# its mor when signed and, for unsigned edges, the prevailing sign of the prior's signed
# edges (+1 when the prior has no signs); this includes self-edges (autoregulation),
# whose direction cannot be estimated from the data. Pairs absent from the prior are
# NA. The value carries direction only, never certainty.
.es_signNetwork <- function(regNet, edgeSign, prior, zScaling) {
  piBar <- if (any(prior$mor != 0)) sign(mean(prior$mor[prior$mor != 0])) else 1
  if (piBar == 0) piBar <- 1
  d <- matrix(NA_real_, nrow(regNet), ncol(regNet), dimnames = dimnames(regNet))
  pr <- prior[prior$tf %in% rownames(d) & prior$target %in% colnames(d), ]
  d[cbind(match(pr$tf, rownames(d)), match(pr$target, colnames(d)))] <- ifelse(pr$mor != 0, pr$mor, piBar)
  if (!is.null(edgeSign)) {
    e <- edgeSign[edgeSign$evidence %in% c("data", "data+prior", "data_vs_prior") &
                    edgeSign$tf %in% rownames(d) & edgeSign$target %in% colnames(d), ]
    d[cbind(match(e$tf, rownames(d)), match(e$target, colnames(d)))] <-
      ifelse(e$label == "conflicting", -e$priorSign, ifelse(e$label == "activating", 1, -1))
  }
  if (zScaling) {
    # Phi^{-1}((1 + Phi(z)) / 2) = -Phi^{-1}(Phi(-z) / 2), computed on the log scale
    magnitude <- -stats::qnorm(stats::pnorm(-regNet, log.p = TRUE) - log(2), log.p = TRUE)
  } else {
    magnitude <- regNet
  }
  d * magnitude
}
