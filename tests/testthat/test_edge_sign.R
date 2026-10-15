# Tests for inferEdgeSign() using SERGIO-simulated expression with a known signed GRN.
# The GRN was sampled from the signed (mor) edges of scorpionTest$tf; see data-raw/sergio/.

fx <- load_sergio()
P <- sergio_prior(fx)
base <- run_sign(fx, P)
tfs <- fx$genes$gene[fx$genes$role != "target"]

# ---------------------------------------------------------------------------
# Positive tests: true signs are recovered
# ---------------------------------------------------------------------------

test_that("signs of true GRN edges are recovered with high accuracy", {
  truth <- true_sign(fx, base)
  called <- data_called(base) & !is.na(truth)
  expect_gte(sum(called), 200)
  expect_gte(mean(data_sign(base)[called] == truth[called]), 0.95)
})

test_that("both activating and repressing edges are recovered", {
  truth <- true_sign(fx, base)
  called <- data_called(base) & !is.na(truth)
  ds <- data_sign(base)
  expect_gte(sum(called & truth > 0 & ds > 0), 100)
  expect_gte(sum(called & truth < 0 & ds < 0), 10)
})

test_that("held-out prior signs are recovered from the data", {
  set.seed(2)
  mask <- fx$prior$inGRN & runif(nrow(P)) < 0.3
  Pm <- P
  Pm$mor[mask] <- 0
  e <- run_sign(fx, Pm)
  held <- e[paste(e$tf, e$target) %in% paste(P$source_genesymbol, P$target_genesymbol)[mask], ]
  rec <- held[held$evidence == "data", ]
  expect_true(all(held$priorSign == 0))
  expect_gte(nrow(rec), 30)
  expect_gte(mean(data_sign(rec) == true_sign(fx, rec)), 0.9)
})

test_that("signs are recovered when TF mRNA is missing from the data", {
  e <- run_sign(fx, P, gex = fx$gex[!rownames(fx$gex) %in% tfs, ])
  truth <- true_sign(fx, e)
  called <- data_called(e) & !is.na(truth)
  expect_gte(sum(called), 200)
  expect_gte(mean(data_sign(e)[called] == truth[called]), 0.95)
})

test_that("signed, partially signed and unsigned priors give the same output", {
  set.seed(8)
  signedIdx <- which(P$mor != 0)
  keep <- unlist(lapply(split(signedIdx, P$source_genesymbol[signedIdx]), function(i) i[sample.int(length(i), min(2, length(i)))]))
  P2 <- P
  P2$mor[-keep] <- 0
  P0 <- P
  P0$mor <- 0
  out <- list(signed = base, partial = run_sign(fx, P2), unsigned = run_sign(fx, P0))
  labels <- c("activating", "repressing", "conflicting", "unresolved", "not_testable")
  for (e in out) {
    expect_identical(colnames(e), colnames(base))
    expect_identical(lapply(e, class), lapply(base, class))
    expect_true(all(e$label %in% labels))
    expect_true(all(e$evidence %in% c("data", "data+prior", "data_vs_prior", "prior_only", "none")))
  }
})

test_that("signs are recovered from partially signed priors", {
  set.seed(8)
  signedIdx <- which(P$mor != 0)
  keep <- unlist(lapply(split(signedIdx, P$source_genesymbol[signedIdx]), function(i) i[sample.int(length(i), min(2, length(i)))]))
  P2 <- P
  P2$mor[-keep] <- 0
  e <- run_sign(fx, P2)
  truth <- true_sign(fx, e)
  called <- data_called(e) & !is.na(truth)
  expect_gte(sum(called), 100)
  expect_gte(mean(data_sign(e)[called] == truth[called]), 0.9)
})

test_that("signs are recovered from unsigned priors", {
  P0 <- P
  P0$mor <- 0
  e <- run_sign(fx, P0)
  truth <- true_sign(fx, e)
  called <- data_called(e) & !is.na(truth)
  expect_true(all(e$evidence %in% c("data", "none")))
  expect_gte(sum(called), 100)
  expect_gte(mean(data_sign(e)[called] == truth[called]), 0.9)
})

test_that("an edge's estimate does not depend on its own prior sign", {
  # Leave-one-target-out: flipping one edge's prior sign leaves its estimate unchanged
  # (other edges of the same TF may change, since their activity uses this sign)
  tested <- base[!is.na(base$estimate) & base$priorSign != 0, ]
  flip <- which(paste(P$source_genesymbol, P$target_genesymbol) == paste(tested$tf[1], tested$target[1]))
  Pf <- P
  Pf$mor[flip] <- -Pf$mor[flip]
  e <- run_sign(fx, Pf)
  k <- paste(P$source_genesymbol, P$target_genesymbol)[flip]
  a <- base[paste(base$tf, base$target) %in% k, ]
  b <- e[paste(e$tf, e$target) %in% k, ]
  expect_gt(sum(!is.na(a$estimate)), 0)
  expect_equal(b$estimate[order(b$group, b$tf, b$target)], a$estimate[order(a$group, a$tf, a$target)])

  # Axis-based TFs: the orientation of edge (i, j) ignores the expected sign of j
  Y <- as.matrix(log1p(fx$gex))
  std <- SCORPION:::.es_standardize(Y, fx$metadata$donor, NULL)
  pr <- SCORPION:::.es_prior(P)
  pr$mor[pr$tf %in% tfs[1:3]] <- 0
  st <- SCORPION:::.es_axes(SCORPION:::.es_structure(pr, rownames(Y), std$ok, 5), std$Z)
  t1 <- tfs[1]
  tg <- st$genes[st$targets[[t1]]]
  o1 <- SCORPION:::.es_orientation(st, rep(t1, length(tg)), tg)
  for (k in seq_along(tg)) {
    st2 <- st
    st2$expSign[[t1]][k] <- -10 * st2$expSign[[t1]][k]
    expect_equal(SCORPION:::.es_orientation(st2, t1, tg[k]), o1[k])
  }
})

test_that("wrong prior signs are flagged as conflicting", {
  set.seed(4)
  flip <- which(fx$prior$inGRN & runif(nrow(P)) < 0.15)
  Pf <- P
  Pf$mor[flip] <- -Pf$mor[flip]
  e <- run_sign(fx, Pf)
  ef <- e[paste(e$tf, e$target) %in% paste(P$source_genesymbol, P$target_genesymbol)[flip] & data_called(e), ]
  expect_gte(nrow(ef), 20)
  expect_gte(mean(ef$label == "conflicting"), 0.9)
  expect_gte(mean(data_sign(ef) == true_sign(fx, ef)), 0.9)
})

# ---------------------------------------------------------------------------
# Negative tests: no signal, no call
# ---------------------------------------------------------------------------

test_that("decoy edges to permuted null genes are rarely called", {
  set.seed(3)
  src <- sample(fx$genes$gene[fx$genes$role == "target"], 60)
  nul <- permute_within_donor(fx$gex[src, ], fx$metadata$donor)
  rownames(nul) <- paste0("NULL_", src)
  dec <- expand.grid(source_genesymbol = tfs, target_genesymbol = rownames(nul), stringsAsFactors = FALSE)
  dec$weight <- 1
  dec$mor <- 0
  e <- run_sign(fx, rbind(P, dec), gex = rbind(fx$gex, nul))
  isDec <- grepl("^NULL_", e$target)
  expect_gte(sum(isDec & !is.na(e$estimate)), 1000)
  expect_lte(mean(data_called(e)[isDec]), 0.03)
})

test_that("permuted expression yields almost no sign calls", {
  set.seed(5)
  e <- run_sign(fx, P, gex = permute_within_donor(fx$gex, fx$metadata$donor))
  tested <- !is.na(e$estimate)
  expect_gte(sum(tested), 1000)
  expect_lte(sum(data_called(e)), 0.01 * sum(tested))
})

test_that("edges are not_testable when a TF has too few targets", {
  ghost <- data.frame(source_genesymbol = "GHOST", target_genesymbol = unique(fx$grn$target)[1:4],
                      weight = 1, mor = 0, stringsAsFactors = FALSE)
  e <- run_sign(fx, rbind(P, ghost))
  g <- e[e$tf == "GHOST", ]
  expect_equal(nrow(g), 8)
  expect_true(all(g$label == "not_testable"))
  expect_true(all(g$reason == "tf_activity_unavailable"))
  expect_true(all(is.na(g$estimate)))
})

test_that("unsigned priors on permuted expression or null genes yield almost no calls", {
  P0 <- P
  P0$mor <- 0
  set.seed(5)
  e <- run_sign(fx, P0, gex = permute_within_donor(fx$gex, fx$metadata$donor))
  expect_gte(sum(!is.na(e$estimate)), 1000)
  expect_lte(sum(e$evidence == "data"), 0.01 * sum(!is.na(e$estimate)))

  set.seed(3)
  src <- sample(fx$genes$gene[fx$genes$role == "target"], 60)
  nul <- permute_within_donor(fx$gex[src, ], fx$metadata$donor)
  rownames(nul) <- paste0("NULL_", src)
  dec <- expand.grid(source_genesymbol = tfs, target_genesymbol = rownames(nul), stringsAsFactors = FALSE)
  dec$weight <- 1
  dec$mor <- 0
  e <- run_sign(fx, rbind(P0, dec), gex = rbind(fx$gex, nul))
  isDec <- grepl("^NULL_", e$target)
  expect_lte(mean(e$evidence[isDec] == "data"), 0.03)
})

test_that("collinear regulators are reported as unresolved, not guessed", {
  clone <- P[P$source_genesymbol == tfs[1], ]
  clone$source_genesymbol <- "CLONE"
  e <- run_sign(fx, rbind(P, clone))
  shared <- e[e$tf %in% c(tfs[1], "CLONE") & e$target %in% clone$target_genesymbol & !is.na(e$estimate), ]
  expect_gt(nrow(shared), 0)
  expect_true(all(shared$reason == "collinear"))
  expect_false(any(data_called(shared)))
})

test_that("prior-only labels are never presented as data evidence", {
  po <- base[base$evidence == "prior_only", ]
  expect_true(all(po$priorSign != 0))
  expect_true(all(po$label == ifelse(po$priorSign > 0, "activating", "repressing")))
  expect_true(all(is.na(po$pAdj) | po$pAdj >= 0.05 | po$signStability < 0.9 | po$vif > 10))
})

# ---------------------------------------------------------------------------
# Method internals
# ---------------------------------------------------------------------------

test_that("Gram-matrix fit equals lm() on explicit leave-one-target-out activities", {
  set.seed(6)
  Y <- as.matrix(log1p(fx$gex[, fx$metadata$cell_type == "A"][, 1:300]))
  donor <- fx$metadata$donor[fx$metadata$cell_type == "A"][1:300]
  std <- SCORPION:::.es_standardize(Y, donor, NULL)
  pr <- SCORPION:::.es_prior(P)
  # Mix of signed and axis-based TFs: drop signs of two TFs
  pr$mor[pr$tf %in% tfs[1:2]] <- 0
  st <- SCORPION:::.es_structure(pr, rownames(Y), std$ok, 5)
  st <- SCORPION:::.es_axes(st, std$Z)
  expect_true(all(c("signed", "axis") %in% st$mode[st$tfs]))
  targets <- unique(pr$target[pr$tf %in% st$tfs & pr$target %in% rownames(Y)[std$ok]])[1:15]
  regs <- lapply(setNames(targets, targets), function(t) match(sort(unique(pr$tf[pr$target == t & pr$tf %in% st$tfs & pr$tf != t])), st$tfs))
  fit <- SCORPION:::.es_fit(std$Z, std$q, st, regs, ridge = 0)
  S <- SCORPION:::.es_activity(std$Z, st)
  for (t in targets) {
    R <- regs[[t]]
    j <- match(t, rownames(std$Z))
    m <- as.numeric(st$M[R, j])
    X <- sapply(seq_along(R), function(k) (S[R[k], ] - m[k] * std$Z[j, ]) / sqrt(st$W[R[k]] - m[k]^2))
    lmFit <- summary(lm(std$Z[j, ] ~ 0 + X))
    f <- fit[fit$target == t, ]
    dfAdj <- sqrt((ncol(Y) - length(R)) / (ncol(Y) - std$q - length(R)))
    expect_equal(f$estimate, unname(lmFit$coefficients[, 1]), tolerance = 1e-6)
    expect_equal(f$SE, unname(lmFit$coefficients[, 2]) * dfAdj, tolerance = 1e-6)
  }
})

test_that("leaving the target out removes the self-inclusion bias", {
  set.seed(7)
  Y <- permute_within_donor(as.matrix(log1p(fx$gex)), fx$metadata$donor)
  std <- SCORPION:::.es_standardize(Y, fx$metadata$donor, NULL)
  pr <- SCORPION:::.es_prior(P)
  st <- SCORPION:::.es_structure(pr, rownames(Y), std$ok, 5)
  S <- SCORPION:::.es_activity(std$Z, st)
  a <- Matrix::summary(st$M)
  a <- a[a$x > 0 & std$ok[a$j], ]
  naive <- mapply(function(i, j) cor(std$Z[j, ], S[i, ]), a$i, a$j)
  loo <- mapply(function(i, j, x) cor(std$Z[j, ], S[i, ] - x * std$Z[j, ]), a$i, a$j, a$x)
  expect_gt(mean(naive), 0.05)
  expect_lt(abs(mean(loo)), 0.02)
})

test_that("TF activity never depends on the TF's own expression", {
  set.seed(10)
  Y <- as.matrix(log1p(fx$gex))
  std <- SCORPION:::.es_standardize(Y, fx$metadata$donor, NULL)
  pr <- SCORPION:::.es_prior(P)
  pr$mor[pr$tf %in% tfs[1:4]] <- 0
  st <- SCORPION:::.es_structure(pr, rownames(Y), std$ok, 5)
  Z2 <- std$Z
  Z2[tfs, ] <- matrix(rnorm(length(tfs) * ncol(Z2)), nrow = length(tfs))
  for (t in st$tfs) {
    # Only the TF's own row differs; its targets are unchanged
    Zt <- std$Z
    Zt[t, ] <- Z2[t, ]
    a1 <- SCORPION:::.es_activity(std$Z, SCORPION:::.es_axes(st, std$Z, orient = TRUE))[t, ]
    a2 <- SCORPION:::.es_activity(Zt, SCORPION:::.es_axes(st, Zt, orient = TRUE))[t, ]
    expect_equal(a1, a2)
  }
})

test_that("results are reproducible with the same seed", {
  e1 <- run_sign(fx, P, seed = 11)
  e2 <- run_sign(fx, P, seed = 11)
  expect_identical(e1, e2)
})

test_that("SCORPION edge scores are attached by group", {
  sc <- data.frame(tf = P$source_genesymbol, target = P$target_genesymbol, A = seq_len(nrow(P)),
                   stringsAsFactors = FALSE)
  e <- inferEdgeSign(fx$gex, P, fx$metadata, groupBy = "cell_type", sampleBy = "donor",
                     gammaValue = 5, nBoot = 2, scorpionEdges = sc, showProgress = FALSE)$edges
  a <- e[e$group == "A", ]
  expect_equal(a$score, sc$A[match(paste(a$tf, a$target), paste(sc$tf, sc$target))])
  expect_true(all(is.na(e$score[e$group == "B"])))
})

# ---------------------------------------------------------------------------
# Input handling
# ---------------------------------------------------------------------------

test_that("invalid inputs are rejected", {
  expect_error(inferEdgeSign(fx$gex, P, fx$metadata, groupBy = "missing_col"), "not found")
  expect_error(inferEdgeSign(fx$gex, P, fx$metadata[-1, ], groupBy = "cell_type"), "one row per column")
  bad <- P
  bad$mor[1] <- 2
  expect_error(inferEdgeSign(fx$gex, bad, fx$metadata, groupBy = "cell_type"), "mor must contain")
  expect_error(inferEdgeSign(fx$gex, P, fx$metadata, groupBy = "cell_type", alpha = 1.5), "alpha")
  expect_error(inferEdgeSign(fx$gex, P, fx$metadata, groupBy = "cell_type", stabilityThreshold = 0.3), "stabilityThreshold")
  expect_error(inferEdgeSign(fx$gex, P, fx$metadata, groupBy = "cell_type", edges = data.frame(a = 1)), "edges must have")
})

test_that("a prior without mor is accepted and treated as unsigned", {
  e <- inferEdgeSign(fx$gex, P[, 1:3], fx$metadata, groupBy = "cell_type", sampleBy = "donor",
                     gammaValue = 5, nBoot = 2, showProgress = FALSE)$edges
  expect_true(all(e$priorSign == 0))
  expect_false(any(e$evidence %in% c("prior_only", "data+prior", "data_vs_prior")))
})

test_that("the edges argument restricts reported edges only", {
  sub <- data.frame(tf = P$source_genesymbol[1:20], target = P$target_genesymbol[1:20])
  e <- run_sign(fx, P, edges = sub)
  expect_true(all(paste(e$tf, e$target) %in% paste(sub$tf, sub$target)))
  b <- base[paste(base$tf, base$target) %in% paste(sub$tf, sub$target), ]
  expect_equal(e$estimate[order(e$group, e$tf, e$target)], b$estimate[order(b$group, b$tf, b$target)])
})
