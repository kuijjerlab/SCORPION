# Tests for the signed regulatory network (regSign) returned by scorpion() and runSCORPION()

data(scorpionTest)
xN <- scorpionTest$metadata$region == "N"
fastSign <- list(nBoot = 5)

run_scorpion <- function(...) {
  set.seed(1)
  scorpion(tfMotifs = scorpionTest$tf, gexMatrix = scorpionTest$gex[, xN], ppiNet = scorpionTest$ppi,
           alphaValue = 0.8, showProgress = FALSE, signParams = fastSign, ...)
}

halfNormal <- function(z) -stats::qnorm(stats::pnorm(-z, log.p = TRUE) - log(2), log.p = TRUE)

base3 <- run_scorpion(outNet = c("regNet", "coregNet", "coopNet"))
withSign <- run_scorpion(samples = scorpionTest$metadata$donor[xN])
# Per-edge details behind regSign (internal; same inputs as scorpion(), so identical results)
edge_details <- function(gex, prior, regNet, samples, signParams) {
  SCORPION:::.scorpion_edgeSign(gex, prior, regNet, TRUE, samples, 25, signParams)
}
signInfo <- edge_details(scorpionTest$gex[, xN], scorpionTest$tf, withSign$regNet,
                         scorpionTest$metadata$donor[xN], fastSign)

test_that("regSign is the only added output of scorpion()", {
  expect_identical(names(withSign), c("regNet", "coregNet", "coopNet", "regSign", "numGenes", "numTFs", "numEdges"))
  expect_false("regSign" %in% names(base3))
  expect_identical(sort(names(attributes(withSign$regSign))), c("dim", "dimnames"))
})

test_that("requesting regSign leaves the PANDA networks unchanged", {
  expect_identical(withSign$regNet, base3$regNet)
  expect_identical(withSign$coregNet, base3$coregNet)
  expect_identical(withSign$coopNet, base3$coopNet)
  expect_identical(withSign[c("numGenes", "numTFs", "numEdges")], base3[c("numGenes", "numTFs", "numEdges")])
})

label_direction <- function(e) {
  ifelse(e$label == "activating", 1, ifelse(e$label == "repressing", -1,
         ifelse(e$label == "conflicting", -e$priorSign, NA)))
}

test_that("regSign is d * Phi^-1((1 + Phi(z)) / 2): called direction, else the prior's expected sign", {
  rs <- withSign$regSign
  expect_identical(dimnames(rs), dimnames(withSign$regNet))
  pr <- SCORPION:::.es_prior(scorpionTest$tf)
  pr <- pr[pr$tf %in% rownames(rs) & pr$target %in% colnames(rs), ]
  piBar <- sign(mean(pr$mor[pr$mor != 0]))
  d <- ifelse(pr$mor != 0, pr$mor, piBar)
  e <- signInfo$edges
  called <- e[e$evidence %in% c("data", "data+prior", "data_vs_prior"), ]
  m <- match(paste(pr$tf, pr$target), paste(called$tf, called$target))
  d[!is.na(m)] <- label_direction(called)[m[!is.na(m)]]
  idx <- cbind(pr$tf, pr$target)
  expect_equal(as.numeric(rs[idx]), d * halfNormal(withSign$regNet[idx]))
})

test_that("every edge is directed, weighted and signed; non-prior edges follow TF activity", {
  rs <- withSign$regSign
  expect_false(anyNA(rs))
  expect_false(any(rs == 0))
  expect_equal(abs(rs), halfNormal(withSign$regNet), check.attributes = FALSE)
  pr <- SCORPION:::.es_prior(scorpionTest$tf)
  piBar <- sign(mean(pr$mor[pr$mor != 0]))
  inPrior <- matrix(FALSE, nrow(rs), ncol(rs), dimnames = dimnames(rs))
  inPrior[cbind(pr$tf, pr$target)[pr$tf %in% rownames(rs) & pr$target %in% colnames(rs), , drop = FALSE]] <- TRUE
  # non-prior edges with an association: its sign; otherwise the prior's expected sign
  a <- matrix(NA_real_, nrow(rs), ncol(rs), dimnames = dimnames(rs))
  sa <- signInfo$association
  a[rownames(sa)[rownames(sa) %in% rownames(rs)], colnames(sa)[colnames(sa) %in% colnames(rs)]] <-
    sa[rownames(sa) %in% rownames(rs), colnames(sa) %in% colnames(rs)]
  useA <- !inPrior & is.finite(a) & a != 0
  expect_gt(sum(useA), 0)
  expect_true(all(sign(rs[useA]) == sign(a[useA])))
  expect_true(all(sign(rs[!inPrior & !useA]) == piBar))
})

test_that("the direction step uses single cells while the network uses gammaValue", {
  # single cells: residual degrees of freedom far above the number of super-cells
  expect_gt(max(signInfo$edges$df, na.rm = TRUE), sum(xN) / 2)
  e10 <- edge_details(scorpionTest$gex[, xN], scorpionTest$tf, withSign$regNet,
                      scorpionTest$metadata$donor[xN], list(nBoot = 2, gammaValue = 10))
  expect_lt(max(e10$edges$df, na.rm = TRUE), sum(xN) / 5)
  set.seed(1)
  r <- scorpion(tfMotifs = scorpionTest$tf, gexMatrix = scorpionTest$gex[, xN], ppiNet = scorpionTest$ppi,
                alphaValue = 0.8, showProgress = FALSE, outNet = c("regNet", "regSign"),
                signParams = list(nBoot = 2, gammaValue = 10))
  expect_identical(r$regNet, base3$regNet)
})

test_that("regSign direction is the called sign whenever an edge is called", {
  e <- signInfo$edges
  called <- e[e$evidence %in% c("data", "data+prior", "data_vs_prior"), ]
  skip_if(nrow(called) == 0)
  expect_true(all(sign(withSign$regSign[cbind(called$tf, called$target)]) == label_direction(called)))
})

test_that("regSign magnitude increases with the edge score", {
  z <- seq(-6, 12, by = 0.5)
  m <- halfNormal(z)
  expect_true(all(diff(m) > 0))
  expect_true(all(m >= 0))
  expect_lt(m[1], 1e-3)
  expect_true(all(is.finite(halfNormal(c(30, 40)))))
})

test_that("regSign is on [-1, 1] when zScaling = FALSE", {
  r <- run_scorpion(outNet = c("regNet", "regSign"), zScaling = FALSE)
  expect_true(all(r$regSign >= -1 & r$regSign <= 1, na.rm = TRUE))
  expect_true(all(abs(r$regSign) <= r$regNet + 1e-12, na.rm = TRUE))
})

test_that("regSign can be requested without regNet and with an unsigned prior", {
  r <- run_scorpion(outNet = "regSign")
  expect_false("regNet" %in% names(r))
  expect_identical(dimnames(r$regSign), dimnames(base3$regNet))
  set.seed(1)
  u <- scorpion(tfMotifs = scorpionTest$tf[, 1:3], gexMatrix = scorpionTest$gex[, xN], ppiNet = scorpionTest$ppi,
                alphaValue = 0.8, showProgress = FALSE, outNet = c("regNet", "regSign"), signParams = fastSign)
  expect_identical(u$regNet, base3$regNet)
  expect_false(anyNA(u$regSign))
})

test_that("runSCORPION() returns regSign rows only when requested", {
  set.seed(1)
  r1 <- runSCORPION(scorpionTest$gex, scorpionTest$tf, scorpionTest$ppi, scorpionTest$metadata,
                    groupBy = "region", alphaValue = 0.8, showProgress = FALSE)
  set.seed(1)
  r2 <- runSCORPION(scorpionTest$gex, scorpionTest$tf, scorpionTest$ppi, scorpionTest$metadata,
                    groupBy = "region", alphaValue = 0.8, showProgress = FALSE,
                    outNet = c("regNet", "regSign"), sampleBy = "donor", signParams = fastSign)
  expect_false("edge_type" %in% colnames(r1))
  expect_setequal(unique(r2$edge_type), c("tf-target", "tf-target-sign"))
  reg <- r2[r2$edge_type == "tf-target", -1]
  rownames(reg) <- NULL
  expect_equal(reg, r1, check.attributes = FALSE)
  sg <- r2[r2$edge_type == "tf-target-sign", ]
  pairs <- function(d) {
    d <- d[, c("tf", "target")]
    rownames(d) <- NULL
    d
  }
  expect_identical(pairs(sg), pairs(r2[r2$edge_type == "tf-target", ]))
  expect_false(anyNA(sg$T))
  expect_false(any(sg$T == 0))
})

test_that("invalid regSign arguments are rejected", {
  expect_error(scorpion(scorpionTest$tf, scorpionTest$gex[, xN], scorpionTest$ppi, outNet = "badNet",
                        showProgress = FALSE), "outNet")
  expect_error(scorpion(scorpionTest$tf, scorpionTest$gex[, xN], scorpionTest$ppi, samples = "D1",
                        showProgress = FALSE), "samples")
  expect_error(runSCORPION(scorpionTest$gex, scorpionTest$tf, scorpionTest$ppi, scorpionTest$metadata,
                           groupBy = "region", outNet = c("regNet", "regSign"), sampleBy = "nope",
                           showProgress = FALSE), "sampleBy")
})

# ---------------------------------------------------------------------------
# Ground truth (SERGIO): regSign directions match the simulated GRN
# ---------------------------------------------------------------------------

fx <- load_sergio()
P <- sergio_prior(fx)
tfsS <- fx$genes$gene[fx$genes$role != "target"]
ppiS <- data.frame(source_genesymbol = tfsS, target_genesymbol = tfsS, weight = 1)
cellsA <- fx$metadata$cell_type == "A"
gexA <- SCORPION:::log_normalize_data(fx$gex[, cellsA])

sergio_regsign <- function(prior, gex = gexA) {
  set.seed(1)
  r <- scorpion(tfMotifs = prior, gexMatrix = gex, ppiNet = ppiS, gammaValue = 5, showProgress = FALSE,
                outNet = c("regNet", "regSign"), samples = fx$metadata$donor[cellsA], signParams = list(nBoot = 30))
  r$edges <- edge_details(gex, prior, r$regNet, fx$metadata$donor[cellsA], list(nBoot = 30))$edges
  r
}

grn_agreement <- function(r) {
  g <- fx$grn[fx$grn$tf %in% rownames(r$regSign) & fx$grn$target %in% colnames(r$regSign), ]
  s <- sign(r$regSign[cbind(g$tf, g$target)])
  c(known = mean(!is.na(s)), agreement = mean(s[!is.na(s)] == g$mor[!is.na(s)]),
    negative = mean(s[!is.na(s)] < 0), majority = mean(g$mor > 0))
}

test_that("regSign recovers simulated directions from signed and unsigned priors", {
  P0 <- P
  P0$mor <- 0
  rS <- sergio_regsign(P)
  rU <- sergio_regsign(P0)
  # Called edges (unsigned priors cannot orient TFs whose targets are mostly repressed)
  for (k in 1:2) {
    r <- list(rS, rU)[[k]]
    e <- r$edges[data_called(r$edges), ]
    tr <- true_sign(fx, e)
    s <- sign(r$regSign[cbind(e$tf, e$target)])
    expect_gte(sum(!is.na(tr)), 50)
    expect_gte(mean(s[!is.na(tr)] == tr[!is.na(tr)]), c(0.95, 0.9)[k])
  }
  # All true edges: signed prior near-perfect; unsigned prior never worse than "all activating"
  expect_gte(grn_agreement(rS)[["agreement"]], 0.95)
  gu <- grn_agreement(rU)
  expect_gte(gu[["agreement"]], gu[["majority"]])
  expect_gt(gu[["negative"]], 0)
})

test_that("regSign carries no direction signal on permuted expression", {
  P0 <- P
  P0$mor <- 0
  set.seed(12)
  gexP <- SCORPION:::log_normalize_data(permute_within_donor(fx$gex[, cellsA], fx$metadata$donor[cellsA]))
  r <- sergio_regsign(P0, gex = gexP)
  expect_lte(sum(data_called(r$edges)), 5)
  # Without signal, unsigned edges keep the prior's expected sign: (almost) nothing is flipped
  expect_lte(grn_agreement(r)[["negative"]], 0.01)
})

test_that("groups of different size from one population do not differ in regSign", {
  # One population split at random into a large and a small group within each donor:
  # unknown directions must not be read as differences, so regSign adds no false
  # positives beyond those of regNet itself.
  meta <- fx$metadata[cellsA, ]
  # Deterministic split (every 5th cell of each donor), independent of the RNG state
  meta$grp <- ave(seq_len(nrow(meta)), meta$donor, FUN = function(i) ifelse(seq_along(i) %% 5 == 0, "small", "big"))
  big <- paste0("D", 1:4, "--big")
  small <- paste0("D", 1:4, "--small")
  for (pr in list(P, transform(P, mor = 0))) {
    set.seed(1)
    net <- runSCORPION(fx$gex[, cellsA], pr, ppiS, meta, groupBy = c("donor", "grp"),
                       outNet = c("regNet", "regSign"), gammaValue = 5, signParams = list(nBoot = 20),
                       showProgress = FALSE)
    sg <- net[net$edge_type == "tf-target-sign", -1]
    rg <- net[net$edge_type == "tf-target", -1]
    ip <- paste(sg$tf, sg$target) %in% paste(pr[, 1], pr[, 2])
    test <- function(d) suppressMessages(testEdges(d[ip, c("tf", "target", big, small)], testType = "two.sample",
                                                   group1 = big, group2 = small, paired = TRUE))
    ds <- test(sg)
    dr <- test(rg)
    fpSign <- paste(ds$tf, ds$target)[which(ds$pAdj < 0.05)]
    fpReg <- paste(dr$tf, dr$target)[which(dr$pAdj < 0.05)]
    expect_true(all(fpSign %in% fpReg))
    expect_lte(length(fpSign), length(fpReg))
  }
})
