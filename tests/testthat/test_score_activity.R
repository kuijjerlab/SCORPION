# Tests for scoreActivity(): per-cell TF activity from a signed network (visualization)

toy <- function() {
  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")  # independent of other tests' RNG state
  a <- rnorm(200)                                   # hidden TF activity per cell
  expr <- rbind(up1 = 5 + a, up2 = 5 + a, up3 = 5 + a, down1 = 5 - a, down2 = 5 - a, other = rnorm(200, 5))
  expr <- expr + matrix(rnorm(length(expr), sd = 0.3), nrow(expr))
  colnames(expr) <- paste0("c", 1:200)
  # Activated and repressed weights balance, so unsigned weights cancel the signal
  W <- matrix(c(2, 1, 1, -2, -2, 0), nrow = 1, dimnames = list("TF", rownames(expr)))
  list(a = a, expr = expr, W = W)
}

test_that("scores equal the signed weighted mean of standardized targets", {
  d <- toy()
  s <- scoreActivity(d$expr, d$W, normalizeData = FALSE, minTargets = 3)
  Z <- t(scale(t(d$expr)))
  expect_equal(as.numeric(s["TF", ]), as.numeric(d$W %*% Z) / sum(abs(d$W)), tolerance = 1e-10)
  expect_identical(colnames(s), colnames(d$expr))
  expect_equal(unname(attr(s, "nTargets")), 5)
})

test_that("signed weights recover activity that unsigned weights blur", {
  d <- toy()
  signedScore <- scoreActivity(d$expr, d$W, normalizeData = FALSE, minTargets = 3)["TF", ]
  unsignedScore <- scoreActivity(d$expr, abs(d$W), normalizeData = FALSE, minTargets = 3)["TF", ]
  expect_gt(cor(signedScore, d$a), 0.95)
  expect_lt(abs(cor(unsignedScore, d$a)), 0.3)
})

test_that("NA and zero weights are ignored; TFs with too few targets are dropped", {
  d <- toy()
  W <- rbind(d$W, TF2 = c(1, NA, 0, 0, 0, 0))
  s <- scoreActivity(d$expr, W, normalizeData = FALSE, minTargets = 3)
  expect_identical(rownames(s), "TF")
  W2 <- d$W; W2[1, "up2"] <- NA
  s2 <- scoreActivity(d$expr, W2, normalizeData = FALSE, minTargets = 3)
  Z <- t(scale(t(d$expr)))
  w <- W2[1, ]; w[is.na(w)] <- 0
  expect_equal(as.numeric(s2["TF", ]), as.numeric(w %*% Z) / sum(abs(w)), tolerance = 1e-10)
})

test_that("a data.frame network (runSCORPION layout) gives the same scores as the matrix", {
  d <- toy()
  df <- data.frame(tf = "TF", target = colnames(d$W), N = as.numeric(d$W[1, ]))
  expect_equal(scoreActivity(d$expr, df, weightColumn = "N", normalizeData = FALSE, minTargets = 3),
               scoreActivity(d$expr, d$W, normalizeData = FALSE, minTargets = 3))
  expect_error(scoreActivity(d$expr, df, normalizeData = FALSE), "weightColumn")
})

test_that("scoreActivity works on scorpion() regSign output", {
  data(scorpionTest)
  cellsN <- scorpionTest$metadata$region == "N"
  set.seed(1)
  net <- scorpion(tfMotifs = scorpionTest$tf, gexMatrix = scorpionTest$gex[, cellsN], ppiNet = scorpionTest$ppi,
                  alphaValue = 0.8, outNet = c("regNet", "regSign"), showProgress = FALSE, signParams = list(nBoot = 2))
  s <- scoreActivity(scorpionTest$gex[, cellsN], net$regSign)
  expect_true(is.matrix(s))
  expect_identical(colnames(s), colnames(scorpionTest$gex)[cellsN])
  expect_true(all(rownames(s) %in% rownames(net$regSign)))
  expect_false(anyNA(s))
})

test_that("invalid inputs are rejected", {
  d <- toy()
  expect_error(scoreActivity(unname(d$expr), d$W), "rownames")
  expect_error(scoreActivity(d$expr, list(a = 1)), "matrix or a data.frame")
  W <- d$W; colnames(W) <- paste0("x", 1:6)
  expect_error(scoreActivity(d$expr, W), "No network targets")
  expect_error(scoreActivity(d$expr, d$W, normalizeData = FALSE, minTargets = 10), "at least 10")
})
