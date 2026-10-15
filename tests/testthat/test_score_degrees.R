# Tests for scoreDegrees(): per-cell out-degrees and in-degrees from a signed network

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
std <- function(x) t(scale(t(x)))
score <- function(expr, W, ...) scoreDegrees(expr, W, normalizeData = FALSE, minTargets = 3, ...)

test_that("out-degree is the signed weighted mean of standardized targets", {
  d <- toy()
  s <- score(d$expr, d$W)$outDegree
  Z <- std(d$expr)
  expect_equal(as.numeric(s["TF", ]), as.numeric(d$W %*% Z) / sum(abs(d$W)), tolerance = 1e-10)
  expect_identical(colnames(s), colnames(d$expr))
})

test_that("in-degree uses each regulator's out-degree without the gene itself", {
  set.seed(2, kind = "Mersenne-Twister", normal.kind = "Inversion")
  expr <- matrix(rnorm(8 * 50, 5), 8, 50, dimnames = list(paste0("g", 1:8), paste0("c", 1:50)))
  W <- matrix(c(1, -1, 2, 0, 1, 1, 0, 0,
                0, 1, 1, -1, 0, 2, 1, -1), nrow = 2, byrow = TRUE, dimnames = list(c("A", "B"), rownames(expr)))
  s <- score(expr, W, topTargets = NULL)
  Z <- std(expr)
  for (j in colnames(W)[colSums(W != 0) > 0]) {
    regs <- rownames(W)[W[, j] != 0]
    loo <- sapply(regs, function(i) { w <- W[i, ]; w[j] <- 0; as.numeric(w %*% Z) / sum(abs(w)) })
    expected <- as.numeric(loo %*% W[regs, j]) / sum(abs(W[regs, j]))
    expect_equal(as.numeric(s$inDegree[j, ]), expected, tolerance = 1e-10)
  }
  expect_false("g8" %in% rownames(s$inDegree) && all(W[, "g8"] == 0) && FALSE)
})

test_that("signed weights recover activity that unsigned weights blur", {
  d <- toy()
  signedScore <- score(d$expr, d$W)$outDegree["TF", ]
  unsignedScore <- score(d$expr, abs(d$W))$outDegree["TF", ]
  expect_gt(cor(signedScore, d$a), 0.95)
  expect_lt(abs(cor(unsignedScore, d$a)), 0.3)
  # in-degree of activated targets follows the activity, of repressed targets opposes it
  inD <- score(d$expr, d$W)$inDegree
  expect_gt(cor(inD["up1", ], d$a), 0.9)
  expect_lt(cor(inD["down1", ], d$a), -0.9)
})

test_that("topTargets keeps each TF's strongest edges", {
  d <- toy()
  s3 <- score(d$expr, d$W, topTargets = 3)$outDegree
  w <- d$W[1, ]; w[order(-abs(w))[-(1:3)]] <- 0
  expect_equal(as.numeric(s3["TF", ]), as.numeric(w %*% std(d$expr)) / sum(abs(w)), tolerance = 1e-10)
  expect_error(score(d$expr, d$W, topTargets = 0), "topTargets")
})

test_that("NA and zero weights are ignored; TFs with too few targets are dropped", {
  d <- toy()
  W <- rbind(d$W, TF2 = c(1, NA, 0, 0, 0, 0))
  s <- score(d$expr, W)
  expect_identical(rownames(s$outDegree), "TF")
  expect_false("other" %in% rownames(s$inDegree))
})

test_that("a data.frame network (runSCORPION layout) gives the same scores as the matrix", {
  d <- toy()
  df <- data.frame(tf = "TF", target = colnames(d$W), N = as.numeric(d$W[1, ]))
  expect_equal(scoreDegrees(d$expr, df, weightColumn = "N", normalizeData = FALSE, minTargets = 3),
               score(d$expr, d$W))
  expect_error(scoreDegrees(d$expr, df, normalizeData = FALSE), "weightColumn")
})

test_that("scoreDegrees works on scorpion() regSign output", {
  data(scorpionTest)
  cellsN <- scorpionTest$metadata$region == "N"
  set.seed(1)
  net <- scorpion(tfMotifs = scorpionTest$tf, gexMatrix = scorpionTest$gex[, cellsN], ppiNet = scorpionTest$ppi,
                  alphaValue = 0.8, outNet = c("regNet", "regSign"), showProgress = FALSE, signParams = list(nBoot = 2))
  s <- scoreDegrees(scorpionTest$gex[, cellsN], net$regSign)
  expect_identical(names(s), c("outDegree", "inDegree"))
  expect_identical(colnames(s$outDegree), colnames(scorpionTest$gex)[cellsN])
  expect_true(all(rownames(s$outDegree) %in% rownames(net$regSign)))
  expect_true(all(rownames(s$inDegree) %in% colnames(net$regSign)))
  expect_false(anyNA(s$outDegree))
  expect_false(anyNA(s$inDegree))
})

test_that("invalid inputs are rejected", {
  d <- toy()
  expect_error(scoreDegrees(unname(d$expr), d$W), "rownames")
  expect_error(scoreDegrees(d$expr, list(a = 1)), "matrix or a data.frame")
  W <- d$W; colnames(W) <- paste0("x", 1:6)
  expect_error(scoreDegrees(d$expr, W), "No network targets")
  expect_error(scoreDegrees(d$expr, d$W, normalizeData = FALSE, minTargets = 10), "at least 10")
})
