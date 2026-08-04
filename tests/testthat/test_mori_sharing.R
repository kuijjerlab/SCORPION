# Tests for optional 'mori' shared-memory acceleration of runSCORPION().
# All tests skip unless mori (and the parallel stack) is installed.

test_that("runSCORPION() with mori matches results without mori", {
  skip_on_cran()
  skip_if_not_installed("mori")
  skip_if_not_installed("furrr")
  skip_if_not_installed("future")

  data(scorpionTest)

  # Groups must each have enough cells to build a network (>= 30 cells).
  groups <- table(scorpionTest$metadata$region)
  skip_if(sum(groups >= 30) < 2, "Need >= 2 groups with >= 30 cells for parallel test")

  old_opt <- getOption("scorpion.use_mori")
  on.exit(options(scorpion.use_mori = old_opt), add = TRUE)

  run <- function() {
    set.seed(1)
    runSCORPION(
      gexMatrix = scorpionTest$gex,
      tfMotifs = scorpionTest$tf,
      ppiNet = scorpionTest$ppi,
      cellsMetadata = scorpionTest$metadata,
      groupBy = "region",
      alphaValue = 0.8,
      nCores = 2L,
      showProgress = FALSE
    )
  }

  options(scorpion.use_mori = FALSE)
  res_plain <- run()

  options(scorpion.use_mori = TRUE)
  res_mori <- run()

  expect_equal(dim(res_plain), dim(res_mori))
  expect_equal(colnames(res_plain), colnames(res_mori))
  expect_equal(res_plain$tf, res_mori$tf)
  expect_equal(res_plain$target, res_mori$target)

  weight_cols <- setdiff(colnames(res_plain), c("tf", "target"))
  for (col in weight_cols) {
    expect_equal(res_plain[[col]], res_mori[[col]],
                 tolerance = 1e-10,
                 label = paste("weights for network", col))
  }
})

test_that("mori::share() emits a compact serialized reference", {
  skip_on_cran()
  skip_if_not_installed("mori")

  # The whole point of mori: a shared object serializes to a tiny reference
  # instead of a full copy, so it is cheap to broadcast to every worker.
  x <- as.data.frame(matrix(rnorm(2e5), ncol = 4))
  shared <- mori::share(x)

  size_plain <- length(serialize(x, NULL))
  size_shared <- length(serialize(shared, NULL))

  expect_lt(size_shared, size_plain)
  expect_lt(size_shared, 10000)  # a few bytes / KB, not the full payload
})

test_that("runSCORPION() parallel performance with vs without mori", {
  skip_on_cran()
  skip_if_not_installed("mori")
  skip_if_not_installed("furrr")
  skip_if_not_installed("future")

  data(scorpionTest)

  groups <- table(scorpionTest$metadata$region)
  skip_if(sum(groups >= 30) < 2, "Need >= 2 groups with >= 30 cells for parallel test")

  old_opt <- getOption("scorpion.use_mori")
  on.exit(options(scorpion.use_mori = old_opt), add = TRUE)

  run <- function() {
    runSCORPION(
      gexMatrix = scorpionTest$gex,
      tfMotifs = scorpionTest$tf,
      ppiNet = scorpionTest$ppi,
      cellsMetadata = scorpionTest$metadata,
      groupBy = "region",
      alphaValue = 0.8,
      nCores = 2L,
      showProgress = FALSE
    )
  }

  # Warm up worker processes so plan startup cost is not attributed to a run.
  invisible(run())

  options(scorpion.use_mori = FALSE)
  t_plain <- system.time(run())[["elapsed"]]

  options(scorpion.use_mori = TRUE)
  t_mori <- system.time(run())[["elapsed"]]

  message(sprintf(
    "runSCORPION parallel timing: without mori = %.2fs, with mori = %.2fs (%.2fx)",
    t_plain, t_mori, t_plain / t_mori
  ))

  # Both paths must complete successfully; timing is reported, not asserted,
  # since wall-clock speedup depends on prior size, worker count and hardware.
  expect_true(is.finite(t_plain) && t_plain >= 0)
  expect_true(is.finite(t_mori) && t_mori >= 0)
})
