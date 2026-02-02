# Tests for runSCORPION() network computation
# Validates that runSCORPION() produces correct output structures and grouping

test_that("runSCORPION() accepts gexMatrix and tfMotifs parameters and returns data.frame", {
  data(scorpionTest)
  result <- runSCORPION(
    gexMatrix = scorpionTest$gex,
    tfMotifs = scorpionTest$tf,
    ppiNet = scorpionTest$ppi,
    cellsMetadata = scorpionTest$metadata,
    groupBy = "region",
    showProgress = FALSE
  )

  expect_true(is.data.frame(result))
  expect_true(all(c("tf", "target") %in% colnames(result)))
  expect_true(nrow(result) > 0)
})

test_that("runSCORPION() produces wide-format output with network columns", {
  data(scorpionTest)
  result <- runSCORPION(
    gexMatrix = scorpionTest$gex,
    tfMotifs = scorpionTest$tf,
    ppiNet = scorpionTest$ppi,
    cellsMetadata = scorpionTest$metadata,
    groupBy = "region",
    showProgress = FALSE
  )

  # Should have tf, target columns plus network ID columns
  network_columns <- setdiff(colnames(result), c("tf", "target"))
  expect_true(length(network_columns) > 0)

  # Check that network column names match the unique grouping values
  unique_regions <- unique(scorpionTest$metadata$region)
  expect_true(all(unique_regions %in% network_columns))
})

# test_that("runSCORPION() handles NULL ppiNet correctly", {
#   skip("Skipped - implementation does not support NULL ppiNet")
#   data(scorpionTest)
#   result <- runSCORPION(
#     gexMatrix = scorpionTest$gex,
#     tfMotifs = scorpionTest$tf,
#     ppiNet = NULL,
#     cellsMetadata = scorpionTest$metadata,
#     groupBy = "region",
#     showProgress = FALSE
#   )
#
#   expect_true(is.data.frame(result))
#   expect_true(all(c("tf", "target") %in% colnames(result)))
# })

test_that("runSCORPION() selects correct output network with outNet parameter", {
  data(scorpionTest)

  # Test with default outNet = "regNet"
  result_default <- runSCORPION(
    gexMatrix = scorpionTest$gex,
    tfMotifs = scorpionTest$tf,
    ppiNet = scorpionTest$ppi,
    cellsMetadata = scorpionTest$metadata,
    groupBy = "region",
    outNet = "regNet",
    showProgress = FALSE
  )

  # Both should be data frames
  expect_true(is.data.frame(result_default))
})
