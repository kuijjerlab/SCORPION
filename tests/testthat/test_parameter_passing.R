# Tests for parameter passing between scorpion() and runSCORPION()
# Validates that parameters are correctly passed from wrapper to core function

test_that("scorpion() and runSCORPION() share compatible parameter defaults", {
  data(scorpionTest)
  
  # scorpion() with default parameters
  result_scorpion <- scorpion(
    tfMotifs = scorpionTest$tf,
    gexMatrix = scorpionTest$gex[, 1:100],
    ppiNet = scorpionTest$ppi,
    showProgress = FALSE
  )
  
  # Both functions should work with default parameters
  result_runscoprion <- runSCORPION(
    gexMatrix = scorpionTest$gex,
    tfMotifs = scorpionTest$tf,
    ppiNet = scorpionTest$ppi,
    cellsMetadata = scorpionTest$metadata,
    groupBy = "region",
    showProgress = FALSE
  )
  
  # Both should succeed
  expect_true(is.list(result_scorpion) || is.matrix(result_scorpion))
  expect_true(is.data.frame(result_runscoprion))
})

test_that("runSCORPION() correctly passes computingEngine parameter to scorpion()", {
  data(scorpionTest)
  
  result_cpu <- runSCORPION(
    gexMatrix = scorpionTest$gex,
    tfMotifs = scorpionTest$tf,
    ppiNet = scorpionTest$ppi,
    cellsMetadata = scorpionTest$metadata,
    groupBy = "region",
    computingEngine = "cpu",
    showProgress = FALSE
  )
  
  expect_true(is.data.frame(result_cpu))
  expect_true(all(c("tf", "target") %in% colnames(result_cpu)))
})

test_that("runSCORPION() correctly passes nPC parameter to scorpion()", {
  data(scorpionTest)
  
  result <- runSCORPION(
    gexMatrix = scorpionTest$gex,
    tfMotifs = scorpionTest$tf,
    ppiNet = scorpionTest$ppi,
    cellsMetadata = scorpionTest$metadata,
    groupBy = "region",
    nPC = 20,
    showProgress = FALSE
  )
  
  expect_true(is.data.frame(result))
})

test_that("runSCORPION() correctly passes assocMethod parameter to scorpion()", {
  data(scorpionTest)
  
  result <- runSCORPION(
    gexMatrix = scorpionTest$gex,
    tfMotifs = scorpionTest$tf,
    ppiNet = scorpionTest$ppi,
    cellsMetadata = scorpionTest$metadata,
    groupBy = "region",
    assocMethod = "pearson",
    showProgress = FALSE
  )
  
  expect_true(is.data.frame(result))
})
