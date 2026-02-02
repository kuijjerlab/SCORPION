# Tests for data preprocessing and parameter effects in runSCORPION()
# Validates normalization, batch correction, and parameter passing

test_that("runSCORPION() parameter passing to scorpion() works with alphaValue", {
  data(scorpionTest)
  result_alpha_low <- runSCORPION(
    gexMatrix = scorpionTest$gex,
    tfMotifs = scorpionTest$tf,
    ppiNet = scorpionTest$ppi,
    cellsMetadata = scorpionTest$metadata,
    groupBy = "region",
    alphaValue = 0.1,
    showProgress = FALSE
  )
  
  result_alpha_high <- runSCORPION(
    gexMatrix = scorpionTest$gex,
    tfMotifs = scorpionTest$tf,
    ppiNet = scorpionTest$ppi,
    cellsMetadata = scorpionTest$metadata,
    groupBy = "region",
    alphaValue = 0.8,
    showProgress = FALSE
  )
  
  # Results should differ with different alphaValue
  expect_false(identical(result_alpha_low, result_alpha_high))
  
  # But should have same structure
  expect_equal(colnames(result_alpha_low), colnames(result_alpha_high))
  expect_equal(nrow(result_alpha_low), nrow(result_alpha_high))
})

test_that("runSCORPION() data normalization affects results", {
  data(scorpionTest)
  result_normalized <- runSCORPION(
    gexMatrix = scorpionTest$gex,
    tfMotifs = scorpionTest$tf,
    ppiNet = scorpionTest$ppi,
    cellsMetadata = scorpionTest$metadata,
    groupBy = "region",
    normalizeData = TRUE,
    showProgress = FALSE
  )
  
  result_unnormalized <- runSCORPION(
    gexMatrix = scorpionTest$gex,
    tfMotifs = scorpionTest$tf,
    ppiNet = scorpionTest$ppi,
    cellsMetadata = scorpionTest$metadata,
    groupBy = "region",
    normalizeData = FALSE,
    showProgress = FALSE
  )
  
  # Results should differ with different normalization
  expect_false(identical(result_normalized, result_unnormalized))
  
  # But should have same structure
  expect_equal(colnames(result_normalized), colnames(result_unnormalized))
})

test_that("runSCORPION() batch effect correction works with valid batch", {
  data(scorpionTest)
  
  result <- runSCORPION(
    gexMatrix = scorpionTest$gex,
    tfMotifs = scorpionTest$tf,
    ppiNet = scorpionTest$ppi,
    cellsMetadata = scorpionTest$metadata,
    groupBy = "region",
    removeBatchEffect = TRUE,
    batch = scorpionTest$metadata$donor,
    showProgress = FALSE
  )
  
  expect_true(is.data.frame(result))
  expect_true(all(c("tf", "target") %in% colnames(result)))
})
