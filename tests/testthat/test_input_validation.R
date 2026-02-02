# Tests for input validation and error handling
# Validates that functions properly validate inputs and provide meaningful errors

test_that("scorpion() alphaValue parameter validation - negative value", {
  data(scorpionTest)
  
  # alphaValue < 0 should error
  expect_error(
    scorpion(
      tfMotifs = scorpionTest$tf,
      gexMatrix = scorpionTest$gex[, 1:50],
      ppiNet = scorpionTest$ppi,
      alphaValue = -0.1,
      showProgress = FALSE
    )
  )
})

test_that("scorpion() alphaValue parameter validation - value > 1", {
  data(scorpionTest)
  
  # alphaValue > 1 should error
  expect_error(
    scorpion(
      tfMotifs = scorpionTest$tf,
      gexMatrix = scorpionTest$gex[, 1:50],
      ppiNet = scorpionTest$ppi,
      alphaValue = 1.5,
      showProgress = FALSE
    )
  )
})

test_that("runSCORPION() input validation checks gexMatrix dimensions", {
  data(scorpionTest)
  
  # Create mismatched metadata
  bad_metadata <- scorpionTest$metadata[1:100, ]
  
  # Should error because gexMatrix has more columns than metadata has rows
  expect_error(
    runSCORPION(
      gexMatrix = scorpionTest$gex,
      tfMotifs = scorpionTest$tf,
      ppiNet = scorpionTest$ppi,
      cellsMetadata = bad_metadata,
      groupBy = "region",
      showProgress = FALSE
    )
  )
})

test_that("runSCORPION() validates groupBy column existence", {
  data(scorpionTest)
  
  # Should error with invalid groupBy column
  expect_error(
    runSCORPION(
      gexMatrix = scorpionTest$gex,
      tfMotifs = scorpionTest$tf,
      ppiNet = scorpionTest$ppi,
      cellsMetadata = scorpionTest$metadata,
      groupBy = "nonexistent_column",
      showProgress = FALSE
    )
  )
})

test_that("runSCORPION() batch effect correction requires batch parameter", {
  data(scorpionTest)
  
  # Should error when removeBatchEffect=TRUE but no batch parameter
  expect_error(
    runSCORPION(
      gexMatrix = scorpionTest$gex,
      tfMotifs = scorpionTest$tf,
      ppiNet = scorpionTest$ppi,
      cellsMetadata = scorpionTest$metadata,
      groupBy = "region",
      removeBatchEffect = TRUE,
      batch = NULL,
      showProgress = FALSE
    )
  )
})
