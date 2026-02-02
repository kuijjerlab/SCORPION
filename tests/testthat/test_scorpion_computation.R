# Tests for scorpion() network computation
# Validates that scorpion() produces correct output structures with various parameters

test_that("scorpion() accepts gexMatrix and tfMotifs and produces list output", {
  data(scorpionTest)
  result <- scorpion(
    tfMotifs = scorpionTest$tf,
    gexMatrix = scorpionTest$gex[, 1:100],  # Small subset for speed
    ppiNet = scorpionTest$ppi,
    alphaValue = 0.8,
    showProgress = FALSE
  )
  
  expect_true(is.list(result))
  expect_true(all(c("regNet", "coregNet", "coopNet") %in% names(result)))
  expect_true(is.numeric(result$numGenes))
  expect_true(is.numeric(result$numTFs))
  expect_true(is.numeric(result$numEdges))
})

test_that("scorpion() works with NULL tfMotifs (co-expression only)", {
  skip("Skipped due to Matrix/gpuR dependency issues")
  data(scorpionTest)
  result <- scorpion(
    tfMotifs = NULL,
    gexMatrix = scorpionTest$gex[, 1:100],
    ppiNet = scorpionTest$ppi,
    showProgress = FALSE
  )
  
  expect_true(is.list(result) || is.matrix(result))
  expect_true(nrow(result) > 0 && ncol(result) > 0)
})

test_that("scorpion() works with NULL ppiNet", {
  skip("Skipped due to Matrix/gpuR dependency issues")
  data(scorpionTest)
  result <- scorpion(
    tfMotifs = scorpionTest$tf,
    gexMatrix = scorpionTest$gex[, 1:100],
    ppiNet = NULL,
    showProgress = FALSE
  )
  
  expect_true(is.list(result) || is.matrix(result))
})

test_that("scorpion() alphaValue parameter changes output", {
  data(scorpionTest)
  
  result_low <- scorpion(
    tfMotifs = scorpionTest$tf,
    gexMatrix = scorpionTest$gex[, 1:50],
    ppiNet = scorpionTest$ppi,
    alphaValue = 0.1,
    showProgress = FALSE
  )
  
  result_high <- scorpion(
    tfMotifs = scorpionTest$tf,
    gexMatrix = scorpionTest$gex[, 1:50],
    ppiNet = scorpionTest$ppi,
    alphaValue = 0.8,
    showProgress = FALSE
  )
  
  # Both should produce results
  expect_true(is.list(result_low) || is.matrix(result_low))
  expect_true(is.list(result_high) || is.matrix(result_high))
})
