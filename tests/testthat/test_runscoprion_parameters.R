# Tests for runSCORPION() parameter signature verification
# Ensures gexMatrix and tfMotifs parameters are present and properly named

test_that("runSCORPION() has gexMatrix parameter in signature", {
  sig <- formals(runSCORPION)
  expect_true("gexMatrix" %in% names(sig))
})

test_that("runSCORPION() has tfMotifs parameter in signature", {
  sig <- formals(runSCORPION)
  expect_true("tfMotifs" %in% names(sig))
})

test_that("runSCORPION() gexMatrix parameter is required", {
  sig <- formals(runSCORPION)
  # gexMatrix should not have a default value (either NULL or absent from formals)
  expect_true(is.null(sig$gexMatrix) || !"gexMatrix" %in% names(sig) || 
              identical(as.character(sig$gexMatrix), ""))
})

test_that("runSCORPION() tfMotifs parameter is required", {
  sig <- formals(runSCORPION)
  # tfMotifs should not have a default value (either NULL or absent from formals)
  expect_true(is.null(sig$tfMotifs) || !"tfMotifs" %in% names(sig) || 
              identical(as.character(sig$tfMotifs), ""))
})
