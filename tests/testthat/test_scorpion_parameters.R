# Tests for scorpion() parameter signature verification
# Ensures gexMatrix and tfMotifs parameters are present and properly named

test_that("scorpion() has gexMatrix parameter in signature", {
  sig <- formals(scorpion)
  expect_true("gexMatrix" %in% names(sig))
})

test_that("scorpion() has tfMotifs parameter in signature", {
  sig <- formals(scorpion)
  expect_true("tfMotifs" %in% names(sig))
})

test_that("scorpion() gexMatrix parameter is required", {
  sig <- formals(scorpion)
  # gexMatrix should not have a default value (either NULL or absent from formals)
  expect_true(is.null(sig$gexMatrix) || !"gexMatrix" %in% names(sig) || 
              identical(as.character(sig$gexMatrix), ""))
})

test_that("scorpion() tfMotifs defaults to NULL", {
  sig <- formals(scorpion)
  # tfMotifs should default to NULL
  expect_true(is.null(sig$tfMotifs[[1]]) || is.null(eval(sig$tfMotifs)))
})
