library(testthat)
library(hdf5r)
library(neuroim2)
library(fmristore)

# Tests for check_same_dims with dimension vectors

test_that("check_same_dims works with numeric dimension vectors", {
  expect_silent(fmristore:::check_same_dims(c(10, 5, 2), c(10, 5, 2)))
})

test_that("check_same_dims fails when numeric dimension vectors differ", {
  expect_error(fmristore:::check_same_dims(c(2, 3, 4), c(2, 3, 5)), "Dimension mismatch")
})
# Tests for validate_same_dims

test_that("validate_same_dims returns NULL on matching dims", {
  res <- fmristore:::validate_same_dims(c(2, 3, 4), c(2, 3, 4))
  expect_null(res)
})

test_that("validate_same_dims reports prefix on mismatch", {
  prefix <- "custom prefix:" 
  res <- fmristore:::validate_same_dims(c(2, 3, 4), c(2, 3, 5), msg = prefix)
  expect_true(!is.null(res))
  expect_true(grepl(prefix, res))
})

# Tests for assert_non_empty_numeric

test_that("assert_non_empty_numeric accepts numeric vector", {
  expect_silent(fmristore:::assert_non_empty_numeric(1:3, arg = "x", fn = "fn"))
})

test_that("assert_non_empty_numeric errors on NULL or empty", {
  expect_error(fmristore:::assert_non_empty_numeric(NULL, arg = "x", fn = "fn"), "x")
  expect_error(fmristore:::assert_non_empty_numeric(numeric(0), arg = "x", fn = "fn"), "x")
})
