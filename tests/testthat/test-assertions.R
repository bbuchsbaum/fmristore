library(testthat)
library(fmristore)

# Tests for assert_non_empty_numeric

test_that("assert_non_empty_numeric errors on invalid input", {
  expect_error(assert_non_empty_numeric(NULL, "x", "fn"),
               "must be a non-empty numeric vector")
  expect_error(assert_non_empty_numeric("a", "x", "fn"),
               "must be a non-empty numeric vector")
  expect_error(assert_non_empty_numeric(numeric(0), "x", "fn"),
               "must be a non-empty numeric vector")
})

# Tests for validate_same_dims

test_that("validate_same_dims returns NULL when dimensions match", {
  a <- array(1, dim = c(2, 3, 4))
  b <- array(1, dim = c(2, 3, 4))
  expect_null(validate_same_dims(a, b))
})

test_that("validate_same_dims returns message when dimensions do not match", {
  a <- array(1, dim = c(2, 3, 4))
  b <- array(1, dim = c(2, 3, 5))
  msg <- validate_same_dims(a, b)
  expect_type(msg, "character")
  expect_match(msg, "Dimension mismatch")
})
