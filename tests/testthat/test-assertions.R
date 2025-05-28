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
