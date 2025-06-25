library(testthat)
library(hdf5r)
library(neuroim2)
library(fmristore)

# Tests for ensure_mask helper

test_that("ensure_mask loads mask from HDF5 when mask is NULL", {
  skip_if_not_installed("hdf5r")
  sp <- NeuroSpace(c(2, 2, 2))
  mask_arr <- array(c(TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE), dim = c(2,2,2))

  tmp <- tempfile(fileext = ".h5")
  on.exit(unlink(tmp), add = TRUE)
  h5 <- H5File$new(tmp, mode = "w")
  mask_int <- array(as.integer(mask_arr), dim = dim(mask_arr))
  h5$create_dataset("mask", robj = mask_int)

  m <- fmristore:::ensure_mask(NULL, h5, sp, path = "/mask")
  expect_s4_class(m, "LogicalNeuroVol")
  expect_equal(as.logical(as.array(m@.Data)), as.logical(mask_arr))
  h5$close_all()
})


test_that("ensure_mask detects dimension mismatch", {
  skip_if_not_installed("neuroim2")
  sp <- NeuroSpace(c(2, 2, 2))
  wrong_sp <- NeuroSpace(c(2, 2, 3))
  mask_wrong <- LogicalNeuroVol(array(TRUE, dim = c(2,2,3)), wrong_sp)

  expect_error(
    fmristore:::ensure_mask(mask_wrong, NULL, sp),
    "Mask dimensions do not match space dimensions"
  )
})


test_that("ensure_mask rejects non LogicalNeuroVol inputs", {
  sp <- NeuroSpace(c(2,2,2))
  arr <- array(1, dim = c(2,2,2))

  expect_error(
    fmristore:::ensure_mask(arr, NULL, sp),
    "Provided 'mask' must be a LogicalNeuroVol object"
  )
})
