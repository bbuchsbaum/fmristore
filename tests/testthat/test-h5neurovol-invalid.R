library(testthat)
library(hdf5r)
library(neuroim2)

# Ensure invalid H5NeuroVol closes handle on error

test_that("invalid H5NeuroVol file throws appropriate error", {
  skip_if_not_installed("hdf5r")
  tmp <- tempfile(fileext = ".h5")
  on.exit(unlink(tmp), add = TRUE)

  # create a malformed file
  h5 <- H5File$new(tmp, mode = "w")
  h5attr(h5, "rtype") <- "WrongType"
  h5$create_group("space")
  h5$create_dataset("space/dim", robj = c(2L,2L,2L))
  h5$close_all()

  # Test that invalid file throws an error (the main functionality we want to test)
  expect_error(H5NeuroVol(tmp), "Invalid HDF5 file")
  
  # Test that the file can still be opened normally after the error
  # (indicating proper cleanup)
  h5_test <- H5File$new(tmp, mode = "r")
  expect_true(h5_test$is_valid)
  h5_test$close_all()
})
