library(testthat)
library(hdf5r)
library(neuroim2)

# Ensure invalid H5NeuroVol closes handle on error

test_that("invalid H5NeuroVol file closes handle", {
  skip_if_not_installed("hdf5r")
  tmp <- tempfile(fileext = ".h5")
  on.exit(unlink(tmp), add = TRUE)

  # create a malformed file
  h5 <- H5File$new(tmp, mode = "w")
  h5attr(h5, "rtype") <- "WrongType"
  h5$create_group("space")
  h5$create_dataset("space/dim", robj = c(2L,2L,2L))
  h5$close_all()

  before <- length(h5validObjects())
  expect_error(H5NeuroVol(tmp))
  after <- length(h5validObjects())
  expect_equal(after, before)
})
