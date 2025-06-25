library(testthat)
library(hdf5r)

# Tests for error handling helpers in h5_utils.R

test_that("assert_h5_path and h5_read handle missing datasets correctly", {
  skip_if_not_installed("hdf5r")

  tmp <- tempfile(fileext = ".h5")
  on.exit(unlink(tmp), add = TRUE)

  # Create file with a simple dataset
  h5w <- H5File$new(tmp, mode = "w")
  h5w$create_dataset("data", robj = 1:5)
  h5w$close_all()

  h5 <- H5File$new(tmp, mode = "r")

  # assert_h5_path returns invisibly when the path exists
  expect_invisible(fmristore:::assert_h5_path(h5, "/data"))

  # Should mention the missing path in the error message
  expect_error(
    fmristore:::assert_h5_path(h5, "/missing"),
    regexp = "/missing"
  )

  # h5_read with missing_ok = TRUE returns NULL
  expect_null(fmristore:::h5_read(h5, "/missing", missing_ok = TRUE))

  # With missing_ok = FALSE we get an informative error
  expect_error(
    fmristore:::h5_read(h5, "/missing", missing_ok = FALSE),
    regexp = "Dataset or group not found"
  )

  h5$close_all()
})
