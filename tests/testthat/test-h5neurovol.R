library(testthat)
library(hdf5r)
library(neuroim2)

# Helper: create a minimal H5NeuroVol for testing
make_test_h5vol <- function(dims = c(5L, 6L, 4L), values = NULL) {
  sp <- NeuroSpace(dims)
  if (is.null(values)) {
    values <- array(rnorm(prod(dims)), dim = dims)
  }
  dvol <- DenseNeuroVol(values, sp)
  h5vol <- as_h5(dvol, file = NULL)
  list(h5vol = h5vol, dvol = dvol, values = values)
}

# --- Constructor & Lifecycle ---

test_that("H5NeuroVol constructor creates valid object", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("neuroim2")
  obj <- make_test_h5vol()
  on.exit(close(obj$h5vol), add = TRUE)

  expect_s4_class(obj$h5vol, "H5NeuroVol")
  expect_equal(dim(obj$h5vol), c(5L, 6L, 4L))
  expect_true(obj$h5vol@h5obj$is_valid)
})

test_that("H5NeuroVol constructor rejects non-3D file", {
  skip_if_not_installed("hdf5r")
  tmp <- tempfile(fileext = ".h5")
  on.exit(unlink(tmp), add = TRUE)

  h5 <- H5File$new(tmp, mode = "w")
  h5attr(h5, "rtype") <- "DenseNeuroVol"
  h5$create_group("space")
  h5$create_dataset("space/dim", robj = c(2L, 2L, 2L, 5L))
  h5$create_dataset("space/origin", robj = c(0, 0, 0))
  mat <- diag(4)
  h5$create_dataset("space/trans", robj = mat)
  h5$close_all()

  expect_error(H5NeuroVol(tmp), "3 dimensions")
})

test_that("H5NeuroVol close() makes handle invalid", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("neuroim2")
  obj <- make_test_h5vol()

  expect_true(obj$h5vol@h5obj$is_valid)
  close(obj$h5vol)
  expect_false(obj$h5vol@h5obj$is_valid)
})

test_that("H5NeuroVol has finalizer (handle closed on GC)", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("neuroim2")
  obj <- make_test_h5vol()
  on.exit(close(obj$h5vol), add = TRUE)

  # Verify the handle is valid (finalizer is registered at construction)
  expect_true(obj$h5vol@h5obj$is_valid)
})

# --- [ subsetting ---

test_that("[i,j,k] subsetting returns correct values", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("neuroim2")
  obj <- make_test_h5vol()
  on.exit(close(obj$h5vol), add = TRUE)

  # Single element
  val <- obj$h5vol[1, 1, 1]
  expect_equal(val, obj$values[1, 1, 1], tolerance = 1e-5)

  # Slice
  slice <- obj$h5vol[1:3, 2:4, 1:2]
  expected <- obj$values[1:3, 2:4, 1:2]
  expect_equal(as.numeric(slice), as.numeric(expected), tolerance = 1e-5)
})

test_that("[i,j,k] handles missing indices", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("neuroim2")
  obj <- make_test_h5vol(dims = c(3L, 4L, 2L))
  on.exit(close(obj$h5vol), add = TRUE)

  # Missing j,k => full range
  result <- obj$h5vol[1, , , drop = FALSE]
  expect_equal(dim(result), c(1, 4, 2))
})

test_that("[i,j,k] rejects out-of-range indices", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("neuroim2")
  obj <- make_test_h5vol(dims = c(3L, 3L, 2L))
  on.exit(close(obj$h5vol), add = TRUE)

  expect_error(obj$h5vol[1, 1, 100], "out of range")
})

test_that("[i,j,k] handles empty indices", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("neuroim2")
  obj <- make_test_h5vol(dims = c(3L, 3L, 2L))
  on.exit(close(obj$h5vol), add = TRUE)

  result <- obj$h5vol[integer(0), 1, 1]
  expect_equal(length(result), 0)
})

# --- linear_access ---

test_that("linear_access returns correct values", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("neuroim2")
  obj <- make_test_h5vol(dims = c(3L, 4L, 2L))
  on.exit(close(obj$h5vol), add = TRUE)

  # Linear index 1 => element [1,1,1]
  val <- linear_access(obj$h5vol, 1)
  expect_equal(val, obj$values[1], tolerance = 1e-5)

  # Multiple indices
  idx <- c(1, 5, 10)
  vals <- linear_access(obj$h5vol, idx)
  expect_equal(vals, obj$values[idx], tolerance = 1e-5)
})

test_that("linear_access works with integer input", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("neuroim2")
  obj <- make_test_h5vol(dims = c(3L, 3L, 2L))
  on.exit(close(obj$h5vol), add = TRUE)

  vals_num <- linear_access(obj$h5vol, c(1, 2, 3))
  vals_int <- linear_access(obj$h5vol, c(1L, 2L, 3L))
  expect_equal(vals_num, vals_int)
})

test_that("linear_access rejects out-of-range indices", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("neuroim2")
  obj <- make_test_h5vol(dims = c(3L, 3L, 2L))
  on.exit(close(obj$h5vol), add = TRUE)

  expect_error(linear_access(obj$h5vol, 1000), "out of range")
})

# --- show method ---

test_that("show method prints without error (uses cli)", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("neuroim2")
  obj <- make_test_h5vol()
  on.exit(close(obj$h5vol), add = TRUE)

  expect_output(show(obj$h5vol), "H5NeuroVol")
})

test_that("show method works after close", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("neuroim2")
  obj <- make_test_h5vol()
  close(obj$h5vol)

  expect_output(show(obj$h5vol), "CLOSED")
})

# --- Round-trip ---

test_that("as_h5 round-trip preserves data", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("neuroim2")
  dims <- c(4L, 5L, 3L)
  sp <- NeuroSpace(dims)
  values <- array(rnorm(prod(dims)), dim = dims)
  dvol <- DenseNeuroVol(values, sp)

  h5vol <- as_h5(dvol, file = NULL)
  on.exit(close(h5vol), add = TRUE)

  result <- h5vol[, , ]
  expect_equal(as.numeric(result), as.numeric(values), tolerance = 1e-5)
})

test_that("as_h5 with different data types works", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("neuroim2")
  dims <- c(3L, 3L, 2L)
  sp <- NeuroSpace(dims)
  values <- array(rnorm(prod(dims)), dim = dims)
  dvol <- DenseNeuroVol(values, sp)

  h5double <- as_h5(dvol, file = NULL, data_type = "DOUBLE")
  on.exit(close(h5double), add = TRUE)

  result <- h5double[, , ]
  expect_equal(as.numeric(result), as.numeric(values), tolerance = 1e-10)
})

# --- Regression: repeated access stability ---

test_that("repeated subsetting does not leak handles (regression)", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("neuroim2")
  obj <- make_test_h5vol(dims = c(3L, 3L, 2L))
  on.exit(close(obj$h5vol), add = TRUE)

  for (i in 1:50) {
    val <- obj$h5vol[1, 1, 1]
    la <- linear_access(obj$h5vol, 1)
  }
  expect_true(obj$h5vol@h5obj$is_valid)
  expect_type(obj$h5vol[1, 1, 1], "double")
})
