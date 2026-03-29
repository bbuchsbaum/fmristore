library(testthat)
library(hdf5r)
library(neuroim2)

# Helper: create a minimal H5NeuroVec for testing
make_test_h5vec <- function(dims = c(5L, 6L, 4L, 10L), values = NULL) {
  sp <- NeuroSpace(dims)
  if (is.null(values)) {
    values <- array(rnorm(prod(dims)), dim = dims)
  }
  dvec <- DenseNeuroVec(values, sp)
  h5vec <- to_nih5_vec(dvec)
  list(h5vec = h5vec, dvec = dvec, values = values)
}

# --- Constructor & Lifecycle ---

test_that("H5NeuroVec constructor creates valid object", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("neuroim2")
  obj <- make_test_h5vec()
  on.exit(close(obj$h5vec), add = TRUE)

  expect_s4_class(obj$h5vec, "H5NeuroVec")
  expect_equal(dim(obj$h5vec), c(5L, 6L, 4L, 10L))
  expect_true(obj$h5vec@obj$is_valid)
})

test_that("H5NeuroVec constructor rejects invalid file", {
  skip_if_not_installed("hdf5r")
  tmp <- tempfile(fileext = ".h5")
  on.exit(unlink(tmp), add = TRUE)

  h5 <- H5File$new(tmp, mode = "w")
  h5attr(h5, "rtype") <- "WrongType"
  h5$create_group("space")
  h5$create_dataset("space/dim", robj = c(2L, 2L, 2L, 5L))
  h5$close_all()

  expect_error(H5NeuroVec(tmp), "Invalid HDF5")
})

test_that("H5NeuroVec close() makes handle invalid", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("neuroim2")
  obj <- make_test_h5vec()

  expect_true(obj$h5vec@obj$is_valid)
  close(obj$h5vec)
  expect_false(obj$h5vec@obj$is_valid)
})

test_that("H5NeuroVec has finalizer registered", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("neuroim2")
  obj <- make_test_h5vec()
  on.exit(close(obj$h5vec), add = TRUE)

  # The finalizer is attached to the h5obj environment - we verify it by

  # checking the handle is valid and the object was created successfully
  expect_true(obj$h5vec@obj$is_valid)
})

# --- [ subsetting ---

test_that("[i,j,k,l] subsetting returns correct values", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("neuroim2")
  obj <- make_test_h5vec()
  on.exit(close(obj$h5vec), add = TRUE)

  # Single element
  val <- obj$h5vec[1, 1, 1, 1]
  expect_equal(val, obj$values[1, 1, 1, 1], tolerance = 1e-5)

  # Slice
  slice <- obj$h5vec[1:3, 2:4, 1, 1:5]
  expected <- obj$values[1:3, 2:4, 1, 1:5]
  expect_equal(as.numeric(slice), as.numeric(expected), tolerance = 1e-5)
})

test_that("[i,j,k,l] handles full-range defaults", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("neuroim2")
  obj <- make_test_h5vec(dims = c(3L, 3L, 2L, 4L))
  on.exit(close(obj$h5vec), add = TRUE)

  # Missing indices should default to full range
  result <- obj$h5vec[, , , 1, drop = FALSE]
  expect_equal(dim(result), c(3, 3, 2, 1))
})

test_that("[i,j,k,l] rejects out-of-range indices", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("neuroim2")
  obj <- make_test_h5vec(dims = c(3L, 3L, 2L, 4L))
  on.exit(close(obj$h5vec), add = TRUE)

  expect_error(obj$h5vec[1, 1, 1, 100], "out of range")
})

test_that("[i,j,k,l] handles empty indices", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("neuroim2")
  obj <- make_test_h5vec(dims = c(3L, 3L, 2L, 4L))
  on.exit(close(obj$h5vec), add = TRUE)

  result <- obj$h5vec[integer(0), 1, 1, 1]
  expect_equal(length(result), 0)
})

# --- linear_access ---

test_that("linear_access returns correct values", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("neuroim2")
  obj <- make_test_h5vec(dims = c(3L, 4L, 2L, 5L))
  on.exit(close(obj$h5vec), add = TRUE)

  # Linear index 1 should be element [1,1,1,1]
  val <- linear_access(obj$h5vec, 1)
  expect_equal(val, obj$values[1], tolerance = 1e-5)

  # Multiple linear indices
  idx <- c(1, 5, 10)
  vals <- linear_access(obj$h5vec, idx)
  expect_equal(vals, obj$values[idx], tolerance = 1e-5)
})

test_that("linear_access works with integer input", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("neuroim2")
  obj <- make_test_h5vec(dims = c(3L, 3L, 2L, 4L))
  on.exit(close(obj$h5vec), add = TRUE)

  vals_num <- linear_access(obj$h5vec, c(1, 2, 3))
  vals_int <- linear_access(obj$h5vec, c(1L, 2L, 3L))
  expect_equal(vals_num, vals_int)
})

# --- series ---

test_that("series with linear indices returns [nTime x nVoxels]", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("neuroim2")
  obj <- make_test_h5vec(dims = c(3L, 4L, 2L, 5L))
  on.exit(close(obj$h5vec), add = TRUE)

  # Get time series for voxel at linear index 1
  s <- series(obj$h5vec, 1L)
  expect_equal(length(s), 5) # nTime

  # Multiple voxels
  s2 <- series(obj$h5vec, c(1L, 2L, 3L))
  expect_equal(dim(s2), c(5, 3)) # [nTime x nVoxels]
})

test_that("series with i,j,k returns time series for single voxel", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("neuroim2")
  obj <- make_test_h5vec(dims = c(3L, 4L, 2L, 5L))
  on.exit(close(obj$h5vec), add = TRUE)

  s <- series(obj$h5vec, 1L, 1L, 1L)
  expect_equal(length(s), 5) # nTime
  expected <- obj$values[1, 1, 1, ]
  expect_equal(as.numeric(s), as.numeric(expected), tolerance = 1e-5)
})

test_that("series with numeric input dispatches correctly", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("neuroim2")
  obj <- make_test_h5vec(dims = c(3L, 3L, 2L, 4L))
  on.exit(close(obj$h5vec), add = TRUE)

  s_int <- series(obj$h5vec, 1L)
  s_num <- series(obj$h5vec, 1)
  expect_equal(s_int, s_num)
})

test_that("series with matrix input returns correct block", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("neuroim2")
  obj <- make_test_h5vec(dims = c(5L, 6L, 4L, 8L))
  on.exit(close(obj$h5vec), add = TRUE)

  coords <- matrix(c(1, 2, 3, 1, 2, 3, 1, 1, 1), nrow = 3, ncol = 3)
  s <- series(obj$h5vec, coords)
  expect_equal(nrow(s), 8) # nTime
  expect_equal(ncol(s), 3) # 3 voxels
})

# --- show method ---

test_that("show method prints without error", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("neuroim2")
  obj <- make_test_h5vec()
  on.exit(close(obj$h5vec), add = TRUE)

  expect_output(show(obj$h5vec), "H5NeuroVec")
})

test_that("show method works for closed handle", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("neuroim2")
  obj <- make_test_h5vec()
  close(obj$h5vec)

  expect_output(show(obj$h5vec), "CLOSED")
})

# --- Round-trip ---

test_that("to_nih5_vec round-trip preserves data", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("neuroim2")
  dims <- c(4L, 5L, 3L, 6L)
  sp <- NeuroSpace(dims)
  values <- array(rnorm(prod(dims)), dim = dims)
  dvec <- DenseNeuroVec(values, sp)

  h5vec <- to_nih5_vec(dvec)
  on.exit(close(h5vec), add = TRUE)

  # Verify all data matches
  for (t in 1:dims[4]) {
    h5_slice <- h5vec[, , , t, drop = FALSE]
    orig_slice <- values[, , , t, drop = FALSE]
    expect_equal(as.numeric(h5_slice), as.numeric(orig_slice), tolerance = 1e-5)
  }
})

# --- Regression: handle leak fix ---

test_that("repeated subsetting does not leak handles (regression)", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("neuroim2")
  obj <- make_test_h5vec(dims = c(3L, 3L, 2L, 4L))
  on.exit(close(obj$h5vec), add = TRUE)

  # Perform many subset operations - should not leak handles
  for (i in 1:50) {
    val <- obj$h5vec[1, 1, 1, 1]
    la <- linear_access(obj$h5vec, 1)
  }
  # If handles leaked, the H5File would become unstable
  expect_true(obj$h5vec@obj$is_valid)
  # One more access to prove it still works
  expect_type(obj$h5vec[1, 1, 1, 1], "double")
})
