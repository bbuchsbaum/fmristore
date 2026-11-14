library(testthat)
library(hdf5r)

# --- close_h5_safely consolidation tests ---

test_that("close_h5_safely closes H5File with close_all", {
  skip_if_not_installed("hdf5r")
  tmp <- tempfile(fileext = ".h5")
  on.exit(unlink(tmp), add = TRUE)

  h5 <- H5File$new(tmp, mode = "w")
  h5$create_group("test_group")
  expect_true(h5$is_valid)

  fmristore:::close_h5_safely(h5)
  expect_false(h5$is_valid)
})

test_that("close_h5_safely closes dataset with close", {
  skip_if_not_installed("hdf5r")
  tmp <- tempfile(fileext = ".h5")
  on.exit(unlink(tmp), add = TRUE)

  h5 <- H5File$new(tmp, mode = "w")
  h5$create_dataset("data", robj = 1:10)
  dset <- h5[["data"]]
  expect_true(dset$is_valid)

  fmristore:::close_h5_safely(dset)
  expect_false(dset$is_valid)

  # File should still be open
  expect_true(h5$is_valid)
  h5$close_all()
})

test_that("close_h5_safely handles NULL gracefully", {
  expect_silent(fmristore:::close_h5_safely(NULL))
})

test_that("close_h5_safely handles already-closed object", {
  skip_if_not_installed("hdf5r")
  tmp <- tempfile(fileext = ".h5")
  on.exit(unlink(tmp), add = TRUE)

  h5 <- H5File$new(tmp, mode = "w")
  h5$close_all()
  expect_false(h5$is_valid)

  # Should not error on double-close
  expect_silent(fmristore:::close_h5_safely(h5))
})

test_that("safe_h5_close delegates to close_h5_safely", {
  skip_if_not_installed("hdf5r")
  tmp <- tempfile(fileext = ".h5")
  on.exit(unlink(tmp), add = TRUE)

  h5 <- H5File$new(tmp, mode = "w")
  expect_true(h5$is_valid)

  fmristore:::safe_h5_close(h5)
  expect_false(h5$is_valid)
})

# --- build_nifti_header tests ---

test_that("build_nifti_header creates correct default fields", {
  skip_if_not_installed("hdf5r")
  quat <- list(quaternion = c(0, 0, 0), qoffset = c(1, 2, 3), qfac = 1)
  tmat <- diag(4)

  hdr <- fmristore:::build_nifti_header(
    dims = c(10L, 10L, 5L, 20L),
    spacing = c(2.0, 2.0, 3.0),
    quat = quat,
    tmat = tmat,
    datatype_code = 16L,
    bitpix = 32L,
    descrip = "test header"
  )

  expect_type(hdr, "list")
  expect_equal(hdr$sizeof_hdr, 348L)
  expect_equal(hdr$dim, c(4L, 10L, 10L, 5L, 20L, 1L, 1L, 1L))
  expect_equal(hdr$datatype, 16L)
  expect_equal(hdr$bitpix, 32L)
  expect_equal(hdr$quatern_b, 0)
  expect_equal(hdr$qoffset_x, 1)
  expect_equal(hdr$qoffset_y, 2)
  expect_equal(hdr$qoffset_z, 3)
  expect_equal(hdr$srow_x, tmat[1, ])
  expect_equal(hdr$descrip, "test header")
  expect_equal(hdr$magic, "n+1")
})

test_that("build_nifti_header applies overrides", {
  quat <- list(quaternion = c(0, 0, 0), qoffset = c(0, 0, 0), qfac = 1)
  tmat <- diag(4)

  hdr <- fmristore:::build_nifti_header(
    dims = c(5L, 5L, 5L, 10L),
    spacing = c(1, 1, 1),
    quat = quat,
    tmat = tmat,
    overrides = list(xyzt_units = 10L, slice_end = 4L)
  )

  expect_equal(hdr$xyzt_units, 10L)
  expect_equal(hdr$slice_end, 4L)
  # Default fields should still be present
  expect_equal(hdr$qform_code, 1L)
})

test_that("build_nifti_header handles NULL quaternion fields", {
  quat <- list(quaternion = NULL, qoffset = NULL, qfac = NULL)
  tmat <- diag(4)

  hdr <- fmristore:::build_nifti_header(
    dims = c(5L, 5L, 5L, 10L),
    spacing = c(1, 1, 1),
    quat = quat,
    tmat = tmat
  )

  expect_equal(hdr$quatern_b, 0)
  expect_equal(hdr$qoffset_x, 0)
  expect_equal(hdr$pixdim[1], 1.0) # qfac defaults to 1.0
})

# --- Verbose gating tests ---

