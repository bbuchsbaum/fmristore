library(testthat)
library(hdf5r)

# Tests for utility helpers in h5_utils.R that currently lack coverage

# -------------------------------------------------------------------------
# open_h5
# -------------------------------------------------------------------------

test_that("open_h5 handles file paths and existing handles", {
  skip_if_not_installed("hdf5r")
  tmp <- tempfile(fileext = ".h5")
  h5 <- H5File$new(tmp, mode = "w")
  h5$create_group("grp")
  h5$close_all()

  res1 <- fmristore:::open_h5(tmp, mode = "r")
  expect_true(inherits(res1$h5, "H5File"))
  expect_true(res1$owns)
  expect_true(res1$h5$is_valid)
  res1$h5$close_all()

  h5b <- H5File$new(tmp, mode = "r")
  res2 <- fmristore:::open_h5(h5b, mode = "r")
  expect_identical(res2$h5, h5b)
  expect_false(res2$owns)
  h5b$close_all()
  unlink(tmp)
})

# -------------------------------------------------------------------------
# h5_dataset_dims / h5_read / h5_read_subset
# -------------------------------------------------------------------------

test_that("h5_read and friends retrieve data correctly", {
  skip_if_not_installed("hdf5r")
  tmp <- tempfile(fileext = ".h5")
  h5w <- H5File$new(tmp, mode = "w")
  mat <- matrix(1:12, nrow = 3, ncol = 4)
  h5w$create_dataset("mat", robj = mat)
  h5w$close_all()

  h5 <- H5File$new(tmp, mode = "r")
  expect_equal(fmristore:::h5_dataset_dims(h5, "mat"), c(3L, 4L))
  expect_equal(fmristore:::h5_read(h5, "mat"), mat)
  sub <- fmristore:::h5_read_subset(h5, "mat", index = list(2:3, 3:4))
  expect_equal(sub, mat[2:3, 3:4, drop = FALSE])
  h5$close_all()
  unlink(tmp)
})

# -------------------------------------------------------------------------
# guess_h5_type / map_dtype
# -------------------------------------------------------------------------

test_that("guess_h5_type and map_dtype cover common types", {
  skip_if_not_installed("hdf5r")
  t_int <- fmristore:::guess_h5_type(1L)
  t_dbl <- fmristore:::guess_h5_type(1.0)
  t_lgl <- fmristore:::guess_h5_type(TRUE)

  expect_equal(t_int$get_class(), hdf5r::h5const$H5T_INTEGER)
  expect_equal(t_dbl$get_class(), hdf5r::h5const$H5T_FLOAT)
  expect_equal(t_lgl$get_size(), 1L)

  expect_equal(fmristore:::map_dtype(hdf5r::h5types$H5T_NATIVE_INT32), c(8L, 32L))
  expect_equal(fmristore:::map_dtype(hdf5r::h5types$H5T_NATIVE_UINT16), c(512L, 16L))
  expect_equal(fmristore:::map_dtype(hdf5r::h5types$H5T_IEEE_F32LE), c(16L, 32L))
})
