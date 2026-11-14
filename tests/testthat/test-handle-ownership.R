# HDF5 shares one file structure between all handles opened on the same path.
# Closing a handle must therefore close only that handle; closing "all objects
# in the file" also invalidates handles owned by unrelated objects.

test_that("close_h5_safely closes an H5File", {
  skip_if_not_installed("hdf5r")
  tmp <- tempfile(fileext = ".h5")
  on.exit(unlink(tmp), add = TRUE)

  h5 <- hdf5r::H5File$new(tmp, mode = "w")
  h5$create_group("test_group")
  expect_true(h5$is_valid)

  fmristore:::close_h5_safely(h5)
  expect_false(h5$is_valid)
})

test_that("close_h5_safely closes dataset with close", {
  skip_if_not_installed("hdf5r")
  tmp <- tempfile(fileext = ".h5")
  on.exit(unlink(tmp), add = TRUE)

  h5 <- hdf5r::H5File$new(tmp, mode = "w")
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

  h5 <- hdf5r::H5File$new(tmp, mode = "w")
  h5$close_all()
  expect_false(h5$is_valid)

  # Should not error on double-close
  expect_silent(fmristore:::close_h5_safely(h5))
})

test_that("closing a file handle releases leaked child handles", {
  path <- tempfile(fileext = ".h5")
  on.exit(unlink(path), add = TRUE)
  h5 <- hdf5r::H5File$new(path, mode = "w")
  h5$create_dataset("x", robj = 1:10)
  close_h5_safely(h5)
  # Truncation fails while any handle on the file is still open.
  reopened <- hdf5r::H5File$new(path, mode = "w")
  expect_true(reopened$is_valid)
  reopened$close_all()
})

test_that("an open HDF5 frame source session survives unrelated reads of its file", {
  observations <- fmridataset::axis_frame(data.frame(.obs_id = paste0("obs-", 1:4)))
  space <- fmridataset::volume_space(c(2L, 2L, 1L), support = 1:4, template = "ownership")
  beta <- matrix(as.double(seq_len(16)), nrow = 4L)
  frame <- fmridataset::fmri_frame(
    assays = list(beta = beta, variance = beta / 10),
    observations = observations,
    features = fmridataset::feature_axis(
      data.frame(.feature_id = fmridataset::feature_ids(space)), space = space
    ),
    active_assay = "beta"
  )
  path <- tempfile(fileext = ".h5")
  on.exit(unlink(path), add = TRUE)
  write_frame_h5(frame, path)
  reopened <- open_frame_h5(path)

  session <- fmridataset::source_open(fmridataset::assay(reopened, "beta")$source)
  on.exit(fmridataset::source_close(session), add = TRUE, after = FALSE)
  invisible(fmridataset::collect_assay(reopened, "variance"))
  invisible(frame_store_info(path))
  expect_equal(fmridataset::source_read(session, 2L, 3L), beta[2, 3, drop = FALSE])
})
