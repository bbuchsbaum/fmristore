# HDF5 shares one file structure between all handles opened on the same path.
# Closing a handle must therefore close only that handle; closing "all objects
# in the file" also invalidates handles owned by unrelated objects.

library(neuroim2)

.ownership_cluster_fixture <- function() {
  dims <- c(4L, 3L, 2L)
  arr <- array(FALSE, dims)
  arr[c(1, 2, 4, 6, 7, 9)] <- TRUE
  mask <- LogicalNeuroVol(arr, NeuroSpace(dims))
  cl <- c(1L, 2L, 1L, 2L, 2L, 1L)
  cvol <- ClusteredNeuroVol(mask, cl)
  full <- outer(seq_along(cl) * 100, seq_len(3L), "+")
  run <- function(name) {
    list(
      scan_name = name, type = "full",
      data = list(
        cluster_1 = full[cl == 1L, , drop = FALSE],
        cluster_2 = full[cl == 2L, , drop = FALSE]
      ),
      metadata = list()
    )
  }
  list(mask = mask, cvol = cvol, run = run)
}

test_that("opening and closing a path does not invalidate other objects on it", {
  path <- tempfile(fileext = ".h5")
  on.exit(unlink(path), add = TRUE)
  vec <- DenseNeuroVec(array(seq_len(4 * 4 * 3 * 2), c(4, 4, 3, 2)), NeuroSpace(c(4, 4, 3, 2)))
  live <- as_h5(vec, file = path)
  on.exit(close(live), add = TRUE, after = FALSE)

  fh <- open_h5(path, mode = "r")
  close_h5_safely(fh$h5)
  expect_true(live@obj$is_valid)

  expect_identical(detect_h5_type(path), "H5NeuroVec")
  expect_true(live@obj$is_valid)
  expect_equal(live[1, 1, 1, 2], vec[1, 1, 1, 2])
})

test_that("closing one parcellated scan leaves another on the same file usable", {
  fx <- .ownership_cluster_fixture()
  path <- tempfile(fileext = ".h5")
  on.exit(unlink(path), add = TRUE)
  write_parcellated_experiment_h5(
    path, fx$mask, fx$cvol, list(fx$run("a"), fx$run("b")),
    overwrite = TRUE, verbose = FALSE
  )
  a <- H5ParcellatedScan(path, "a", fx$mask, fx$cvol)
  b <- H5ParcellatedScan(path, "b", fx$mask, fx$cvol)
  close(b)
  expect_true(a@obj$is_valid)
  expect_equal(as.numeric(series(a, 1L)), c(101, 102, 103))
  close(a)
})

test_that("closing a file handle releases leaked child handles", {
  path <- tempfile(fileext = ".h5")
  on.exit(unlink(path), add = TRUE)
  h5 <- hdf5r::H5File$new(path, mode = "w")
  h5_write(h5, "/x", 1:10)
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

test_that("LatentNeuroVec write failures raise and remove the partial file", {
  lvec <- fmristore:::create_minimal_LatentNeuroVec(space_dims = c(3L, 3L, 2L), n_time = 4L, n_comp = 2L)
  lvec@offset <- rnorm(sum(lvec@mask) + 5) # wrong length; bypasses validity
  base <- tempfile()
  expect_error(write_vec(lvec, base), "to_h5_latentvec")
  expect_false(file.exists(paste0(base, ".lv.h5")))
})

test_that("as_h5() on a LatentNeuroVec returns an open handle the caller owns", {
  lvec <- fmristore:::create_minimal_LatentNeuroVec(space_dims = c(3L, 3L, 2L), n_time = 4L, n_comp = 2L)
  path <- tempfile(fileext = ".lv.h5")
  on.exit(unlink(path), add = TRUE)
  h5 <- as_h5(lvec, file = path)
  expect_true(h5$is_valid)
  close_h5_safely(h5)
  expect_true(validate_latent_file(path))
  # validate_latent_file() must not leave the file open
  reopened <- hdf5r::H5File$new(path, mode = "w")
  expect_true(reopened$is_valid)
  reopened$close_all()
})

test_that("load_data() on a LatentNeuroVec file leaves the file closed", {
  lvec <- suppressWarnings(fmristore:::create_minimal_LatentNeuroVec(
    space_dims = c(3L, 3L, 2L), n_time = 4L, n_comp = 2L
  ))
  base <- tempfile()
  write_vec(lvec, base)
  path <- paste0(base, ".lv.h5")
  on.exit(unlink(path), add = TRUE)
  loaded <- suppressMessages(load_data(fmristore:::LatentNeuroVecSource(path)))
  expect_s4_class(loaded, "LatentNeuroVec")
  # Truncation fails while any handle on the file is still open.
  reopened <- hdf5r::H5File$new(path, mode = "w")
  expect_true(reopened$is_valid)
  reopened$close_all()
})
