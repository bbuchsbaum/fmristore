library(testthat)
library(hdf5r)
library(neuroim2)

# Regression tests for H5NeuroVol, H5NeuroVec and LatentNeuroVec HDF5 I/O.
# hdf5r also exports an h5file() function, so fmristore's generic is called
# with its namespace below.

# --- Helpers ---

make_regression_lv <- function() {
  set.seed(2)
  mask_arr <- array(FALSE, c(4, 4, 3))
  mask_arr[2:3, 2:4, 1:2] <- TRUE
  sp4 <- NeuroSpace(c(4, 4, 3, 6), spacing = c(2, 2, 2), origin = c(-3, -3, -2))
  mask <- LogicalNeuroVol(mask_arr, drop_dim(sp4))
  nv <- sum(mask_arr)
  basis <- matrix(rnorm(6 * 2), 6, 2)
  loadings <- matrix(rnorm(nv * 2), nv, 2)
  lv <- suppressWarnings(
    LatentNeuroVec(basis, loadings, sp4, mask, offset = rnorm(nv), label = "runA")
  )
  list(lv = lv, sp4 = sp4, loadings = loadings)
}

write_regression_lv <- function(lv) {
  base <- tempfile("latent_")
  suppressWarnings(write_vec(lv, base))
  paste0(base, ".lv.h5")
}

load_regression_lv <- function(file) {
  suppressMessages(load_data(fmristore:::LatentNeuroVecSource(file)))
}

# --- LatentNeuroVec ---

test_that("write_vec/load_data round-trip keeps the LatentNeuroVec origin", {
  fx <- make_regression_lv()
  f <- write_regression_lv(fx$lv)
  on.exit(unlink(f), add = TRUE)

  lv2 <- load_regression_lv(f)
  sp2 <- space(lv2)
  expect_equal(origin(sp2), origin(fx$sp4))
  expect_equal(spacing(sp2), spacing(fx$sp4))
  expect_equal(trans(sp2), trans(fx$sp4))
  expect_equal(as.array(lv2), as.array(fx$lv),
    tolerance = 1e-5, ignore_attr = TRUE
  )

  # A flipped, rotated transform survives the quaternion round-trip too
  th <- 0.3
  tm <- diag(c(-2, 2, 3, 1))
  tm[1:3, 1:3] <- rbind(c(cos(th), -sin(th), 0), c(sin(th), cos(th), 0), c(0, 0, 1)) %*%
    tm[1:3, 1:3]
  tm[1:3, 4] <- c(90, -126, -72)
  sp_obl <- NeuroSpace(c(4, 4, 3, 6), trans = tm)
  lv_obl <- suppressWarnings(LatentNeuroVec(
    as.matrix(fx$lv@basis), fx$loadings, sp_obl,
    LogicalNeuroVol(fx$lv@mask@.Data, drop_dim(sp_obl))
  ))
  f_obl <- write_regression_lv(lv_obl)
  on.exit(unlink(f_obl), add = TRUE)
  sp_obl2 <- space(load_regression_lv(f_obl))
  expect_equal(trans(sp_obl2), tm)
  expect_equal(origin(sp_obl2), c(90, -126, -72))
})

test_that("load_data and validate_latent_file accept a single-slice (Z = 1) mask", {
  mask_arr <- array(TRUE, c(4, 4, 1))
  sp4 <- NeuroSpace(c(4, 4, 1, 5))
  mask <- LogicalNeuroVol(mask_arr, drop_dim(sp4))
  set.seed(3)
  lv <- suppressWarnings(LatentNeuroVec(
    matrix(rnorm(10), 5, 2), matrix(rnorm(32), 16, 2), sp4, mask
  ))
  f <- write_regression_lv(lv)
  on.exit(unlink(f), add = TRUE)

  expect_true(validate_latent_file(f))
  lv2 <- load_regression_lv(f)
  expect_equal(dim(lv2), c(4, 4, 1, 5))
  expect_equal(as.array(lv2), as.array(lv),
    tolerance = 1e-5, ignore_attr = TRUE
  )
})

test_that("load_data accepts a mask with an interior singleton dimension", {
  mask_arr <- array(c(TRUE, FALSE, TRUE, TRUE), c(4, 1, 3))
  sp4 <- NeuroSpace(c(4, 1, 3, 5))
  mask <- LogicalNeuroVol(mask_arr, drop_dim(sp4))
  nv <- sum(mask_arr)
  set.seed(4)
  lv <- suppressWarnings(LatentNeuroVec(
    matrix(rnorm(10), 5, 2), matrix(rnorm(nv * 2), nv, 2), sp4, mask
  ))
  f <- write_regression_lv(lv)
  on.exit(unlink(f), add = TRUE)

  expect_true(validate_latent_file(f))
  lv2 <- load_regression_lv(f)
  expect_equal(dim(lv2), c(4, 1, 3, 5))
  expect_equal(lv2@mask@.Data, mask_arr, ignore_attr = TRUE)
})

test_that("load_data reads a CSR-stored sparse basis", {
  fx <- make_regression_lv()
  f <- write_regression_lv(fx$lv)
  on.exit(unlink(f), add = TRUE)

  # Replace the dense basis with a sparsified CSR (k x nVox) copy.
  loadings <- fx$loadings
  loadings[abs(loadings) < 0.5] <- 0
  h <- H5File$new(f, mode = "r+")
  h$link_delete("basis/basis_matrix")
  csr <- methods::as(Matrix::Matrix(t(loadings), sparse = TRUE), "RsparseMatrix")
  g <- h$create_group("basis/basis_matrix_sparse")
  h5attr(g, "storage") <- "csr"
  h5attr(g, "shape") <- dim(csr)
  g[["data"]] <- csr@x
  g[["indices"]] <- as.integer(csr@j)
  g[["indptr"]] <- as.integer(csr@p)
  g$close()
  h$close_all()

  lv2 <- load_regression_lv(f)
  expect_equal(dim(lv2@loadings), dim(loadings))
  expect_equal(as.matrix(lv2@loadings), loadings,
    tolerance = 1e-5, ignore_attr = TRUE
  )
  expect_equal(as.matrix(lv2@basis), as.matrix(fx$lv@basis),
    tolerance = 1e-5, ignore_attr = TRUE
  )
})

# --- H5NeuroVol ---

test_that("H5NeuroVol [i, , k] honours k and [i] indexes linearly", {
  arr <- array(as.numeric(seq_len(4 * 5 * 6)), c(4, 5, 6))
  vol <- NeuroVol(arr, NeuroSpace(c(4, 5, 6)))
  h <- as_h5(vol, tempfile(fileext = ".h5"), data_type = "DOUBLE")
  on.exit(close(h), add = TRUE)

  # neuroim2's DenseNeuroVol is the reference behaviour
  expect_equal(h[1:2, , 3], vol[1:2, , 3])
  expect_equal(h[2, , 4], vol[2, , 4])
  expect_equal(h[2, , 4, drop = FALSE], vol[2, , 4, drop = FALSE])
  expect_equal(h[2L, , 4L], vol[2L, , 4L])
  expect_equal(h[c(3, 1), , c(6, 2)], vol[c(3, 1), , c(6, 2)])
  expect_equal(h[1:2, , ], vol[1:2, , ])
  expect_equal(h[1, , , drop = FALSE], vol[1, , , drop = FALSE])

  expect_equal(h[2], vol[2])
  expect_equal(h[1:3], vol[1:3])
  expect_equal(h[c(2, 50)], vol[c(2, 50)])
  expect_equal(h[3L], vol[3L])

  expect_error(h[2, , 7], "out of range")
  expect_error(h[length(arr) + 1], "out of range")
})

test_that("h5file() and write_dataset() work for H5NeuroVol", {
  arr <- array(as.numeric(1:24), c(2, 3, 4))
  vol <- NeuroVol(arr, NeuroSpace(c(2, 3, 4)))
  src <- tempfile(fileext = ".h5")
  out <- tempfile(fileext = ".h5")
  on.exit(unlink(c(src, out)), add = TRUE)

  hv <- as_h5(vol, src, data_type = "DOUBLE")
  on.exit(close(hv), add = TRUE, after = FALSE)

  expect_equal(normalizePath(fmristore::h5file(hv)), normalizePath(src))
  expect_message(write_dataset(hv, out), "Copying H5NeuroVol")
  expect_true(file.exists(out))

  copy <- H5NeuroVol(out)
  on.exit(close(copy), add = TRUE, after = FALSE)
  expect_equal(copy[1:2, 1:3, 1:4], arr)
})

# --- H5NeuroVec ---

make_regression_h5vec <- function(dims = c(4, 4, 3, 5)) {
  set.seed(1)
  arr <- array(rnorm(prod(dims)), dims)
  dvec <- DenseNeuroVec(arr, NeuroSpace(dims))
  h <- to_nih5_vec(dvec, tempfile(fileext = ".h5"), data_type = "DOUBLE")
  list(h = h, dvec = dvec, arr = arr)
}

test_that("series(H5NeuroVec, matrix) follows the caller's row order", {
  fx <- make_regression_h5vec()
  on.exit(close(fx$h), add = TRUE)

  # nrow equals the bounding-box voxel count, which used to skip reordering
  reversed <- rbind(c(2, 1, 1), c(1, 1, 1))
  expect_equal(series(fx$h, reversed), series(fx$dvec, reversed))

  dupes <- rbind(c(1, 1, 1), c(1, 1, 1), c(2, 2, 1), c(1, 2, 1))
  expect_equal(series(fx$h, dupes), series(fx$dvec, dupes))

  shuffled <- rbind(c(3, 2, 2), c(1, 4, 1), c(2, 3, 3), c(4, 1, 2))
  expect_equal(series(fx$h, shuffled), series(fx$dvec, shuffled))

  one <- rbind(c(3, 2, 2))
  expect_equal(series(fx$h, one), series(fx$dvec, one))
})

test_that("linear_access(H5NeuroVec) rejects out-of-range indices", {
  fx <- make_regression_h5vec()
  on.exit(close(fx$h), add = TRUE)
  tot <- length(fx$arr)

  expect_error(linear_access(fx$h, 0), "out of range")
  expect_error(linear_access(fx$h, tot + 1), "out of range")
  expect_error(linear_access(fx$h, -1), "out of range")
  expect_error(linear_access(fx$h, c(1, tot + 5)), "out of range")
  expect_identical(linear_access(fx$h, integer(0)), numeric(0))
  expect_identical(linear_access(fx$h, numeric(0)), numeric(0))
  expect_equal(linear_access(fx$h, c(tot, 1, 7)), fx$arr[c(tot, 1, 7)])
})

test_that("H5NeuroVec() opens files written by to_nih5_vec() by default", {
  dims <- c(3, 3, 2, 4)
  set.seed(5)
  arr <- array(rnorm(prod(dims)), dims)
  dvec <- DenseNeuroVec(arr, NeuroSpace(dims))
  f <- tempfile(fileext = ".h5")
  on.exit(unlink(f), add = TRUE)
  close(to_nih5_vec(dvec, f, data_type = "DOUBLE"))

  h <- H5NeuroVec(f)
  on.exit(close(h), add = TRUE, after = FALSE)
  expect_equal(h@dataset_name, "data/elements")
  expect_equal(h[1, 1, 1, 1], arr[1, 1, 1, 1])
  expect_equal(h[, , , 3], arr[, , , 3])
})

test_that("H5NeuroVec() still reads a 4-D dataset stored directly at /data", {
  dims <- c(3L, 2L, 2L, 3L)
  arr <- array(as.numeric(seq_len(prod(dims))), dims)
  f <- tempfile(fileext = ".h5")
  on.exit(unlink(f), add = TRUE)

  h5 <- H5File$new(f, mode = "w")
  h5attr(h5, "rtype") <- "DenseNeuroVec"
  h5$create_group("space")
  h5[["space/dim"]] <- dims
  h5[["space/origin"]] <- c(0, 0, 0)
  h5[["space/trans"]] <- diag(4)
  h5[["data"]] <- arr
  h5$close_all()

  h <- H5NeuroVec(f)
  on.exit(close(h), add = TRUE, after = FALSE)
  expect_equal(h@dataset_name, "data")
  expect_equal(h[2, 1, 2, 3], arr[2, 1, 2, 3])
})

test_that("h5file() and write_dataset() work for H5NeuroVec", {
  dims <- c(3, 3, 2, 4)
  set.seed(6)
  arr <- array(rnorm(prod(dims)), dims)
  dvec <- DenseNeuroVec(arr, NeuroSpace(dims))
  src <- tempfile(fileext = ".h5")
  out <- tempfile(fileext = ".h5")
  on.exit(unlink(c(src, out)), add = TRUE)
  close(to_nih5_vec(dvec, src, data_type = "DOUBLE"))

  h <- H5NeuroVec(src)
  on.exit(close(h), add = TRUE, after = FALSE)

  expect_equal(normalizePath(fmristore::h5file(h)), normalizePath(src))
  expect_message(write_dataset(h, out), "Copying H5NeuroVec")
  expect_true(file.exists(out))

  copy <- H5NeuroVec(out)
  on.exit(close(copy), add = TRUE, after = FALSE)
  expect_equal(copy[, , , ], arr)
})

test_that("read_dataset() detects and reads LatentNeuroVec files", {
  fx <- make_regression_lv()
  f <- write_regression_lv(fx$lv)
  on.exit(unlink(f), add = TRUE)

  expect_identical(detect_h5_type(f), "LatentNeuroVec")
  auto <- suppressMessages(read_dataset(f))
  typed <- suppressMessages(read_dataset(f, type = "LatentNeuroVec"))
  for (lv in list(auto, typed)) {
    expect_s4_class(lv, "LatentNeuroVec")
    expect_equal(dim(lv), dim(fx$lv))
    expect_equal(origin(space(lv)), origin(fx$sp4))
    expect_equal(as.array(lv), as.array(fx$lv), tolerance = 1e-5, ignore_attr = TRUE)
  }
})

test_that("H5NeuroVol and H5NeuroVec materialize like in-memory arrays", {
  # Compare with @.Data: as.array() on an in-memory NeuroVec returns the S4
  # object itself when the installed package runs after other test files.
  vol <- DenseNeuroVol(array(as.double(1:24), c(2, 3, 4)), NeuroSpace(c(2, 3, 4), origin = c(1, 2, 3)))
  h5vol <- as_h5(vol, file = tempfile(fileext = ".h5"))
  on.exit(close(h5vol), add = TRUE)
  expect_equal(h5vol[, , ], vol@.Data, ignore_attr = TRUE)
  expect_identical(dim(h5vol[, , ]), c(2L, 3L, 4L))
  expect_identical(dim(h5vol[, , 2]), c(2L, 3L))
  expect_equal(as.array(h5vol), vol@.Data, ignore_attr = TRUE)
  dense_vol <- as(h5vol, "DenseNeuroVol")
  expect_s4_class(dense_vol, "DenseNeuroVol")
  expect_equal(dense_vol@.Data, vol@.Data, ignore_attr = TRUE)
  expect_equal(origin(space(dense_vol)), origin(space(vol)))

  vec <- DenseNeuroVec(array(as.double(1:48), c(2, 3, 4, 2)), NeuroSpace(c(2, 3, 4, 2)))
  h5vec <- as_h5(vec, file = tempfile(fileext = ".h5"))
  on.exit(close(h5vec), add = TRUE)
  expect_equal(as.array(h5vec), vec@.Data, ignore_attr = TRUE)
  dense_vec <- as(h5vec, "DenseNeuroVec")
  expect_s4_class(dense_vec, "DenseNeuroVec")
  expect_equal(dense_vec@.Data, vec@.Data, ignore_attr = TRUE)
})
