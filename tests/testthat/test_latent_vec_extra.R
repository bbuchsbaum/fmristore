library(testthat)
library(hdf5r)
library(neuroim2)
library(Matrix)
library(fmristore)

# Helper to create a simple deterministic LatentNeuroVec
make_simple_lvec <- function() {
  dims <- c(3, 3, 2, 4)  # small volume
  sp <- NeuroSpace(dims)
  mask <- array(FALSE, dim = dims[1:3])
  mask[2:3, , ] <- TRUE
  mask_vol <- LogicalNeuroVol(mask, drop_dim(sp))

  k <- 2
  nt <- dims[4]
  nv <- sum(mask)

  basis <- Matrix(matrix(seq_len(nt * k), nrow = nt, ncol = k))
  loadings <- Matrix(matrix(seq_len(nv * k) / 10, nrow = nv, ncol = k))
  offset <- seq_len(nv)

  LatentNeuroVec(basis, loadings, sp, mask_vol, offset)
}


test_that("as.array reconstructs the full 4D data correctly", {
  lvec <- make_simple_lvec()
  arr <- as.array(lvec)
  expect_equal(dim(arr), dim(space(lvec)))

  basis_mat <- as.matrix(lvec@basis)
  loadings_mat <- as.matrix(lvec@loadings)
  offset_vec <- lvec@offset
  mask_idx <- which(as.logical(as.array(lvec@mask)))
  dims <- dim(lvec)

  expected <- array(0, dim = dims)
  for (t in seq_len(dims[4])) {
    vals <- as.vector(tcrossprod(basis_mat[t, , drop = FALSE], loadings_mat)) + offset_vec
    vol <- array(0, dim = dims[1:3])
    vol[mask_idx] <- vals
    expected[,,, t] <- vol
  }

  expect_equal(arr, expected, tolerance = 1e-12)
})


test_that("linear_access returns correct values and errors out-of-range", {
  lvec <- make_simple_lvec()
  full_arr <- as.array(lvec)
  nels <- prod(dim(full_arr))
  set.seed(123)
  idx <- sample(nels, 10)

  expect_equal(linear_access(lvec, idx), full_arr[idx], tolerance = 1e-12)
  expect_error(linear_access(lvec, 0), "positive integers")
  expect_error(linear_access(lvec, nels + 1), "exceed the total number of elements")
})

