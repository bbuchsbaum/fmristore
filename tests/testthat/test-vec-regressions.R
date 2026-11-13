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
    LogicalNeuroVol(as.array(fx$lv@mask), drop_dim(sp_obl))
  ))
  f_obl <- write_regression_lv(lv_obl)
  on.exit(unlink(f_obl), add = TRUE)
  sp_obl2 <- space(load_regression_lv(f_obl))
  expect_equal(trans(sp_obl2), tm)
  expect_equal(origin(sp_obl2), c(90, -126, -72))
})
