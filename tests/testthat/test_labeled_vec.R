# test-labeledVolumeSet.R

library(testthat)
library(hdf5r)
library(neuroim2)  # or your package providing NeuroVec, DenseNeuroVol, LabeledVolumeSet
# (Make sure that you have your LabeledVolumeSet class, read_labeled_vec, write_labeled_vec defined/loaded)

test_that("Write and read a small LabeledVolumeSet", {
  # 1) Create small 4D data in memory
  set.seed(42)
  X <- 4; Y <- 5; Z <- 3; nVol <- 2  # small shape
  arr_data <- array(rnorm(X*Y*Z*nVol), dim=c(X,Y,Z,nVol))

  # 2) Create a NeuroSpace
  spc <- NeuroSpace(dim=c(X,Y,Z,nVol), spacing=c(1,1,1))

  # 3) Create a DenseNeuroVec
  #    Or if you have your own constructor: e.g. DenseNeuroVec(arr_data, space=spc)
  vec_obj <- DenseNeuroVec(arr_data, space=spc)

  # 4) Create a 3D mask (some pattern of 1/0)
  #    We'll mark ~half as valid
  mask_array <- array(sample(c(TRUE,FALSE), X*Y*Z, replace=TRUE, prob=c(0.6,0.4)), dim=c(X,Y,Z))
  mask_vol   <- LogicalNeuroVol(mask_array, space=drop_dim(spc))

  # 5) Create label names
  labels <- c("volume1", "volume2")

  # 6) Write to a temporary HDF5 file
  tmpfile <- tempfile(fileext=".h5")
  write_labeled_vec(vec=vec_obj, mask=mask_vol, labels=labels,
                    file=tmpfile, compression=4)

  # 7) Read it back as LabeledVolumeSet
  lvs <- read_labeled_vec(tmpfile, memoise=FALSE)

  # 8) Check class
  expect_s4_class(lvs, "LabeledVolumeSet")
  expect_true(inherits(lvs, "NeuroVec"))
  expect_length(lvs@labels, 2)
  expect_equal(lvs@labels, labels)

  # 9) Check that sub-volumes match
  #    We'll load each volume and compare it to arr_data[,,, i] masked
  for (i in seq_along(labels)) {
    # Single volume via [[
    vol_i <- lvs[[i]]  # or lvs[[labels[i]]]
    expect_s4_class(vol_i, "DenseNeuroVol")
    # Compare data in non-zero mask
    mask_idx <- which(mask_array)
    orig_1d  <- arr_data[,,,i][mask_idx]
    read_1d  <- as.array(vol_i)[mask_idx]
    # They should match up to a small tolerance
    expect_equal(orig_1d, read_1d, tolerance=1e-8)
  }

  # 10) Test 4D subsetting with "["
  #    e.g. lvs[i_range, j_range, k_range, l_range]
  i_range <- 2:3
  j_range <- 1:4
  k_range <- 1:2
  l_range <- 1:2

  sub_arr <- lvs[i_range, j_range, k_range, l_range, drop=FALSE]
  # shape = [2,4,2,2]
  expect_equal(dim(sub_arr), c(length(i_range), length(j_range), length(k_range), length(l_range)))
  # We can cross-check values
  # Let original subarr = arr_data[i_range, j_range, k_range, l_range]
  # Then compare
  expected_sub <- arr_data[i_range, j_range, k_range, l_range, drop=FALSE]
  expect_equal(sub_arr, expected_sub, tolerance=1e-8)

  # 11) 1D linear_access style check
  #     Suppose we want indices c(1, X*Y*Z, X*Y + 2, etc.).
  #     We'll do a small set to confirm correctness
  # Not all tests show partial usage, but here's the gist:
  # (Requires 'linear_access' method for LabeledVolumeSet)
  if (any(grepl("linear_access", as.character(getMethods("linear_access"))))) {
    # We'll pick some random linear indices in [1..X*Y*Z*nVol]
    # total voxels = X*Y*Z, times nVol => length = X*Y*Z*nVol
    tot <- X*Y*Z*nVol
    idx_samp <- sample.int(tot, 5)
    lv <- linear_access(lvs, idx_samp)
    # Re-check by creating a big 4D array and subsetting
    big_4d <- arr_data
    # => big_4d is [X, Y, Z, nVol]
    # We compare big_4d[ idx_samp ] in column-major sense
    # We do arrayInd:
    subs <- arrayInd(idx_samp, .dim=dim(big_4d))
    revals <- numeric(length(idx_samp))
    for (r in seq_along(idx_samp)) {
      revals[r] <- big_4d[ subs[r,1], subs[r,2], subs[r,3], subs[r,4] ]
    }
    expect_equal(lv, revals, tolerance=1e-8)
  }

  # 12) Done
  # if we want to keep tmpfile for debug, we skip unlink
  unlink(tmpfile)
})
