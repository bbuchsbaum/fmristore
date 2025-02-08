# test-labeledVolumeSet.R

library(testthat)
library(hdf5r)
library(neuroim2)


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

  diff_array <- sub_arr - expected_sub
  diff_array[ !mask_array[i_range,j_range,k_range] ] <- 0
  expect_true( all(abs(diff_array) < 1e-8) )


  # 11) 1D linear_access style check
  #     Suppose we want indices c(1, X*Y*Z, X*Y + 2, etc.).
  #     We'll do a small set to confirm correctness
  # Not all tests show partial usage, but here's the gist:
  # (Requires 'linear_access' method for LabeledVolumeSet)
  if (any(grepl("linear_access", as.character(findMethods("linear_access"))))) {
    # Total number of elements in the 4D array (dimensions: X, Y, Z, nVol)
    tot <- X * Y * Z * nVol

    # Pick 5 random linear indices in the full 4D space.
    idx_samp <- sample.int(tot, 5)

    # Obtain the values using linear_access on the LabeledVolumeSet.
    lv <- linear_access(lvs, idx_samp)

    # The original array 'arr_data' is unmasked.
    # However, the LabeledVolumeSet only stores data for voxels where the mask is TRUE.
    # Thus, we construct an "expected" full 4D array by applying the mask.
    # Repeat the 3D mask for each volume:
    mask_4d <- array(rep(as.logical(mask_array), nVol), dim = c(X, Y, Z, nVol))
    expected_full <- arr_data * mask_4d

    # Convert the sampled linear indices to subscripts.
    subs <- arrayInd(idx_samp, .dim = dim(expected_full))

    # For each sampled index, get the corresponding value from expected_full.
    revals <- numeric(length(idx_samp))
    for (r in seq_along(idx_samp)) {
      revals[r] <- expected_full[ subs[r, 1], subs[r, 2], subs[r, 3], subs[r, 4] ]
    }

    # Compare the output of linear_access with the expected values.
    expect_equal(lv, revals, tolerance = 1e-8)
  }

  # 12) Done
  # if we want to keep tmpfile for debug, we skip unlink
  unlink(tmpfile)
})
