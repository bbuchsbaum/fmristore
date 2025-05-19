# test-labeledVolumeSet.R

library(testthat)
library(hdf5r)
library(neuroim2)

# -------------------------------------------------------------------------
# helper: build a tiny 4-D NeuroVec + mask for round-trip tests
make_toy_vec <- function(nx = 4, ny = 3, nz = 2, nvol = 3) {
  arr  <- array(rnorm(nx * ny * nz * nvol), dim = c(nx, ny, nz, nvol))
  spc  <- neuroim2::NeuroSpace(c(nx, ny, nz, nvol), spacing = c(1, 1, 1))
  mask <- neuroim2::LogicalNeuroVol(array(runif(nx * ny * nz) > .4,
                                          dim = c(nx, ny, nz)),
                                    space = neuroim2::drop_dim(spc))
  vec  <- neuroim2::DenseNeuroVec(arr, spc)

  list(vec = vec, mask = mask)
}


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

  # Calculate expected result
  # Equivalent to arr_data[i_range, j_range, k_range, l_range] with mask applied
  expected_sub <- array(0, dim=c(length(i_range), length(j_range), length(k_range), length(l_range)))
  for (l_pos in seq_along(l_range)) {
    l_val <- l_range[l_pos]
    # Apply mask to this volume's data
    vol_data <- arr_data[,,,l_val] * as.numeric(mask_array)
    # Extract the requested subset
    expected_sub[,,,l_pos] <- vol_data[i_range, j_range, k_range]
  }

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


test_that("LabeledVolumeSet array subsetting with [] works correctly", {
  # 1) Create small predictable 4D data
  set.seed(123)
  X <- 3; Y <- 3; Z <- 3; nVol <- 3 # 3x3x3x3 cube
  arr_data <- array(1:(X*Y*Z*nVol), dim=c(X,Y,Z,nVol))

  # 2) Create a NeuroSpace
  spc <- NeuroSpace(dim=c(X,Y,Z,nVol), spacing=c(1,1,1))

  # 3) Create a DenseNeuroVec
  vec_obj <- DenseNeuroVec(arr_data, space=spc)

  # 4) Create a 3D mask (simple checkerboard-like pattern)
  mask_array <- array(FALSE, dim=c(X,Y,Z))
  mask_array[seq(1, X*Y*Z, by=2)] <- TRUE # Mask roughly half
  mask_vol   <- LogicalNeuroVol(mask_array, space=drop_dim(spc))

  # 5) Create label names
  labels <- paste0("vol", 1:nVol)

  # 6) Write to a temporary HDF5 file
  tmpfile <- tempfile(fileext=".h5")
  on.exit(unlink(tmpfile), add = TRUE) # Ensure cleanup
  write_labeled_vec(vec=vec_obj, mask=mask_vol, labels=labels,
                    file=tmpfile, compression=0) # No compression for simplicity

  # 7) Read it back
  lvs <- read_labeled_vec(tmpfile, memoise=FALSE)

  # 8) Define Expected Data (Original data * Mask)
  #    The LVS stores 0 where the mask is FALSE.
  mask_3d_numeric <- array(as.numeric(mask_array), dim=c(X,Y,Z))
  expected_masked_data <- array(0, dim=dim(arr_data))
  for(v in 1:nVol) {
    expected_masked_data[,,,v] <- arr_data[,,,v] * mask_3d_numeric
  }

  # 9) Perform Subsetting Tests

  # Case 1: Single element access
  expect_equal(lvs[1, 1, 1, 1, drop=FALSE], expected_masked_data[1, 1, 1, 1, drop=FALSE])
  expect_equal(lvs[1, 1, 1, 1, drop=TRUE],  expected_masked_data[1, 1, 1, 1]) # Drop to scalar
  expect_equal(lvs[2, 2, 2, 2, drop=FALSE], expected_masked_data[2, 2, 2, 2, drop=FALSE])

  # Case 2: Slice access (drop=FALSE)
  expect_equal(lvs[1, , 1, 1, drop=FALSE], expected_masked_data[1, , 1, 1, drop=FALSE])
  expect_equal(lvs[, 2, , 2, drop=FALSE], expected_masked_data[, 2, , 2, drop=FALSE])
  expect_equal(lvs[, , 3, 3, drop=FALSE], expected_masked_data[, , 3, 3, drop=FALSE])
  expect_equal(lvs[1, 1, 1, , drop=FALSE], expected_masked_data[1, 1, 1, , drop=FALSE])

  # Case 3: Slice access (drop=TRUE)
  expect_equal(lvs[1, , 1, 1, drop=TRUE], expected_masked_data[1, , 1, 1]) # Drops 1st, 3rd, 4th dim
  expect_equal(lvs[, 2, 3, 2, drop=TRUE], expected_masked_data[, 2, 3, 2]) # Drops 2nd, 3rd, 4th dim
  expect_equal(lvs[1, 2, 3, , drop=TRUE], expected_masked_data[1, 2, 3, ]) # Drops 1st, 2nd, 3rd dim

  # Case 4: Range access
  expect_equal(lvs[1:2, 1:2, 1, 1:2, drop=FALSE], expected_masked_data[1:2, 1:2, 1, 1:2, drop=FALSE])
  expect_equal(lvs[1:X, Y, 1:Z, nVol, drop=FALSE], expected_masked_data[1:X, Y, 1:Z, nVol, drop=FALSE])

  # Case 5: Non-contiguous access
  idx_i <- c(1, X)
  idx_j <- c(1, Y)
  idx_k <- c(1, Z)
  idx_l <- c(1, nVol)
  expect_equal(lvs[idx_i, idx_j, idx_k, idx_l, drop=FALSE],
               expected_masked_data[idx_i, idx_j, idx_k, idx_l, drop=FALSE])

  # Case 6: Full dimension access (missing arguments)
  expect_equal(lvs[,,,1, drop=FALSE], expected_masked_data[,,,1, drop=FALSE]) # Full 3D volume
  expect_equal(lvs[1,,,, drop=FALSE], expected_masked_data[1,,,, drop=FALSE]) # Full Y, Z, L slice

  # Case 7: Test names() method added previously
  expect_equal(names(lvs), labels)

})

test_that("Write and read with special characters in labels", {
  # 1) Create small 4D data in memory
  set.seed(42)
  X <- 4; Y <- 5; Z <- 3; nVol <- 3  # small shape
  arr_data <- array(rnorm(X*Y*Z*nVol), dim=c(X,Y,Z,nVol))

  # 2) Create a NeuroSpace
  spc <- NeuroSpace(dim=c(X,Y,Z,nVol), spacing=c(1,1,1))

  # 3) Create a DenseNeuroVec
  vec_obj <- DenseNeuroVec(arr_data, space=spc)

  # 4) Create a mask
  mask_array <- array(sample(c(TRUE,FALSE), X*Y*Z, replace=TRUE, prob=c(0.6,0.4)), dim=c(X,Y,Z))
  mask_vol   <- LogicalNeuroVol(mask_array, space=drop_dim(spc))

  # 5) Create label names with special characters
  labels <- c("volume/1", "volume 2", "volume-3")

  # 6) Write to a temporary HDF5 file
  tmpfile <- tempfile(fileext=".h5")
  on.exit(unlink(tmpfile), add = TRUE)
  
  write_labeled_vec(vec=vec_obj, mask=mask_vol, labels=labels,
                    file=tmpfile, compression=4)
                    
  # 7) Read it back
  lvs <- read_labeled_vec(tmpfile, memoise=FALSE)
  
  # 8) Check that original labels are preserved
  expect_equal(lvs@labels, labels)
  
  # 9) Check actual data
  for (i in seq_along(labels)) {
    # Single volume via [[
    vol_i <- lvs[[i]]
    expect_s4_class(vol_i, "DenseNeuroVol")
    # Compare data in non-zero mask
    mask_idx <- which(mask_array)
    orig_1d  <- arr_data[,,,i][mask_idx]
    read_1d  <- as.array(vol_i)[mask_idx]
    # They should match up to a small tolerance
    expect_equal(orig_1d, read_1d, tolerance=1e-8)
  }
})

# Test for duplicate labels error
test_that("Duplicate labels are detected and reported", {
  # 1) Create small 4D data
  set.seed(42)
  X <- 4; Y <- 5; Z <- 3; nVol <- 2
  arr_data <- array(rnorm(X*Y*Z*nVol), dim=c(X,Y,Z,nVol))
  
  # 2) Create a NeuroSpace
  spc <- NeuroSpace(dim=c(X,Y,Z,nVol), spacing=c(1,1,1))
  
  # 3) Create a DenseNeuroVec
  vec_obj <- DenseNeuroVec(arr_data, space=spc)
  
  # 4) Create a mask
  mask_array <- array(sample(c(TRUE,FALSE), X*Y*Z, replace=TRUE, prob=c(0.6,0.4)), dim=c(X,Y,Z))
  mask_vol   <- LogicalNeuroVol(mask_array, space=drop_dim(spc))
  
  # 5) Create duplicate labels
  labels <- c("volume1", "volume1")
  
  # 6) Write to a temporary HDF5 file - should error
  tmpfile <- tempfile(fileext=".h5")
  on.exit(unlink(tmpfile), add = TRUE)
  
  # Expect this to error with duplicate labels
  expect_error(
    write_labeled_vec(vec=vec_obj, mask=mask_vol, labels=labels,
                      file=tmpfile, compression=4),
    "Duplicate labels detected"
  )
})

# Test for empty mask handling - should now error
test_that("Empty mask causes write_labeled_vec to error", {
  # 1) Create small 4D data
  set.seed(42)
  X <- 4; Y <- 5; Z <- 3; nVol <- 2
  arr_data <- array(rnorm(X*Y*Z*nVol), dim=c(X,Y,Z,nVol))
  
  # 2) Create a NeuroSpace
  spc <- NeuroSpace(dim=c(X,Y,Z,nVol), spacing=c(1,1,1))
  
  # 3) Create a DenseNeuroVec
  vec_obj <- DenseNeuroVec(arr_data, space=spc)

  # 4) Create an empty mask
  mask_array <- array(FALSE, dim=c(X,Y,Z))
  mask_vol   <- LogicalNeuroVol(mask_array, space=drop_dim(spc))

  # 5) Create labels
  labels <- c("volume1", "volume2")

  # 6) Write to a temporary HDF5 file - should now error
  tmpfile <- tempfile(fileext=".h5")
  on.exit(unlink(tmpfile), add = TRUE)

  # Expect an error because the mask is empty
  expect_error(
    write_labeled_vec(vec=vec_obj, mask=mask_vol, labels=labels,
                      file=tmpfile, compression=4),
    regexp = "Mask is empty"
  )

  # 7) Verify the file was NOT created (or is empty if error occurred after creation)
  expect_false(file.exists(tmpfile) && file.info(tmpfile)$size > 0)
})


test_that("write_labeled_vec → read_labeled_vec round-trip is loss-less", {
  skip_if_not_installed("hdf5r")
  tmp <- withr::local_tempfile(fileext = ".h5")

  toy <- make_toy_vec()
  vec   <- toy$vec
  mask  <- toy$mask
  lbls  <- paste0("vol_", seq_len(dim(vec)[4]))

  # ---------- write ----------
  
  write_labeled_vec(vec, mask, lbls, tmp,
                      compression = 6,       # make sure deflate path is used
                      chunk_size  = 5)       # exercise chunk logic
  

  # ---------- read ----------
  lset <- read_labeled_vec(tmp, memoise = TRUE)
  expect_s4_class(lset, "LabeledVolumeSet")
  expect_identical(lset@labels, lbls)

  # header / space
  expect_equal(dim(space(lset@mask)),
               dim(neuroim2::space(vec))[1:3])
  expect_equal(neuroim2::spacing(space(lset@mask)),
               neuroim2::spacing(neuroim2::space(vec))[1:3])

  # mask round-trip
  expect_identical(as.logical(as.array(lset@mask)),
                   as.logical(as.array(mask)))

  # each volume's voxels (only inside mask for speed)
  idx <- which(as.logical(as.array(mask)))
  for (v in seq_along(lbls)) {
    orig <- vec[idx + (v - 1) * prod(dim(mask))]
    read <- linear_access(lset,
                          idx + (v - 1) * prod(dim(mask)))
    expect_equal(read, orig, tolerance = 1e-10)
  }
})

test_that("label sanitisation & empty-mask corner cases behave correctly", {
  skip_if_not_installed("hdf5r")

  toy   <- make_toy_vec(nx = 2, ny = 2, nz = 1, nvol = 2)
  vec   <- toy$vec
  mask  <- toy$mask

  # ------------------------------------------------------------------
  # (a) duplicate after sanitisation  →  hard error
  dup_labels <- c("bad/label", "bad label")   # both sanitise to "bad_label"
  path <- paste0(withr::local_tempfile(), ".h5")
  expect_error(
    write_labeled_vec(vec, mask, dup_labels, path),
    "Duplicate labels",
    fixed = TRUE
  )

  # ------------------------------------------------------------------
  # (b) empty mask still produces valid zero volumes round-trip
  empty_mask <- neuroim2::LogicalNeuroVol(array(FALSE,
                                                dim = dim(mask)),
                                          space = space(mask))
  tmp <- withr::local_tempfile()
  tmp <- paste0(tmp, ".h5")
 
  expect_error(write_labeled_vec(vec, empty_mask,
                      labels = c("volA", "volB"), file = tmp),
               regexp = "Mask is empty")
  
})

# New test case for local_tempfile path is respected
test_that("local_tempfile path is respected", {
  skip_if_not_installed("hdf5r")

  # Create a DenseNeuroVec and an empty mask
  vec <- DenseNeuroVec(array(1:8, dim=c(2,2,1,2)), space=spc)
  mask_arr <- array(FALSE, dim=c(2,2,1))
  mask <- LogicalNeuroVol(mask_arr, drop_dim(spc))

  tmp <- withr::local_tempfile()
  tmp <- paste0(tmp, ".h5")
  # Expect error when writing with empty mask
  expect_error(
    write_labeled_vec(vec = vec,
                      mask = mask,
                      labels = c("volA", "volB"),
                      file = tmp),
    regexp = "Mask is empty"
  )

  # Verify the file wasn't created
  expect_false(file.exists(tmp))
})

