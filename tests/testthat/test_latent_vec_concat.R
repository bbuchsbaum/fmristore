# ------------------------------------------------------------------------------
# test_latent_vec_concat.R - Test suite for concat() method on LatentNeuroVec
# ------------------------------------------------------------------------------
library(neuroim2)
library(fmristore)

context("LatentNeuroVec :: concat()")

# Helper to create test objects with configurable time points
create_test_latent_vecs <- function(t1 = 10, t2 = 5, compatible = TRUE) {
  # Create spatial components - same for both objects to ensure compatibility
  dims_3d <- c(6, 6, 4)
  dims_4d_1 <- c(dims_3d, t1)
  dims_4d_2 <- c(dims_3d, t2)

  sp1 <- NeuroSpace(dims_4d_1)
  sp2 <- NeuroSpace(dims_4d_2)

  # Create same mask for both objects
  mask_arr <- array(TRUE, dim = dims_3d)
  # Make some voxels FALSE to test mask handling
  mask_arr[1:2, 1:2, 1] <- FALSE
  mask_vol <- LogicalNeuroVol(mask_arr, NeuroSpace(dims_3d))

  n_vox <- sum(mask_vol)
  k <- 3 # Number of components

  # Create basis matrices with appropriate time dimensions
  basis1 <- Matrix(rnorm(t1 * k), nrow = t1, ncol = k)
  basis2 <- Matrix(rnorm(t2 * k), nrow = t2, ncol = k)

  # Use same loadings for compatibility
  loadings <- Matrix(rnorm(n_vox * k), nrow = n_vox, ncol = k)

  # Create a second loadings matrix for incompatible test
  loadings2 <- if (compatible) loadings else Matrix(rnorm(n_vox * k), nrow = n_vox, ncol = k)

  # Create offset
  offset <- rnorm(n_vox)

  # Create LatentNeuroVec objects
  lvec1 <- LatentNeuroVec(
    basis = basis1,
    loadings = loadings,
    space = sp1,
    mask = mask_vol,
    offset = offset,
    label = "test_lvec_1"
  )

  lvec2 <- LatentNeuroVec(
    basis = basis2,
    loadings = loadings2, # Use either same or different loadings
    space = sp2,
    mask = mask_vol,
    offset = offset,
    label = "test_lvec_2"
  )

  list(lvec1 = lvec1, lvec2 = lvec2)
}

test_that("concat produces LatentNeuroVec when spatial components match", {
  # Create compatible test objects
  test_vecs <- create_test_latent_vecs(t1 = 10, t2 = 5, compatible = TRUE)
  lvec1 <- test_vecs$lvec1
  lvec2 <- test_vecs$lvec2

  # Concatenate and check result
  result <- concat(lvec1, lvec2)

  # Should return a LatentNeuroVec
  expect_s4_class(result, "LatentNeuroVec")

  # Check dimensions - time should be sum of input times
  expect_equal(dim(result)[4], 15) # 10 + 5 = 15

  # Check basis shape
  expect_equal(dim(result@basis), c(15, 3)) # 15 timepoints, 3 components

  # Check structure is preserved
  expect_equal(dim(space(result))[1:3], dim(space(lvec1))[1:3]) # Spatial dims match
  expect_equal(as.array(mask(result)), as.array(mask(lvec1))) # Mask unchanged

  # Check basis concatenation - first part should match lvec1
  expect_equal(
    as.matrix(result@basis[1:10, ]),
    as.matrix(lvec1@basis),
    tolerance = 1e-12
  )

  # Second part should match lvec2
  expect_equal(
    as.matrix(result@basis[11:15, ]),
    as.matrix(lvec2@basis),
    tolerance = 1e-12
  )

  # Loadings should be identical to the first object
  expect_equal(
    as.matrix(result@loadings),
    as.matrix(lvec1@loadings),
    tolerance = 1e-12
  )
})

test_that("concat falls back to NeuroVecSeq when spatial components don't match", {
  # Create incompatible test objects
  test_vecs <- create_test_latent_vecs(t1 = 8, t2 = 7, compatible = FALSE)
  lvec1 <- test_vecs$lvec1
  lvec2 <- test_vecs$lvec2

  # Concatenate and check result
  result <- concat(lvec1, lvec2)

  # Should return a NeuroVecSeq
  expect_s4_class(result, "NeuroVecSeq")

  # Check that both objects are stored
  expect_equal(length(result@vecs), 2)
})

test_that("concat handles mixed object types by falling back to NeuroVecSeq", {
  # Create a LatentNeuroVec
  test_vecs <- create_test_latent_vecs(t1 = 8, t2 = 7, compatible = TRUE)
  lvec1 <- test_vecs$lvec1

  # Create a different type of object (e.g., DenseNeuroVol)
  dims_3d <- c(6, 6, 4)
  sp <- NeuroSpace(dims_3d)
  vol_data <- array(rnorm(prod(dims_3d)), dim = dims_3d)
  dense_vol <- DenseNeuroVol(vol_data, sp)

  # Concatenate and check result
  result <- concat(lvec1, dense_vol)

  # Should return a NeuroVecSeq
  expect_s4_class(result, "DenseNeuroVec")
})

test_that("concat handles three or more objects correctly", {
  # Create three compatible test objects
  dims_3d <- c(6, 6, 4)
  sp1 <- NeuroSpace(c(dims_3d, 5))
  sp2 <- NeuroSpace(c(dims_3d, 3))
  sp3 <- NeuroSpace(c(dims_3d, 4))

  mask_arr <- array(TRUE, dim = dims_3d)
  mask_vol <- LogicalNeuroVol(mask_arr, NeuroSpace(dims_3d))

  n_vox <- sum(mask_vol)
  k <- 2 # Number of components

  loadings <- Matrix(rnorm(n_vox * k), nrow = n_vox, ncol = k)
  offset <- rnorm(n_vox)

  basis1 <- Matrix(rnorm(5 * k), nrow = 5, ncol = k)
  basis2 <- Matrix(rnorm(3 * k), nrow = 3, ncol = k)
  basis3 <- Matrix(rnorm(4 * k), nrow = 4, ncol = k)

  lvec1 <- LatentNeuroVec(basis1, loadings, sp1, mask_vol, offset, "lvec1")
  lvec2 <- LatentNeuroVec(basis2, loadings, sp2, mask_vol, offset, "lvec2")
  lvec3 <- LatentNeuroVec(basis3, loadings, sp3, mask_vol, offset, "lvec3")

  # Concatenate all three
  result <- concat(lvec1, lvec2, lvec3)

  # Should return a LatentNeuroVec
  expect_s4_class(result, "LatentNeuroVec")

  # Check dimensions - time should be sum of input times
  expect_equal(dim(result)[4], 12) # 5 + 3 + 4 = 12

  # Check basis shape
  expect_equal(dim(result@basis), c(12, 2)) # 12 timepoints, 2 components
})

test_that("concat falls back to NeuroVecSeq when masks are incompatible", {
  # Create two objects with different masks but otherwise compatible
  dims_3d <- c(6, 6, 4)
  k <- 3
  n_vox_total <- prod(dims_3d)

  # Mask 1 (more sparse)
  mask_arr1 <- array(TRUE, dim = dims_3d)
  mask_arr1[1:3, 1:3, 1:2] <- FALSE
  mask_vol1 <- LogicalNeuroVol(mask_arr1, NeuroSpace(dims_3d))
  n_vox1 <- sum(mask_vol1)

  # Mask 2 (less sparse)
  mask_arr2 <- array(TRUE, dim = dims_3d)
  mask_arr2[1, 1, 1] <- FALSE
  mask_vol2 <- LogicalNeuroVol(mask_arr2, NeuroSpace(dims_3d))
  n_vox2 <- sum(mask_vol2)

  # Create loadings compatible with *their respective* masks
  loadings1 <- Matrix(rnorm(n_vox1 * k), nrow = n_vox1, ncol = k)
  loadings2 <- Matrix(rnorm(n_vox2 * k), nrow = n_vox2, ncol = k) # Different loadings due to mask size

  # Create basis matrices
  t1 <- 5
  t2 <- 4
  basis1 <- Matrix(rnorm(t1 * k), nrow = t1, ncol = k)
  basis2 <- Matrix(rnorm(t2 * k), nrow = t2, ncol = k)

  # Create space objects
  sp1 <- NeuroSpace(c(dims_3d, t1))
  sp2 <- NeuroSpace(c(dims_3d, t2))

  # Create LatentNeuroVec objects
  lvec1 <- LatentNeuroVec(basis1, loadings1, sp1, mask_vol1, offset = numeric(0), "lvec_mask1")
  lvec2 <- LatentNeuroVec(basis2, loadings2, sp2, mask_vol2, offset = numeric(0), "lvec_mask2")

  # Concatenate - should fall back because masks (and thus loadings rows) differ
  result <- concat(lvec1, lvec2)

  # Should return a NeuroVecSeq
  expect_s4_class(result, "NeuroVecSeq")

  # Check that both original objects are stored within the sequence
  expect_equal(length(result@vecs), 2)
  expect_identical(result@vecs[[1]], lvec1)
  expect_identical(result@vecs[[2]], lvec2)
})
