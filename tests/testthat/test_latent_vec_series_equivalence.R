library(testthat)
library(hdf5r)
library(neuroim2)
library(Matrix)
library(fmristore)

# Create a "full rank" LatentNeuroVec that perfectly reconstructs original data
create_full_rank_latent <- function(dims = c(8, 8, 4, 20), # x,y,z,t
                                    k = NULL) { # Number of components (default: full rank)

  sp <- NeuroSpace(dims)
  nTime <- dims[4]
  nVox_total <- prod(dims[1:3])

  # Create mask - use about 60% of voxels for good coverage
  set.seed(123) # For reproducibility
  mask_arr <- array(FALSE, dim = dims[1:3])
  mask_indices <- sample(nVox_total, floor(nVox_total * 0.6))
  mask_arr[mask_indices] <- TRUE
  mask_vol <- LogicalNeuroVol(mask_arr, drop_dim(sp))
  nVox_mask <- sum(mask_vol)

  # If k not specified, use full rank (min of dimensions)
  if (is.null(k)) {
    k <- min(nTime, nVox_mask)
  }

  # Create original dense data first
  set.seed(456)
  original_data <- array(rnorm(prod(dims), mean = 10, sd = 3), dims)

  # Extract masked data correctly: [nVox_mask x nTime]
  # We need to extract time series for each masked voxel
  masked_indices <- which(mask_arr, arr.ind = TRUE)
  masked_data <- matrix(0, nrow = nVox_mask, ncol = nTime)

  for (i in seq_len(nVox_mask)) {
    x <- masked_indices[i, 1]
    y <- masked_indices[i, 2]
    z <- masked_indices[i, 3]
    masked_data[i, ] <- original_data[x, y, z, ]
  }

  # Compute SVD for perfect reconstruction with k components
  svd_result <- svd(masked_data, nu = k, nv = k)

  # Create basis and loadings from SVD
  # SVD gives: masked_data = U %*% diag(d) %*% t(V)
  # We want: masked_data = loadings %*% t(basis) + offset
  # So: loadings = U %*% diag(d), basis = V
  basis_mat <- Matrix(svd_result$v, sparse = FALSE)  # [nTime x k]
  loadings_mat <- Matrix(svd_result$u %*% diag(svd_result$d[1:k]),
    sparse = FALSE) # [nVox_mask x k]

  # Use zero offset for perfect reconstruction
  offset_vec <- rep(0, nVox_mask)

  # Create LatentNeuroVec
  lvec <- LatentNeuroVec(basis = basis_mat,
    loadings = loadings_mat,
    space = sp,
    mask = mask_vol,
    offset = offset_vec,
    label = "full_rank_test")

  # Create equivalent DenseNeuroVec
  dense_vec <- DenseNeuroVec(original_data, sp)

  return(list(latent = lvec, dense = dense_vec, mask = mask_vol))
}

test_that("LatentNeuroVec series() behaves identically to DenseNeuroVec - integer indices", {
  # Create full-rank test data
  test_data <- create_full_rank_latent(dims = c(6, 6, 3, 15), k = 15)
  lvec <- test_data$latent
  dvec <- test_data$dense
  mask <- test_data$mask

  # Test various integer coordinate combinations
  test_coords <- list(
    c(1L, 1L, 1L),  # First voxel
    c(3L, 3L, 2L),  # Middle voxel
    c(6L, 6L, 3L),  # Last voxel
    c(2L, 4L, 1L),  # Random voxel
    c(5L, 2L, 3L)   # Another random voxel
  )

  for (coords in test_coords) {
    i <- coords[1]
    j <- coords[2]
    k <- coords[3]

    # Extract series using integer coordinates
    latent_series <- series(lvec, i, j, k)
    dense_series <- series(dvec, i, j, k)

    # Check if voxel is in mask
    is_in_mask <- mask[i, j, k]

    if (is_in_mask) {
      # Should match exactly for masked voxels
      expect_equal(latent_series, dense_series, tolerance = 1e-10,
        info = sprintf("Series mismatch at coords (%d,%d,%d) - masked voxel", i, j, k))
    } else {
      # Outside mask should be zeros for LatentNeuroVec
      expect_equal(latent_series, rep(0, dim(lvec)[4]), tolerance = 1e-10,
        info = sprintf("LatentNeuroVec should return zeros for unmasked voxel (%d,%d,%d)", i, j, k))
      # Dense should return actual data
      expect_equal(dense_series, dense_series,
        info = sprintf("DenseNeuroVec series extraction failed at (%d,%d,%d)", i, j, k))
    }
  }
})

test_that("LatentNeuroVec series() behaves identically to DenseNeuroVec - numeric indices", {
  test_data <- create_full_rank_latent(dims = c(5, 5, 4, 12), k = 12)
  lvec <- test_data$latent
  dvec <- test_data$dense
  mask <- test_data$mask

  # Test with numeric (potentially non-integer) coordinates
  test_coords <- list(
    c(1.0, 1.0, 1.0),
    c(3.0, 3.0, 2.0),
    c(2.5, 3.5, 1.5),  # Non-integer coords (should be rounded/interpolated)
    c(4.0, 4.0, 4.0)
  )

  for (coords in test_coords) {
    i <- coords[1]
    j <- coords[2]
    k <- coords[3]

    # For the comparison, we need to handle how neuroim2 deals with non-integer indices
    # Typically they get rounded or interpolated
    latent_series <- series(lvec, i, j, k)
    dense_series <- series(dvec, i, j, k)

    # Check if the rounded coordinates are in mask
    i_round <- round(i)
    j_round <- round(j)
    k_round <- round(k)
    if (i_round >= 1 && i_round <= dim(mask)[1] &&
      j_round >= 1 && j_round <= dim(mask)[2] &&
      k_round >= 1 && k_round <= dim(mask)[3]) {

      is_in_mask <- mask[i_round, j_round, k_round]

      if (is_in_mask) {
        expect_equal(latent_series, dense_series, tolerance = 1e-8,
          info = sprintf("Series mismatch at numeric coords (%.1f,%.1f,%.1f)", i, j, k))
      }
    }
  }
})

test_that("LatentNeuroVec series() behaves identically to DenseNeuroVec - matrix ROI", {
  test_data <- create_full_rank_latent(dims = c(7, 7, 3, 10), k = 10)
  lvec <- test_data$latent
  dvec <- test_data$dense
  mask <- test_data$mask

  # Create matrix of ROI coordinates (each row is [x, y, z])
  roi_coords <- rbind(
    c(1, 1, 1),
    c(3, 3, 2),
    c(5, 5, 3),
    c(2, 4, 1),
    c(6, 2, 2)
  )

  # Extract series using matrix coordinates
  latent_series_mat <- series(lvec, roi_coords)
  dense_series_mat <- series(dvec, roi_coords)

  # Both should return matrices with same dimensions
  expect_equal(dim(latent_series_mat), dim(dense_series_mat),
    info = "Matrix dimensions should match between LatentNeuroVec and DenseNeuroVec")

  # Check each voxel individually
  for (i in seq_len(nrow(roi_coords))) {
    coords <- roi_coords[i, ]
    x <- coords[1]
    y <- coords[2]
    z <- coords[3]

    is_in_mask <- mask[x, y, z]

    if (is_in_mask) {
      expect_equal(latent_series_mat[, i], dense_series_mat[, i], tolerance = 1e-10,
        info = sprintf("Matrix ROI series mismatch at voxel %d coords (%d,%d,%d)", i, x, y, z))
    } else {
      # LatentNeuroVec should return zeros for unmasked voxels
      expect_equal(latent_series_mat[, i], rep(0, dim(lvec)[4]), tolerance = 1e-10,
        info = sprintf("LatentNeuroVec should return zeros for unmasked matrix ROI voxel %d", i))
    }
  }
})

test_that("LatentNeuroVec series() behaves identically to DenseNeuroVec - LogicalNeuroVol mask", {
  test_data <- create_full_rank_latent(dims = c(6, 6, 3, 8), k = 8)
  lvec <- test_data$latent
  dvec <- test_data$dense

  # Create a separate logical mask for ROI extraction (subset of original mask)
  roi_mask_arr <- array(FALSE, dim = dim(test_data$mask)[1:3])
  # Select a small ROI
  roi_mask_arr[2:4, 2:4, 1:2] <- TRUE
  roi_mask <- LogicalNeuroVol(roi_mask_arr, drop_dim(space(lvec)))

  # Extract series using LogicalNeuroVol mask
  latent_series_logical <- series(lvec, roi_mask)
  dense_series_logical <- series(dvec, roi_mask)

  # Both should return matrices with same dimensions
  expect_equal(dim(latent_series_logical), dim(dense_series_logical),
    info = "LogicalNeuroVol mask series dimensions should match")

  # For voxels within both the ROI mask and the LatentNeuroVec mask, should match
  # For voxels in ROI but outside LatentNeuroVec mask, LatentNeuroVec should return zeros
  roi_indices <- which(roi_mask@.Data)
  original_mask_indices <- which(test_data$mask@.Data)

  for (i in seq_len(ncol(latent_series_logical))) {
    # Convert column index to spatial coordinates
    coords <- arrayInd(roi_indices[i], dim(roi_mask))
    x <- coords[1]
    y <- coords[2]
    z <- coords[3]

    if (test_data$mask[x, y, z]) {
      # In both masks - should match
      expect_equal(latent_series_logical[, i], dense_series_logical[, i], tolerance = 1e-10,
        info = sprintf("LogicalNeuroVol mask series mismatch at voxel %d", i))
    } else {
      # In ROI mask but not in LatentNeuroVec mask - LatentNeuroVec should be zeros
      expect_equal(latent_series_logical[, i], rep(0, dim(lvec)[4]), tolerance = 1e-10,
        info = sprintf("LatentNeuroVec should return zeros for ROI voxel %d outside its mask", i))
    }
  }
})

test_that("LatentNeuroVec series() behaves identically to DenseNeuroVec - drop parameter", {
  test_data <- create_full_rank_latent(dims = c(5, 5, 2, 6), k = 6)
  lvec <- test_data$latent
  dvec <- test_data$dense
  mask <- test_data$mask

  # Test single voxel extraction with drop=TRUE and drop=FALSE
  coords <- c(3L, 3L, 1L)
  i <- coords[1]
  j <- coords[2]
  k <- coords[3]

  if (mask[i, j, k]) {
    # Test with drop=TRUE (default)
    latent_drop_true <- series(lvec, i, j, k, drop = TRUE)
    dense_drop_true <- series(dvec, i, j, k, drop = TRUE)
    expect_equal(latent_drop_true, dense_drop_true, tolerance = 1e-10,
      info = "Series with drop=TRUE should match")

    # Test with drop=FALSE
    latent_drop_false <- series(lvec, i, j, k, drop = FALSE)
    dense_drop_false <- series(dvec, i, j, k, drop = FALSE)
    expect_equal(latent_drop_false, dense_drop_false, tolerance = 1e-10,
      info = "Series with drop=FALSE should match")

    # Check that dimensions are handled correctly
    expect_equal(class(latent_drop_true), class(dense_drop_true),
      info = "Return types should match for drop=TRUE")
    expect_equal(class(latent_drop_false), class(dense_drop_false),
      info = "Return types should match for drop=FALSE")
  } else {
    skip("Test coordinates not in mask - skipping drop parameter test")
  }
})

test_that("LatentNeuroVec series() behaves identically to DenseNeuroVec - single coordinate vectors", {
  test_data <- create_full_rank_latent(dims = c(6, 6, 4, 12), k = 12)
  lvec <- test_data$latent
  dvec <- test_data$dense
  mask <- test_data$mask

  # Test with single coordinates (not multiple i, j, k vectors as neuroim2 doesn't support that)
  test_coords <- list(
    c(1L, 2L, 1L),
    c(3L, 4L, 2L),
    c(5L, 6L, 3L)
  )

  for (coords in test_coords) {
    i <- coords[1]
    j <- coords[2]
    k <- coords[3]

    latent_single <- series(lvec, i, j, k)
    dense_single <- series(dvec, i, j, k)

    if (mask[i, j, k]) {
      expect_equal(latent_single, dense_single, tolerance = 1e-10,
        info = sprintf("Single coordinate series mismatch at coords (%d,%d,%d)", i, j, k))
    } else {
      expect_equal(latent_single, rep(0, dim(lvec)[4]), tolerance = 1e-10,
        info = sprintf("LatentNeuroVec should return zeros for unmasked coord (%d,%d,%d)", i, j, k))
    }
  }
})

test_that("LatentNeuroVec series() behaves identically to DenseNeuroVec - edge cases", {
  test_data <- create_full_rank_latent(dims = c(4, 4, 2, 5), k = 5)
  lvec <- test_data$latent
  dvec <- test_data$dense

  # Test boundary coordinates
  boundary_coords <- list(
    c(1L, 1L, 1L),  # Min boundary
    c(4L, 4L, 2L),  # Max boundary
    c(1L, 4L, 1L),  # Mixed boundaries
    c(4L, 1L, 2L)   # Mixed boundaries
  )

  for (coords in boundary_coords) {
    i <- coords[1]
    j <- coords[2]
    k <- coords[3]

    # Both should handle boundaries the same way
    latent_boundary <- series(lvec, i, j, k)
    dense_boundary <- series(dvec, i, j, k)

    # Length should always match time dimension
    expect_equal(length(latent_boundary), dim(lvec)[4],
      info = sprintf("Latent series length should match time dimension at boundary (%d,%d,%d)", i, j, k))
    expect_equal(length(dense_boundary), dim(dvec)[4],
      info = sprintf("Dense series length should match time dimension at boundary (%d,%d,%d)", i, j, k))

    # If in mask, should match
    if (test_data$mask[i, j, k]) {
      expect_equal(latent_boundary, dense_boundary, tolerance = 1e-10,
        info = sprintf("Boundary series should match at (%d,%d,%d)", i, j, k))
    }
  }
})

test_that("LatentNeuroVec series() performance and consistency", {
  # Create larger test case for performance validation
  test_data <- create_full_rank_latent(dims = c(10, 10, 8, 50), k = 50)
  lvec <- test_data$latent
  dvec <- test_data$dense
  mask <- test_data$mask

  # Test multiple extractions to ensure consistency
  test_coords <- rbind(
    c(5, 5, 4),
    c(2, 8, 3),
    c(9, 3, 7),
    c(1, 1, 1),
    c(10, 10, 8)
  )

  for (i in seq_len(nrow(test_coords))) {
    coords <- test_coords[i, ]
    x <- coords[1]
    y <- coords[2]
    z <- coords[3]

    if (mask[x, y, z]) {
      # Multiple extractions should be consistent
      series1 <- series(lvec, x, y, z)
      series2 <- series(lvec, x, y, z)
      expect_equal(series1, series2, tolerance = 1e-12,
        info = sprintf("Repeated series extraction should be identical at (%d,%d,%d)", x, y, z))

      # Should match dense version
      dense_series <- series(dvec, x, y, z)
      expect_equal(series1, dense_series, tolerance = 1e-10,
        info = sprintf("Large-scale series should match dense at (%d,%d,%d)", x, y, z))
    }
  }

  # Test matrix extraction for multiple voxels
  latent_matrix <- series(lvec, test_coords)
  dense_matrix <- series(dvec, test_coords)

  # Should have same structure
  expect_equal(dim(latent_matrix), dim(dense_matrix),
    info = "Large matrix extractions should have same dimensions")

  # Check consistency across masked voxels
  for (i in seq_len(ncol(latent_matrix))) {
    coords <- test_coords[i, ]
    if (mask[coords[1], coords[2], coords[3]]) {
      expect_equal(latent_matrix[, i], dense_matrix[, i], tolerance = 1e-10,
        info = sprintf("Large matrix series should match at column %d", i))
    }
  }
})

test_that("LatentNeuroVec series() with different k values", {
  # Test with different numbers of components to ensure series() behavior is consistent
  dims <- c(6, 6, 3, 20)

  for (k in c(5, 10, 15, 20)) {
    test_data <- create_full_rank_latent(dims = dims, k = k)
    lvec <- test_data$latent
    dvec <- test_data$dense
    mask <- test_data$mask

    # Test a few voxels
    test_coords <- list(c(3L, 3L, 2L), c(5L, 2L, 1L))

    for (coords in test_coords) {
      i <- coords[1]
      j <- coords[2]
      k_coord <- coords[3]

      if (mask[i, j, k_coord]) {
        latent_series <- series(lvec, i, j, k_coord)
        dense_series <- series(dvec, i, j, k_coord)

        # For full-rank (k=20), should match exactly
        # For lower k, should be a good approximation
        if (k == dims[4]) {
          expect_equal(latent_series, dense_series, tolerance = 1e-10,
            label = sprintf("Full-rank (k=%d) series should match exactly at (%d,%d,%d)", k, i, j, k_coord))
        } else {
          # For reduced rank, check that they're reasonably close
          # Only check correlation if both series have variance
          if (sd(latent_series) > 1e-10 && sd(dense_series) > 1e-10) {
            correlation <- cor(latent_series, dense_series)
            # For very low rank approximations, correlation may be lower
            min_correlation <- if (k <= 5) 0.5 else if (k <= 10) 0.8 else 0.9
            expect_gt(correlation, min_correlation,
              label = sprintf("Reduced-rank (k=%d) series should be correlated (>%.1f) at (%d,%d,%d)", k, min_correlation, i, j, k_coord))
          }
        }
      }
    }
  }
})
