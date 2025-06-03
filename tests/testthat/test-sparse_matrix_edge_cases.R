# Test suite for sparse matrix edge cases in LatentNeuroVec
# This addresses gaps in testing extreme sparsity and numerical edge cases

library(fmristore)

test_that("LatentNeuroVec handles extremely sparse data correctly", {
  skip_if_not_installed("Matrix")
  skip_if_not_installed("neuroim2")
  
  # Test 1: Extremely sparse loadings (>99% zeros)
  dims <- c(10L, 10L, 10L)
  n_time <- 50L
  n_comp <- 20L
  
  # Create mask with only 10% of voxels
  mask_array <- array(FALSE, dim = dims)
  n_mask_voxels <- 100L  # Out of 1000 total
  mask_indices <- sample(prod(dims), n_mask_voxels)
  mask_array[mask_indices] <- TRUE
  mask <- neuroim2::LogicalNeuroVol(mask_array, neuroim2::NeuroSpace(dims))
  
  # Create extremely sparse loadings - only 1% non-zero
  n_loadings_elements <- n_mask_voxels * n_comp
  n_nonzero <- ceiling(n_loadings_elements * 0.01)
  
  # Create sparse matrix directly
  i_indices <- sample(n_mask_voxels, n_nonzero, replace = TRUE)
  j_indices <- sample(n_comp, n_nonzero, replace = TRUE)
  values <- rnorm(n_nonzero, mean = 0, sd = 10)  # Large values to compensate for sparsity
  
  loadings <- Matrix::sparseMatrix(
    i = i_indices,
    j = j_indices,
    x = values,
    dims = c(n_mask_voxels, n_comp)
  )
  
  # Create basis
  basis <- matrix(rnorm(n_time * n_comp), nrow = n_time, ncol = n_comp)
  
  # Create LatentNeuroVec
  space_4d <- neuroim2::NeuroSpace(c(dims, n_time))
  lnv <- LatentNeuroVec(
    basis = basis,
    loadings = loadings,
    space = space_4d,
    mask = mask
  )
  
  # Test that operations work with extreme sparsity
  expect_s4_class(lnv, "LatentNeuroVec")
  expect_equal(dim(lnv), c(dims, n_time))
  
  # Test series extraction - should work even with mostly zeros
  series_data <- series(lnv, mask_indices[1:5])
  expect_equal(dim(series_data), c(n_time, 5))
  
  # Most values should be near zero due to sparsity
  expect_true(mean(abs(series_data)) < 1)
  
  # Test 2: Write and read back extremely sparse data
  temp_file <- tempfile(fileext = ".h5")
  on.exit(unlink(temp_file), add = TRUE)
  
  write_vec(lnv, temp_file, compression = 9)
  
  # Check file size - should be small due to sparsity
  file_size_mb <- file.info(temp_file)$size / 1024^2
  expect_true(file_size_mb < 1, 
              info = sprintf("Sparse file should be <1MB, got %.2f MB", file_size_mb))
  
  # Read back and verify
  h5file <- hdf5r::H5File$new(temp_file, mode = "r")
  
  # Check that sparse format was used
  expect_true(h5file$exists("basis/basis_matrix_sparse"),
              info = "Should use sparse storage for extremely sparse data")
  
  # Read sparse components
  sparse_grp <- h5file[["basis/basis_matrix_sparse"]]
  sparse_data <- sparse_grp[["data"]][]
  sparse_indices <- sparse_grp[["indices"]][]
  sparse_indptr <- sparse_grp[["indptr"]][]
  
  expect_equal(length(sparse_data), n_nonzero, tolerance = n_nonzero * 0.1,
               info = "Sparse data should have approximately same number of non-zeros")
  
  h5file$close()
  
  # Test 3: Edge case with all-zero components
  loadings_zero_col <- loadings
  loadings_zero_col[, 1] <- 0  # First component all zeros
  loadings_zero_col[, n_comp] <- 0  # Last component all zeros
  
  lnv_zero <- LatentNeuroVec(
    basis = basis,
    loadings = loadings_zero_col,
    space = space_4d,
    mask = mask
  )
  
  # Should handle all-zero components gracefully
  comp_list <- components(lnv_zero)
  expect_length(comp_list, n_comp)
  
  # First and last components should be all zeros
  expect_true(all(comp_list[[1]]@data == 0))
  expect_true(all(comp_list[[n_comp]]@data == 0))
})

test_that("LatentNeuroVec handles numerical edge cases", {
  skip_if_not_installed("Matrix")
  skip_if_not_installed("neuroim2")
  
  # Setup
  dims <- c(5L, 5L, 4L)
  n_time <- 10L
  n_comp <- 5L
  mask <- fmristore:::create_minimal_LogicalNeuroVol(dims)
  n_mask_voxels <- sum(mask@.Data)
  
  # Test 1: Inf and -Inf values in basis
  basis_inf <- matrix(rnorm(n_time * n_comp), nrow = n_time, ncol = n_comp)
  basis_inf[1, 1] <- Inf
  basis_inf[2, 2] <- -Inf
  
  loadings_normal <- Matrix::Matrix(
    matrix(rnorm(n_mask_voxels * n_comp), nrow = n_mask_voxels, ncol = n_comp),
    sparse = FALSE
  )
  
  space_4d <- neuroim2::NeuroSpace(c(dims, n_time))
  lnv_inf <- LatentNeuroVec(
    basis = basis_inf,
    loadings = loadings_normal,
    space = space_4d,
    mask = mask
  )
  
  # Operations should propagate Inf correctly
  series_inf <- series(lnv_inf, 1)
  expect_true(any(is.infinite(series_inf)))
  
  # Test 2: NaN values in loadings
  loadings_nan <- loadings_normal
  loadings_nan[1, 1] <- NaN
  loadings_nan[5, 3] <- NaN
  
  basis_normal <- matrix(rnorm(n_time * n_comp), nrow = n_time, ncol = n_comp)
  
  lnv_nan <- LatentNeuroVec(
    basis = basis_normal,
    loadings = loadings_nan,
    space = space_4d,
    mask = mask
  )
  
  # NaN should propagate through reconstruction
  series_nan <- series(lnv_nan, 1)  # First voxel has NaN in loadings
  expect_true(any(is.nan(series_nan)))
  
  # Test 3: Very large and very small values (numerical precision)
  basis_extreme <- basis_normal
  basis_extreme[, 1] <- basis_extreme[, 1] * 1e15  # Very large
  basis_extreme[, 2] <- basis_extreme[, 2] * 1e-15  # Very small
  
  loadings_extreme <- loadings_normal
  loadings_extreme[, 1] <- loadings_extreme[, 1] * 1e-15  # Compensate
  loadings_extreme[, 2] <- loadings_extreme[, 2] * 1e15   # Compensate
  
  lnv_extreme <- LatentNeuroVec(
    basis = basis_extreme,
    loadings = loadings_extreme,
    space = space_4d,
    mask = mask,
    offset = rep(1e10, n_mask_voxels)  # Large offset too
  )
  
  # Should handle extreme values without overflow
  series_extreme <- series(lnv_extreme, 1:5)
  expect_false(any(is.nan(series_extreme)),
               info = "Should not produce NaN from extreme but valid values")
  
  # Test 4: Write/read with special values
  temp_file <- tempfile(fileext = ".h5")
  on.exit(unlink(temp_file), add = TRUE)
  
  # Create data with mix of special values
  basis_special <- matrix(0, nrow = n_time, ncol = n_comp)
  basis_special[1:5, 1] <- c(Inf, -Inf, NaN, 1e20, -1e-20)
  basis_special[, 2:n_comp] <- rnorm((n_comp-1) * n_time)
  
  loadings_special <- Matrix::Matrix(
    matrix(rnorm(n_mask_voxels * n_comp), nrow = n_mask_voxels, ncol = n_comp)
  )
  loadings_special[1, 1] <- NaN
  
  lnv_special <- LatentNeuroVec(
    basis = basis_special,
    loadings = loadings_special,
    space = space_4d,
    mask = mask
  )
  
  # Should write without error
  expect_silent(write_vec(lnv_special, temp_file))
  
  # Should read back with special values preserved
  h5file <- hdf5r::H5File$new(temp_file, mode = "r")
  basis_read <- h5file[["scans/embedding"]][]
  h5file$close()
  
  expect_true(is.infinite(basis_read[1, 1]) && basis_read[1, 1] > 0)
  expect_true(is.infinite(basis_read[2, 1]) && basis_read[2, 1] < 0)
  expect_true(is.nan(basis_read[3, 1]))
})

test_that("Sparse matrix format conversions and memory efficiency", {
  skip_if_not_installed("Matrix")
  skip_if_not_installed("neuroim2")
  
  # Test different sparse matrix formats
  dims <- c(20L, 20L, 10L)
  n_time <- 100L
  n_comp <- 50L
  
  # Create mask with 1000 voxels
  mask_array <- array(FALSE, dim = dims)
  mask_indices <- sample(prod(dims), min(1000L, prod(dims)))
  mask_array[mask_indices] <- TRUE
  mask <- neuroim2::LogicalNeuroVol(mask_array, neuroim2::NeuroSpace(dims))
  n_mask_voxels <- sum(mask@.Data)
  
  # Test 1: CSC format (column sparse)
  loadings_csc <- Matrix::sparseMatrix(
    i = sample(n_mask_voxels, 100, replace = TRUE),
    j = sample(n_comp, 100, replace = TRUE),
    x = rnorm(100),
    dims = c(n_mask_voxels, n_comp),
    repr = "C"  # CSC format
  )
  
  basis <- matrix(rnorm(n_time * n_comp), nrow = n_time, ncol = n_comp)
  space_4d <- neuroim2::NeuroSpace(c(dims, n_time))
  
  lnv_csc <- LatentNeuroVec(basis = basis, loadings = loadings_csc, space = space_4d, mask = mask)
  expect_s4_class(lnv_csc, "LatentNeuroVec")
  
  # Test 2: CSR format (row sparse) - should be converted internally
  loadings_csr <- as(loadings_csc, "RsparseMatrix")
  
  lnv_csr <- LatentNeuroVec(basis = basis, loadings = loadings_csr, space = space_4d, mask = mask)
  expect_s4_class(lnv_csr, "LatentNeuroVec")
  
  # Results should be identical
  series_csc <- series(lnv_csc, 1:10)
  series_csr <- series(lnv_csr, 1:10)
  expect_equal(series_csc, series_csr, tolerance = 1e-10)
  
  # Test 3: Triplet format
  loadings_triplet <- as(loadings_csc, "TsparseMatrix")
  
  lnv_triplet <- LatentNeuroVec(basis = basis, loadings = loadings_triplet, space = space_4d, mask = mask)
  expect_s4_class(lnv_triplet, "LatentNeuroVec")
  
  # Test 4: Memory efficiency for very sparse data
  # Create increasingly sparse matrices and check memory usage
  sparsity_levels <- c(0.001, 0.01, 0.1)  # 0.1%, 1%, 10% non-zero
  
  for (sparsity in sparsity_levels) {
    n_nonzero <- ceiling(n_mask_voxels * n_comp * sparsity)
    
    loadings_sparse <- Matrix::sparseMatrix(
      i = sample(n_mask_voxels, n_nonzero, replace = TRUE),
      j = sample(n_comp, n_nonzero, replace = TRUE),
      x = rnorm(n_nonzero),
      dims = c(n_mask_voxels, n_comp)
    )
    
    # Check memory usage of sparse vs dense
    sparse_size <- object.size(loadings_sparse)
    dense_size <- object.size(as.matrix(loadings_sparse))
    
    expect_true(
      sparse_size < dense_size * 0.5,
      info = sprintf("Sparse matrix with %.1f%% non-zeros should use less than half the memory of dense",
                     sparsity * 100)
    )
    
    # Test write efficiency
    temp_file <- tempfile(fileext = ".h5")
    on.exit(unlink(temp_file), add = TRUE)
    
    lnv_sparse <- LatentNeuroVec(basis = basis, loadings = loadings_sparse, space = space_4d, mask = mask)
    write_vec(lnv_sparse, temp_file, compression = 6)
    
    # Verify sparse storage was used for very sparse data
    if (sparsity <= 0.01) {  # 1% or less
      h5file <- hdf5r::H5File$new(temp_file, mode = "r")
      expect_true(
        h5file$exists("basis/basis_matrix_sparse"),
        info = sprintf("Should use sparse storage for %.1f%% sparsity", sparsity * 100)
      )
      h5file$close()
    }
    
    unlink(temp_file)
  }
  
  # Test 5: Conversion during read - dense to sparse threshold
  # Write a moderately sparse matrix as dense, ensure it can still be read
  loadings_moderate <- Matrix::Matrix(0, nrow = n_mask_voxels, ncol = n_comp, sparse = FALSE)
  # Add ~30% non-zeros
  n_elements <- ceiling(n_mask_voxels * n_comp * 0.3)
  indices <- sample(length(loadings_moderate), n_elements)
  loadings_moderate[indices] <- rnorm(n_elements)
  
  lnv_moderate <- LatentNeuroVec(basis = basis, loadings = loadings_moderate, space = space_4d, mask = mask)
  
  temp_file <- tempfile(fileext = ".h5")
  on.exit(unlink(temp_file), add = TRUE)
  write_vec(lnv_moderate, temp_file)
  
  # Should be stored as dense (30% is above typical sparse threshold)
  h5file <- hdf5r::H5File$new(temp_file, mode = "r")
  expect_true(h5file$exists("basis/basis_matrix"))
  expect_false(h5file$exists("basis/basis_matrix_sparse"))
  h5file$close()
})