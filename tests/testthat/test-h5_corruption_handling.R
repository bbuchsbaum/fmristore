# Test suite for HDF5 file corruption and invalid format handling
# This addresses a critical gap in testing robustness against malformed files

library(fmristore)

test_that("H5NeuroVol handles corrupted HDF5 files gracefully", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("neuroim2")
  
  # Create a valid H5NeuroVol file first
  temp_file <- tempfile(fileext = ".h5")
  on.exit(unlink(temp_file), add = TRUE)
  
  # Create valid file
  dims <- c(10L, 10L, 5L)
  vol_data <- array(rnorm(prod(dims)), dim = dims)
  space_obj <- neuroim2::NeuroSpace(dims)
  vol <- neuroim2::NeuroVol(vol_data, space_obj)
  
  # Write valid H5NeuroVol
  h5_vol <- as_h5(vol, file = temp_file)
  close(h5_vol)
  
  # Test 1: Missing required dataset
  h5file <- hdf5r::H5File$new(temp_file, mode = "r+")
  h5file$link_delete("data/elements")  # Remove critical data
  h5file$close()
  
  expect_error(
    H5NeuroVol(temp_file),
    regexp = "data/elements.*not found|does not exist",
    info = "Should fail when data/elements dataset is missing"
  )
  
  # Test 2: Wrong data type in space dimensions
  h5file <- hdf5r::H5File$new(temp_file, mode = "r+")
  if (h5file$exists("space/dim")) {
    h5file$link_delete("space/dim")
  }
  # Create dimension with wrong type (string instead of integer)
  h5file[["space/dim"]] <- c("10", "10", "5")
  h5file$close()
  
  expect_error(
    H5NeuroVol(temp_file),
    info = "Should fail when space/dim has wrong data type"
  )
  
  # Test 3: Inconsistent dimensions between space and data
  temp_file2 <- tempfile(fileext = ".h5")
  on.exit(unlink(temp_file2), add = TRUE)
  
  h5file <- hdf5r::H5File$new(temp_file2, mode = "w")
  h5file$create_group("space")
  h5file$create_group("data")
  
  # Write inconsistent dimensions
  h5file[["space/dim"]] <- c(10L, 10L, 5L)
  h5file[["space/origin"]] <- c(0, 0, 0)
  h5file[["space/trans"]] <- diag(4)
  
  # Write data with different dimensions
  wrong_data <- array(rnorm(20*20*10), dim = c(20L, 20L, 10L))
  h5file[["data/elements"]] <- wrong_data
  hdf5r::h5attr(h5file, "rtype") <- "DenseNeuroVol"
  h5file$close()
  
  expect_error(
    H5NeuroVol(temp_file2),
    info = "Should fail when data dimensions don't match space dimensions"
  )
})

test_that("H5ClusterExperiment handles invalid cluster configurations", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("neuroim2")
  
  temp_file <- tempfile(fileext = ".h5")
  on.exit(unlink(temp_file), add = TRUE)
  
  # Create a basic valid structure with enough voxels for 3 clusters
  # Create a mask with more true voxels
  mask_array <- array(FALSE, dim = c(5L, 5L, 5L))
  # Set 20 voxels to TRUE
  mask_indices <- sample(prod(c(5L, 5L, 5L)), 20)
  mask_array[mask_indices] <- TRUE
  mask <- neuroim2::LogicalNeuroVol(mask_array, neuroim2::NeuroSpace(c(5L, 5L, 5L)))
  clusters <- fmristore:::create_minimal_ClusteredNeuroVol(mask, num_clusters = 3L)
  
  # Write initial valid file
  h5file <- hdf5r::H5File$new(temp_file, mode = "w")
  
  # Write header
  header_grp <- h5file$create_group("header")
  header_grp[["dim"]] <- c(4L, 5L, 5L, 5L, 10L, 1L, 1L, 1L)
  header_grp[["pixdim"]] <- c(0, 3, 3, 3, 2, 0, 0, 0)
  
  # Write mask
  h5file[["mask"]] <- as.integer(mask@.Data)
  
  # Write cluster_map with values
  cluster_ids <- clusters@clusters
  voxel_indices <- which(mask@.Data)
  cluster_map_full <- integer(length(mask@.Data))
  cluster_map_full[voxel_indices] <- cluster_ids
  h5file[["cluster_map"]] <- cluster_map_full[voxel_indices]
  
  # Write scans group with inconsistent cluster IDs
  scans_grp <- h5file$create_group("scans")
  run1_grp <- scans_grp$create_group("run1")
  clusters_grp <- run1_grp$create_group("clusters")
  
  # Create clusters that don't match cluster_map
  # cluster_map has IDs 1,2,3 but we'll create datasets for 4,5,6
  clusters_grp[["cluster_4"]] <- matrix(rnorm(10 * 5), nrow = 5, ncol = 10)
  clusters_grp[["cluster_5"]] <- matrix(rnorm(8 * 5), nrow = 5, ncol = 10)
  clusters_grp[["cluster_6"]] <- matrix(rnorm(7 * 5), nrow = 5, ncol = 10)
  
  h5file$close()
  
  # Test: Should detect mismatch between cluster_map and actual cluster datasets
  expect_error(
    H5ClusterExperiment(temp_file),
    info = "Should fail when cluster IDs in data don't match cluster_map"
  )
  
  # Test 2: Corrupt cluster metadata
  temp_file2 <- tempfile(fileext = ".h5")
  on.exit(unlink(temp_file2), add = TRUE)
  
  # Create valid initial structure
  fmristore:::create_minimal_h5_for_H5ClusterExperiment(
    file_path = temp_file2,
    mask_dims = c(5L, 5L, 5L),
    num_clusters = 3L,
    n_time_run1 = 10L
  )
  
  # Corrupt the cluster metadata
  h5file <- hdf5r::H5File$new(temp_file2, mode = "r+")
  if (h5file$exists("clusters")) {
    clusters_grp <- h5file[["clusters"]]
    if (clusters_grp$exists("cluster_meta")) {
      clusters_grp$link_delete("cluster_meta")
      # Write invalid metadata structure
      clusters_grp[["cluster_meta"]] <- list(
        invalid = "structure",
        not_a = "data.frame"
      )
    }
  }
  h5file$close()
  
  # Should handle gracefully (warning but not error, since cluster_meta is optional)
  expect_warning(
    H5ClusterExperiment(temp_file2),
    regexp = "cluster_meta",
    info = "Should warn about invalid cluster metadata structure"
  )
})

test_that("Concurrent HDF5 access and resource cleanup", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("neuroim2")
  
  temp_file <- tempfile(fileext = ".h5")
  on.exit(unlink(temp_file), add = TRUE)
  
  # Create a test H5NeuroVec file
  dims <- c(5L, 5L, 5L, 20L)
  vec_data <- array(rnorm(prod(dims)), dim = dims)
  space_obj <- neuroim2::NeuroSpace(dims)
  vec <- neuroim2::NeuroVec(vec_data, space_obj)
  h5_vec <- as_h5(vec, file = temp_file)
  close(h5_vec)
  
  # Test 1: Multiple read handles to same file
  h5_vec1 <- H5NeuroVec(temp_file)
  h5_vec2 <- H5NeuroVec(temp_file)
  h5_vec3 <- H5NeuroVec(temp_file)
  
  # All should be able to read
  expect_equal(dim(h5_vec1), dims)
  expect_equal(dim(h5_vec2), dims)
  expect_equal(dim(h5_vec3), dims)
  
  # Access data from different handles using linear_access
  data1 <- linear_access(h5_vec1, 1:5)
  data2 <- linear_access(h5_vec2, 6:10)
  data3 <- linear_access(h5_vec3, 11:15)
  
  expect_length(data1, 5)
  expect_length(data2, 5)
  expect_length(data3, 5)
  
  # Close in different order
  close(h5_vec2)
  close(h5_vec1)
  close(h5_vec3)
  
  # Test 2: Resource cleanup after error
  h5_vec <- H5NeuroVec(temp_file)
  
  # Force an error during operation
  expect_error({
    # Try to access invalid indices
    linear_access(h5_vec, prod(dims) + 100)
  })
  
  # File handle should still be valid and closeable
  expect_true(h5_vec@obj$is_valid)
  close(h5_vec)
  
  # Test 3: Cleanup in error conditions during write
  temp_file2 <- tempfile(fileext = ".h5")
  on.exit(unlink(temp_file2), add = TRUE)
  
  # Create a situation that will fail during write
  bad_vec <- vec
  # Inject NaN values that might cause issues with certain compression
  bad_vec@.Data[1:100] <- NaN
  
  # Should handle NaN values appropriately
  h5_bad <- as_h5(bad_vec, file = temp_file2, compression = 9)
  close(h5_bad)
  
  # Verify file was written and can be read back
  h5_check <- H5NeuroVec(temp_file2)
  check_data <- h5_check[1, 1, 1, 1:10]
  expect_true(is.nan(check_data[1]))
  close(h5_check)
  
  # Test 4: Handle cleanup with gc()
  # Create handles but don't explicitly close them
  for (i in 1:3) {
    h5_temp <- H5NeuroVec(temp_file)
    # Intentionally not closing - relying on finalizers
  }
  
  # Force garbage collection
  gc()
  Sys.sleep(0.1)  # Give finalizers time to run
  
  # Should still be able to open the file
  h5_final <- H5NeuroVec(temp_file)
  expect_equal(dim(h5_final), dims)
  close(h5_final)
})