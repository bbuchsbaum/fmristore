# Test suite for HDF5 file corruption and invalid format handling
# This addresses a critical gap in testing robustness against malformed files

library(fmristore)

test_that("H5ParcellatedMultiScan handles invalid cluster configurations", {
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
  clusters_grp[["cluster_4"]] <- matrix(rnorm(10 * 5), nrow = 10, ncol = 5)
  clusters_grp[["cluster_5"]] <- matrix(rnorm(8 * 5), nrow = 8, ncol = 5)
  clusters_grp[["cluster_6"]] <- matrix(rnorm(7 * 5), nrow = 7, ncol = 5)

  h5file$close()

  # Test: Should detect mismatch between cluster_map and actual cluster datasets
  expect_error(
    H5ParcellatedMultiScan(temp_file),
    info = "Should fail when cluster IDs in data don't match cluster_map"
  )

  # Test 2: Corrupt cluster metadata
  temp_file2 <- tempfile(fileext = ".h5")
  on.exit(unlink(temp_file2), add = TRUE)

  # Create valid initial structure
  fmristore:::create_minimal_h5_for_H5ParcellatedMultiScan(
    file_path = temp_file2,
    master_mask_dims = c(5L, 5L, 5L), # Fixed: use correct parameter name
    num_master_clusters = 3L, # Fixed: use correct parameter name
    n_time_run1 = 10L
  )

  # Create /clusters group and add invalid cluster_meta
  h5file <- hdf5r::H5File$new(temp_file2, mode = "r+")
  # Create /clusters group if it doesn't exist
  if (!h5file$exists("clusters")) {
    h5file$create_group("clusters")
  }
  clusters_grp <- h5file[["clusters"]]
  # Create cluster_meta as an H5Group with inconsistent dataset lengths
  # This should trigger the "inconsistent lengths" warning
  meta_grp <- clusters_grp$create_group("cluster_meta")
  meta_grp$create_dataset("id", robj = 1:3, dtype = hdf5r::h5types$H5T_NATIVE_INT)
  meta_grp$create_dataset("name", robj = c("A", "B"), dtype = hdf5r::H5T_STRING$new(size = 10)) # Wrong length!
  h5file$close()

  # Should fail with error about invalid cluster_metadata slot (list instead of data.frame)
  # This happens because the inconsistent lengths cause it to return a list
  expect_error(
    H5ParcellatedMultiScan(temp_file2),
    regexp = "invalid object for slot.*cluster_metadata.*got class.*list.*should be.*data.frame",
    info = "Should fail when cluster metadata has inconsistent lengths"
  )
})

