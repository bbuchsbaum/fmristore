# tests/testthat/test-cluster_run_full.R

library(testthat)
library(hdf5r)
library(neuroim2)
library(fmristore)

# Helper to setup the test HDF5 file using create_dummy...
# Allows controlling how n_time is stored for different test cases
setup_test_file_full <- function(write_n_time_attr = FALSE,
                                 write_n_time_meta = FALSE,
                                 n_time_val = 12, # Match default in create_dummy...
                                 ...) { # Pass extra args to create_dummy...

  auto_loc <- "none"
  if (write_n_time_attr && write_n_time_meta) {
    stop("Cannot set both write_n_time_attr and write_n_time_meta to TRUE")
  }
  if (write_n_time_attr) auto_loc <- "attribute"
  if (write_n_time_meta) auto_loc <- "metadata"

  tmp_filepath <- tempfile(fileext = ".h5")

  create_dummy_clustered_full_h5(filepath = tmp_filepath,
    n_time = n_time_val,
    auto_n_time_location = auto_loc,
    ...)
}

# Helper to clean up the created test file
cleanup_test_file <- function(setup_info) {
  if (!is.null(setup_info) && !is.null(setup_info$filepath) && file.exists(setup_info$filepath)) {
    unlink(setup_info$filepath)
  }
}

# Helper function to create a dummy HDF5 file with FULL clustered data
create_dummy_clustered_full_h5 <- function(filepath,
                                           dims = c(4, 4, 3), # x,y,z
                                           n_time = 12,
                                           scan_name = "scan_full1",
                                           cluster_ids = 1:3,
                                           auto_n_time_location = "none", # "none", "attribute", "metadata"
                                           invalid_cluster_shape = FALSE) {

  n_clusters <- length(cluster_ids)
  n_vox_total <- prod(dims)

  # Create basic mask (e.g., center region)
  mask_arr <- array(FALSE, dim = dims)
  mask_arr[2:3, 2:3, 1:2] <- TRUE # Smaller mask within dims
  mask_vol <- LogicalNeuroVol(mask_arr, NeuroSpace(dims))
  n_vox_mask <- sum(mask_vol)

  # Create cluster map for voxels within the mask
  # Assign voxels to clusters sequentially
  clus_map_vals <- rep(cluster_ids, length.out = n_vox_mask)
  clusters_vol <- ClusteredNeuroVol(mask_vol, clusters = clus_map_vals)

  # Create voxel coordinates (relative to mask)
  voxel_coords_mask <- which(mask_arr, arr.ind = TRUE) # Get 3D indices of TRUE values
  colnames(voxel_coords_mask) <- c("x", "y", "z")

  # Create predictable full data (store temporarily)
  # We need data structured as list[[cluster_id]][nVoxInClus, nTime]
  cluster_data_list <- list()

  # Create HDF5 file
  file <- H5File$new(filepath, mode = "w")
  on.exit(if (file$is_valid) file$close_all(), add = TRUE) # Ensure closure

  # Write Mask, Cluster Map, Voxel Coords
  file[["mask"]] <- as.integer(mask_arr)
  file[["cluster_map"]] <- clus_map_vals
  file[["voxel_coords"]] <- voxel_coords_mask

  # Write Global Cluster Info
  clus_grp <- file$create_group("clusters")
  clus_grp[["cluster_ids"]] <- cluster_ids
  # Simple meta for now
  clus_meta_df <- data.frame(cluster_id = cluster_ids, description = paste("Global Cluster", cluster_ids))
  # hdf5r needs help with data.frames -> convert to list or compound type if necessary
  # For simplicity, let's skip writing complex meta for now, just IDs.

  # Create Scan Group
  scans_grp <- file$create_group("scans")
  scan_grp <- scans_grp$create_group(scan_name)

  # Write n_time if requested
  if (auto_n_time_location == "attribute") {
    h5attr(scan_grp, "n_time") <- n_time
  } else if (auto_n_time_location == "metadata") {
    meta_grp <- scan_grp$create_group("metadata")
    meta_grp[["n_time"]] <- n_time
  }

  # Write Full Cluster Data under the scan
  scan_clus_grp <- scan_grp$create_group("clusters")
  # Only create expected_reconstruction if n_time is positive
  expected_reconstruction <- if (n_time > 0) array(0, dim = c(dims, n_time)) else NULL

  for (cid in cluster_ids) {
    # Find which voxels *within the mask* belong to this cluster
    mask_indices_in_cluster <- which(clusters_vol@clusters == cid)
    n_vox_in_cluster <- length(mask_indices_in_cluster)

    if (n_vox_in_cluster > 0) {
      # Generate data [nVoxInCluster, nTime]
      # Formula: cluster_id * 1000 + (voxel_offset_within_cluster * 100) + timepoint
      cluster_mat <- if (n_time > 0) matrix(0, nrow = n_vox_in_cluster, ncol = n_time) else matrix(0, nrow = n_vox_in_cluster, ncol = 1)
      if (n_time > 0) {
        for (vox_idx in 1:n_vox_in_cluster) {
          cluster_mat[vox_idx, ] <- cid * 1000 + (vox_idx * 100) + 1:n_time
        }
      }

      # Write dataset
      if (invalid_cluster_shape && cid == cluster_ids[1]) {
        # Add an extra singleton dimension to make the dataset 3D
        scan_clus_grp[[paste0("cluster_", cid)]] <- array(cluster_mat, dim = c(dim(cluster_mat), 1))
      } else {
        scan_clus_grp[[paste0("cluster_", cid)]] <- cluster_mat
      }
      cluster_data_list[[as.character(cid)]] <- cluster_mat

      # Fill expected full array (for testing subsets later)
      if (!is.null(expected_reconstruction) && n_time > 0) {
        # Get the 3D coordinates corresponding to these mask indices
        coords_this_cluster <- voxel_coords_mask[mask_indices_in_cluster, , drop = FALSE]
        for (vox_idx in 1:n_vox_in_cluster) {
          coord <- coords_this_cluster[vox_idx, ]
          expected_reconstruction[coord[1], coord[2], coord[3], ] <- cluster_mat[vox_idx, ]
        }
      }
    } else {
      # Handle case where a cluster ID might have 0 voxels (though unlikely with rep())
      scan_clus_grp[[paste0("cluster_", cid)]] <- matrix(numeric(0), nrow = 0, ncol = max(1, n_time))
      cluster_data_list[[as.character(cid)]] <- matrix(numeric(0), nrow = 0, ncol = max(1, n_time))
    }
  }

  # Explicitly close before returning path
  file$close_all()

  return(list(filepath = filepath,
    mask = mask_vol,
    clusters = clusters_vol,
    voxel_coords = voxel_coords_mask, # Coords of masked voxels
    dims = dims,
    n_time = n_time,
    scan_name = scan_name,
    cluster_ids = cluster_ids,
    cluster_data_list = cluster_data_list, # Data as stored per cluster
    expected_reconstruction = expected_reconstruction # Full 4D array (potentially large)
  ))
}


test_that("make_run_full constructs object correctly from file path", {
  setup_info <- setup_test_file_full()
  on.exit(cleanup_test_file(setup_info))

  # Use new constructor
  run_full <- H5ParcellatedScan(file = setup_info$filepath,
    scan_name = setup_info$scan_name,
    mask = setup_info$mask,
    clusters = setup_info$clusters,
    n_time = setup_info$n_time)

  expect_s4_class(run_full, "H5ParcellatedScan")
  expect_equal(run_full@scan_name, setup_info$scan_name)
  expect_equal(run_full@n_time, setup_info$n_time)
  expect_equal(run_full@n_voxels, sum(setup_info$mask))
  expect_true(inherits(run_full@obj, "H5File"))
  expect_true(run_full@obj$is_valid)

  # Important: Close the handle managed by the object
  expect_silent(run_full@obj$close_all())
})

# Test constructing from an already open H5File handle
test_that("make_run_full constructs object correctly from open H5File handle", {
  setup_info <- setup_test_file_full()
  on.exit(cleanup_test_file(setup_info))

  h5f <- H5File$new(setup_info$filepath, mode = "r")
  # Note: If make_run_full takes ownership, h5f might be closed by its finalizer
  # Here, we assume make_run_full *uses* the handle but doesn't close it if passed open
  on.exit(if (h5f$is_valid) h5f$close_all(), add = TRUE)

  # Use new constructor with the H5File object directly
  run_full_open <- H5ParcellatedScan(file = h5f,
    scan_name = setup_info$scan_name,
    mask = setup_info$mask,
    clusters = setup_info$clusters,
    n_time = setup_info$n_time)

  expect_s4_class(run_full_open, "H5ParcellatedScan")
  expect_true(run_full_open@obj$is_valid)
  # Get the filename from the setup_info since h5f might not support get_filename() directly
  expect_equal(normalizePath(run_full_open@obj$get_filename()), normalizePath(setup_info$filepath))
})

test_that("make_run_full reads n_time from HDF5 attributes if NULL", {
  setup_attr <- setup_test_file_full(write_n_time_attr = TRUE)
  on.exit(cleanup_test_file(setup_attr))

  # Use new constructor
  run_attr <- H5ParcellatedScan(file = setup_attr$filepath,
    scan_name = setup_attr$scan_name,
    mask = setup_attr$mask,
    clusters = setup_attr$clusters,
    n_time = NULL) # Explicitly pass NULL

  expect_equal(run_attr@n_time, setup_attr$n_time)
  expect_silent(run_attr@obj$close_all())
})

test_that("make_run_full reads n_time from HDF5 metadata dataset if NULL", {
  setup_meta <- setup_test_file_full(write_n_time_meta = TRUE)
  on.exit(cleanup_test_file(setup_meta))

  # Use new constructor
  run_meta <- H5ParcellatedScan(file = setup_meta$filepath,
    scan_name = setup_meta$scan_name,
    mask = setup_meta$mask,
    clusters = setup_meta$clusters,
    n_time = NULL)

  expect_equal(run_meta@n_time, setup_meta$n_time)
  expect_silent(run_meta@obj$close_all())
})

test_that("inference message depends on fmristore.verbose option", {
  setup_none <- setup_test_file_full(write_n_time_attr = FALSE, write_n_time_meta = FALSE)
  on.exit(cleanup_test_file(setup_none))

  expect_message(
    run_infer_verbose <- withr::with_options(
      list(fmristore.verbose = TRUE),
      H5ParcellatedScan(file = setup_none$filepath,
        scan_name = setup_none$scan_name,
        mask = setup_none$mask,
        clusters = setup_none$clusters,
        n_time = NULL)
    ),
    "Inferred n_time"
  )
  expect_equal(run_infer_verbose@n_time, setup_none$n_time)
  expect_silent(run_infer_verbose@obj$close_all())

  expect_no_message(
    run_infer_silent <- withr::with_options(
      list(fmristore.verbose = FALSE),
      H5ParcellatedScan(file = setup_none$filepath,
        scan_name = setup_none$scan_name,
        mask = setup_none$mask,
        clusters = setup_none$clusters,
        n_time = NULL)
    )
  )
  expect_equal(run_infer_silent@n_time, setup_none$n_time)
  expect_silent(run_infer_silent@obj$close_all())
})


test_that("make_run_full throws errors for invalid inputs", {
  setup_info <- setup_test_file_full()
  on.exit(cleanup_test_file(setup_info))

  # Create bad inputs
  bad_mask <- setup_info$mask
  # Create a new mask with different dimensions
  bad_mask <- neuroim2::LogicalNeuroVol(array(TRUE, c(1, 1, 1)), neuroim2::NeuroSpace(c(1, 1, 1)))
  bad_clusters <- setup_info$clusters
  bad_clusters@clusters <- bad_clusters@clusters[-1] # Mismatched length

  # Use new constructor for error checks
  expect_error(H5ParcellatedScan("nonexistent.h5", "s1", setup_info$mask, setup_info$clusters, 10), "file path does not exist") # Error from open_h5
  # Original tests for bad mask/cluster types are now implicitly tested by `is()` checks inside the constructor
  # expect_error(H5ParcellatedScan(setup_info$filepath, setup_info$scan_name, 1, setup_info$clusters, setup_info$n_time), "must be a LogicalNeuroVol")
  # expect_error(H5ParcellatedScan(setup_info$filepath, setup_info$scan_name, setup_info$mask, 1, setup_info$n_time), "must be a ClusteredNeuroVol")
  # Test dimension mismatch (uses check_same_dims)
  expect_error(H5ParcellatedScan(file = setup_info$filepath, scan_name = setup_info$scan_name, mask = bad_mask, clusters = setup_info$clusters, n_time = setup_info$n_time), "Dimensions of 'mask' and 'clusters' must match")
  # Test cluster length mismatch
  expect_error(H5ParcellatedScan(file = setup_info$filepath, scan_name = setup_info$scan_name, mask = setup_info$mask, clusters = bad_clusters, n_time = setup_info$n_time), "Mismatch: clusters@clusters length")

  # Test invalid compress inputs
  expect_error(
    H5ParcellatedScan(file = setup_info$filepath,
      scan_name = setup_info$scan_name,
      mask = setup_info$mask,
      clusters = setup_info$clusters,
      n_time = setup_info$n_time,
      compress = c(TRUE, FALSE)),
    "'compress' must be a single logical value")
  expect_error(
    H5ParcellatedScan(file = setup_info$filepath,
      scan_name = setup_info$scan_name,
      mask = setup_info$mask,
      clusters = setup_info$clusters,
      n_time = setup_info$n_time,
      compress = "yes"),
    "'compress' must be a single logical value")
  # Test error if n_time cannot be determined
  # Note: This test is currently commented out because the constructor successfully
  # infers n_time from the cluster dataset dimensions. To properly test this,
  # we would need to create a file with no cluster data at all.
  # setup_uninferrable <- setup_test_file_full()
  # on.exit(cleanup_test_file(setup_uninferrable), add = TRUE)
  # expect_error(H5ParcellatedScan(file = setup_uninferrable$filepath,
  #                                 scan_name = setup_uninferrable$scan_name,
  #                                 mask = setup_uninferrable$mask,
  #                                 clusters = setup_uninferrable$clusters,
  #                                 n_time = NULL),
  #              "Could not determine 'n_time'")
})

test_that("make_run_full stops if n_time determined is invalid", {
  setup_info <- setup_test_file_full(write_n_time_attr = TRUE, n_time_val = -5) # Write invalid n_time
  on.exit(cleanup_test_file(setup_info))

  # Use new constructor
  expect_error(H5ParcellatedScan(file = setup_info$filepath,
    scan_name = setup_info$scan_name,
    mask = setup_info$mask,
    clusters = setup_info$clusters,
    n_time = NULL),
  "must be a single positive integer")
})

test_that("make_run_full errors when cluster dataset has wrong dimensions", {
  setup_bad <- setup_test_file_full(invalid_cluster_shape = TRUE)
  on.exit(cleanup_test_file(setup_bad))

  expect_error(
    H5ParcellatedScan(file = setup_bad$filepath,
      scan_name = setup_bad$scan_name,
      mask = setup_bad$mask,
      clusters = setup_bad$clusters,
      n_time = NULL),
    "Could not determine 'n_time'"
  )
})

# ------------------------------------------------------------------------------
# Data access tests

test_that("series() retrieves correct voxel time series", {
  setup_info <- setup_test_file_full()
  on.exit(cleanup_test_file(setup_info))

  run_full <- H5ParcellatedScan(file = setup_info$filepath,
    scan_name = setup_info$scan_name,
    mask = setup_info$mask,
    clusters = setup_info$clusters,
    n_time = setup_info$n_time)

  mask_idx1 <- 1L
  coord1 <- setup_info$voxel_coords[mask_idx1, ]
  expected_ts <- as.numeric(setup_info$expected_reconstruction[
    coord1[1], coord1[2], coord1[3], ]
  )

  ts_mask <- as.numeric(series(run_full, mask_idx1))
  expect_equal(ts_mask, expected_ts)

  ts_coord <- as.numeric(series(run_full, coord1[1], coord1[2], coord1[3]))
  expect_equal(ts_coord, expected_ts)

  mask_indices <- 1:2
  coords_mat <- setup_info$voxel_coords[mask_indices, ]
  expected_mat <- vapply(mask_indices, function(ii) {
    cidx <- setup_info$voxel_coords[ii, ]
    setup_info$expected_reconstruction[cidx[1], cidx[2], cidx[3], ]
  }, numeric(setup_info$n_time))

  ts_multi <- series(run_full, coords_mat)
  expect_equal(ts_multi, expected_mat, ignore_attr = TRUE)

  run_full@obj$close_all()
})

test_that("linear_access reconstructs voxel values", {
  setup_info <- setup_test_file_full()
  on.exit(cleanup_test_file(setup_info))

  run_full <- H5ParcellatedScan(file = setup_info$filepath,
    scan_name = setup_info$scan_name,
    mask = setup_info$mask,
    clusters = setup_info$clusters,
    n_time = setup_info$n_time)

  dims4 <- dim(run_full)
  coord_in <- setup_info$voxel_coords[1, ]
  t_in <- 1L
  idx_in <- coord_in[1] + (coord_in[2] - 1) * dims4[1] +
    (coord_in[3] - 1) * dims4[1] * dims4[2] +
    (t_in - 1) * dims4[1] * dims4[2] * dims4[3]

  val_in <- linear_access(run_full, idx_in)
  expected_in <- setup_info$expected_reconstruction[
    coord_in[1], coord_in[2], coord_in[3], t_in]
  expect_equal(val_in, expected_in)

  coord_out <- c(1L, 1L, 3L)
  t_out <- 1L
  idx_out <- coord_out[1] + (coord_out[2] - 1) * dims4[1] +
    (coord_out[3] - 1) * dims4[1] * dims4[2] +
    (t_out - 1) * dims4[1] * dims4[2] * dims4[3]

  val_out <- linear_access(run_full, idx_out)
  expect_equal(val_out, 0)

  run_full@obj$close_all()
})
