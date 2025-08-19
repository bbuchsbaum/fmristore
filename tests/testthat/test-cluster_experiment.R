# tests/testthat/test-cluster_experiment.R

library(testthat)
library(hdf5r)
library(neuroim2)
library(fmristore) # Make sure package functions are available



test_that("writer + reader round-trip", {
  tf <- tempfile(fileext = ".h5")
  on.exit(unlink(tf), add = TRUE)

  msp <- NeuroSpace(c(4, 4, 2), c(1, 1, 1))
  # Create dummy mask and clusters
  mask <- LogicalNeuroVol(array(c(rep(FALSE, 16 * 1), rep(TRUE, 16 * 1)), dim(msp)), msp)
  expect_equal(sum(mask), 16)
  set.seed(1)
  clus <- ClusteredNeuroVol(mask, sample(1:3, sum(mask), TRUE))

  # Generate some dummy data
  # Need nVoxels = sum(mask) = 16
  # cluster 1: ? voxels, cluster 2: ? voxels, cluster 3: ? voxels
  cluster_counts <- table(clus@clusters)
  `%||%` <- fmristore:::`%||%`
  n_vox_c1 <- cluster_counts[["1"]] %||% 0
  n_vox_c2 <- cluster_counts[["2"]] %||% 0
  n_vox_c3 <- cluster_counts[["3"]] %||% 0
  total_n_vox <- n_vox_c1 + n_vox_c2 + n_vox_c3
  expect_equal(total_n_vox, sum(mask)) # Sanity check

  n_time1 <- 10
  n_time2 <- 12

  # Create full run data list (matrices [nVoxInClus, nTime])
  full_data_list <- list()
  if (n_vox_c1 > 0) full_data_list$cluster_1 <- matrix(rnorm(n_vox_c1 * n_time1), n_vox_c1, n_time1)
  if (n_vox_c2 > 0) full_data_list$cluster_2 <- matrix(rnorm(n_vox_c2 * n_time1), n_vox_c2, n_time1)
  if (n_vox_c3 > 0) full_data_list$cluster_3 <- matrix(rnorm(n_vox_c3 * n_time1), n_vox_c3, n_time1)

  # Create summary run data matrix [nTime, nClusters]
  summ_mat  <- matrix(rnorm(n_time2 * 3), n_time2, 3)

  # Create cluster metadata data.frame
  clus_meta_df <- data.frame(cluster_id = 1:3,
    desc = paste("Cluster", 1:3),
    size = c(n_vox_c1, n_vox_c2, n_vox_c3))

  # Create the runs_data list for the writer
  runs <- list(
    list(scan_name = "run1",
      type = "full",
      data = full_data_list,
      metadata = list(TR = 2.0, Task = "RestingState", SubjectID = "Sub01")),
    list(scan_name = "run2",
      type = "summary",
      data = summ_mat,
      metadata = list(TR = 2.5, Task = "FingerTapping", SubjectID = "Sub01"))
  )

  # Write the file
  write_parcellated_experiment_h5(tf, mask, clus, runs,
    cluster_metadata = clus_meta_df, overwrite = TRUE, verbose = FALSE)

  # --- Test Reading ---
  # Can we open it?
  exp <- H5ParcellatedMultiScan(tf, keep_handle_open = TRUE)
  on.exit(try(exp$h5file$close_all(), silent = TRUE), add = TRUE) # Ensure cleanup

  expect_s4_class(exp, "H5ParcellatedMultiScan")
  expect_equal(n_scans(exp), 2)
  expect_equal(scan_names(exp), c("run1", "run2"))

  # Check run types
  expect_s4_class(exp@runs[[1]], "H5ParcellatedScan")
  expect_s4_class(exp@runs[[2]], "H5ParcellatedScanSummary")

  # Check loaded scan metadata
  expect_equal(exp@scan_metadata$run1$TR, 2.0)
  expect_equal(exp@scan_metadata$run1$Task, "RestingState")
  expect_equal(exp@scan_metadata$run1$SubjectID, "Sub01")
  expect_equal(exp@scan_metadata$run2$TR, 2.5)
  expect_equal(exp@scan_metadata$run2$Task, "FingerTapping")
  expect_equal(exp@scan_metadata$run2$SubjectID, "Sub01")

  # Check loaded cluster metadata
  expect_s3_class(exp@cluster_metadata, "data.frame")
  # Order might change on read/write depending on dataset listing order, so test carefully
  expect_setequal(names(exp@cluster_metadata), names(clus_meta_df))
  expect_equal(nrow(exp@cluster_metadata), nrow(clus_meta_df))
  # Check content after sorting by cluster_id to ensure alignment
  clus_meta_read_sorted <- exp@cluster_metadata[order(exp@cluster_metadata$cluster_id), ]
  clus_meta_orig_sorted <- clus_meta_df[order(clus_meta_df$cluster_id), ]
  expect_equal(clus_meta_read_sorted$cluster_id, clus_meta_orig_sorted$cluster_id)
  expect_equal(clus_meta_read_sorted$desc, clus_meta_orig_sorted$desc)
  expect_equal(clus_meta_read_sorted$size, clus_meta_orig_sorted$size)

  # Check shared components
  expect_true(identical(space(mask(exp)), space(mask)))
  expect_true(identical(mask(exp)@.Data, mask@.Data)) # Check data matches
  expect_true(identical(space(clusters(exp)), space(clus)))
  expect_true(identical(clusters(exp)@clusters, clus@clusters)) # Check data matches

  # Check handles and linkage within loaded object
  expect_true(identical(mask(exp), mask(exp@runs[[1]])))
  expect_true(identical(mask(exp), mask(exp@runs[[2]])))
  expect_true(identical(clusters(exp), clusters(exp@runs[[1]])))
  # Summary run might have NULL clusters slot if not explicitly loaded, check object identity if non-NULL
  if (!is.null(clusters(exp@runs[[2]]))) {
    expect_true(identical(clusters(exp), clusters(exp@runs[[2]])))
  }
  expect_true(identical(exp@runs[[1]]@obj, exp@runs[[1]]@obj))
  expect_true(identical(exp@runs[[1]]@obj, exp@runs[[2]]@obj))

  # --- Test Data Concatenation ---
  # Voxel concatenation (series_concat) - should match the input full data
  vox_indices <- 1:total_n_vox # Get all voxels in mask order

  # Reconstruct the expected full matrix [nTime, nVoxelsInMask]
  expected_full_t_vox <- matrix(NA_real_, nrow = n_time1, ncol = total_n_vox)
  mask_indices_global <- which(as.array(mask))

  for (cid_str in names(cluster_counts)) {
    cid <- as.integer(cid_str)
    current_run_data <- full_data_list[[paste0("cluster_", cid)]]
    if (!is.null(current_run_data)) {
      # Find which mask indices belong to this cluster
      vox_in_cluster_mask_indices <- which(clus@clusters == cid)
      # These indices correspond to columns in the final concatenated matrix
      expect_equal(nrow(current_run_data), length(vox_in_cluster_mask_indices)) # sanity check
      expected_full_t_vox[, vox_in_cluster_mask_indices] <- t(current_run_data)
    }
  }

  # Call series_concat for the first run (which is the full run)
  concatenated_voxels <- series_concat(exp, vox_indices, run_indices = 1)
  expect_equal(nrow(concatenated_voxels), n_time1)
  expect_equal(ncol(concatenated_voxels), total_n_vox)
  expect_equal(concatenated_voxels, expected_full_t_vox, tolerance = 1e-6)

  # Matrix concatenation (matrix_concat) - should match the input summary data
  # Call matrix_concat for the second run (which is the summary run)
  concatenated_summary <- matrix_concat(exp, run_indices = 2)
  expect_equal(nrow(concatenated_summary), n_time2)
  expect_equal(ncol(concatenated_summary), 3)
  expect_equal(concatenated_summary, summ_mat, tolerance = 1e-5, ignore_attr = TRUE)

  # Test trying to concat wrong types
  expect_error(series_concat(exp, vox_indices, run_indices = 2), ".*not an H5ParcellatedScan object.*")
  expect_error(matrix_concat(exp, run_indices = 1), ".*not an H5ParcellatedScanSummary object.*")
})

# context("H5ParcellatedMultiScan Constructor Validation")

test_that("writer errors if mask or clusters are NULL", {
  # Create minimal valid mask/clusters
  msp <- NeuroSpace(c(2, 2, 2), c(1, 1, 1))
  mask_val <- LogicalNeuroVol(array(TRUE, dim = c(2, 2, 2)), msp)
  clus_val <- ClusteredNeuroVol(mask_val, 1:8)
  runs_val <- list(list(scan_name = "s1", type = "summary", data = matrix(1:10, 5, 2)))
  tf <- tempfile(fileext = ".h5")
  on.exit(unlink(tf), add = TRUE)

  expect_error(write_parcellated_experiment_h5(tf, mask = NULL, clusters = clus_val, runs_data = runs_val), ".*must be.*LogicalNeuroVol.*")
  expect_error(write_parcellated_experiment_h5(tf, mask = mask_val, clusters = NULL, runs_data = runs_val), ".*must be.*ClusteredNeuroVol.*")
})

test_that("reader default path (mask=NULL, clusters=NULL) works", {
  tf <- tempfile(fileext = ".h5")
  on.exit(unlink(tf), add = TRUE)
  msp <- NeuroSpace(c(2, 2, 1), c(1, 1, 1))
  mask <- LogicalNeuroVol(array(c(FALSE, TRUE, TRUE, FALSE), dim = c(2, 2, 1)), msp)
  clus <- ClusteredNeuroVol(mask, 1:2)
  runs <- list(list(scan_name = "run1", type = "summary", data = matrix(1:10, 5, 2)))
  write_parcellated_experiment_h5(tf, mask, clus, runs, overwrite = TRUE, verbose = FALSE)

  # Should load successfully with mask=NULL, clusters=NULL
  exp <- H5ParcellatedMultiScan(tf)
  expect_s4_class(exp, "H5ParcellatedMultiScan")
  expect_s4_class(mask(exp), "LogicalNeuroVol")
  expect_s4_class(clusters(exp), "ClusteredNeuroVol")
  expect_equal(sum(mask(exp)), sum(mask)) # Check if loaded mask is correct
  # Check space consistency implicitly via identical check below
  expect_true(identical(space(mask(exp)), space(mask)))
})

test_that("reader validation works for provided mask and clusters", {
  tf <- tempfile(fileext = ".h5")
  on.exit(unlink(tf), add = TRUE)
  msp <- NeuroSpace(c(2, 2, 2), c(1, 1, 1))
  mask_orig <- LogicalNeuroVol(array(c(rep(FALSE, 4), rep(TRUE, 4)), dim = c(2, 2, 2)), msp)
  # Ensure we have both clusters represented
  cluster_ids <- c(1, 1, 2, 2)  # Guarantees both clusters are present
  clus_orig <- ClusteredNeuroVol(mask_orig, cluster_ids)
  runs <- list(list(scan_name = "run1", type = "summary", data = matrix(1:10, 5, 2)))
  # Write cluster metadata as well
  clus_meta_df <- data.frame(cluster_id = 1:2, desc = c("One", "Two"))
  write_parcellated_experiment_h5(tf, mask_orig, clus_orig, runs, cluster_metadata = clus_meta_df, overwrite = TRUE, verbose = FALSE)

  # 1. Provide valid mask and clusters - should succeed
  expect_message(exp_valid <- H5ParcellatedMultiScan(tf, mask = mask_orig, clusters = clus_orig), ".*validation successful.*", all = TRUE)
  expect_s4_class(exp_valid, "H5ParcellatedMultiScan")
  expect_true(identical(mask(exp_valid), mask_orig)) # Checks space and data
  expect_true(identical(clusters(exp_valid), clus_orig)) # Checks space and data

  # 2. Provide mask with different space - currently only warns, doesn't fail
  mask_bad_space <- LogicalNeuroVol(mask_orig@.Data, NeuroSpace(dim = c(2, 2, 2), spacing = c(2, 2, 2)))
  expect_warning(H5ParcellatedMultiScan(tf, mask = mask_bad_space), ".*NeuroSpace object does not match.*")

  # 3. Provide mask with different voxel pattern - should fail
  mask_bad_pattern <- mask_orig
  mask_bad_pattern[1, 1, 1] <- TRUE # Original was FALSE, this adds an extra TRUE voxel
  # This will fail at the count check before pattern check
  expect_error(H5ParcellatedMultiScan(tf, mask = mask_bad_pattern), ".*TRUE voxel count.*does not match.*")

  # 4. Provide mask with wrong sum (caught by length check against cluster_map)
  mask_bad_sum <- LogicalNeuroVol(array(TRUE, dim = c(2, 2, 2)), NeuroSpace(c(2, 2, 2), c(1, 1, 1))) # sum=8, but cluster_map length=4
  expect_error(H5ParcellatedMultiScan(tf, mask = mask_bad_sum), ".*TRUE voxel count.*does not match length.*cluster_map.*")

  # 5. Provide clusters with different space - should fail
  # Create a ClusteredNeuroVol with the wrong space mask
  clus_bad_space <- ClusteredNeuroVol(mask_bad_space, clus_orig@clusters)
  # Need to pass the original mask here, error happens when comparing clusters space to mask space
  expect_error(H5ParcellatedMultiScan(tf, mask = mask_orig, clusters = clus_bad_space), ".*clusters object.*NeuroSpace does not match the mask.*")

  # 6. Provide clusters with wrong length - should fail
  # Can't create ClusteredNeuroVol with wrong length - it would fail at construction
  # So we'll skip this test as it's not a valid test case
  # Test removed - invalid test case as ClusteredNeuroVol constructor would fail

})

test_that("summary_preference defaults based on /scans attribute", {
  tf1 <- tempfile(fileext = ".h5")
  tf2 <- tempfile(fileext = ".h5")
  on.exit(
    {
      unlink(tf1)
      unlink(tf2)
    },
    add = TRUE)

  mask <- fmristore:::create_minimal_LogicalNeuroVol(dims = c(2, 2, 2))
  clus <- fmristore:::create_minimal_ClusteredNeuroVol(mask_vol = mask, num_clusters = 2)

  runs_sum <- list(list(scan_name = "s1", type = "summary", data = matrix(1:4, 2, 2)))
  write_parcellated_experiment_h5(tf1, mask, clus, runs_sum, overwrite = TRUE, verbose = FALSE)

  exp_sum <- H5ParcellatedMultiScan(tf1)
  expect_s4_class(exp_sum@runs[[1]], "H5ParcellatedScanSummary")

  runs_full <- list(list(scan_name = "f1", type = "full", data = list(cluster_1 = matrix(1:4, 2, 2))))
  write_parcellated_experiment_h5(tf2, mask, clus, runs_full, overwrite = TRUE, verbose = FALSE)

  exp_full <- H5ParcellatedMultiScan(tf2)
  expect_s4_class(exp_full@runs[[1]], "H5ParcellatedScan")
})

test_that("warning for mismatched summary_only attribute", {
  tf <- tempfile(fileext = ".h5")
  on.exit(unlink(tf), add = TRUE)

  mask <- fmristore:::create_minimal_LogicalNeuroVol(dims = c(2, 2, 2))
  clus <- fmristore:::create_minimal_ClusteredNeuroVol(mask_vol = mask, num_clusters = 2)

  runs_full <- list(list(scan_name = "f1", type = "full", data = list(cluster_1 = matrix(1:4, 2, 2))))
  write_parcellated_experiment_h5(tf, mask, clus, runs_full, overwrite = TRUE, verbose = FALSE)

  h5f <- hdf5r::H5File$new(tf, mode = "r+")
  h5attr(h5f[["scans"]], "summary_only") <- TRUE
  h5f$close_all()

  expect_warning(exp_warn <- H5ParcellatedMultiScan(tf), "/scans@summary_only")
  expect_s4_class(exp_warn@runs[[1]], "H5ParcellatedScan")
})



test_that("dataset cluster_meta is read as data.frame", {
  tf <- tempfile(fileext = ".h5")
  on.exit(unlink(tf), add = TRUE)

  msp <- NeuroSpace(c(2, 2, 1), c(1, 1, 1))
  mask <- LogicalNeuroVol(array(TRUE, dim = c(2, 2, 1)), msp)
  clus <- ClusteredNeuroVol(mask, c(1L, 2L, 1L, 2L))

  full_data_list <- list(
    cluster_1 = matrix(rnorm(2 * 3), 2, 3),
    cluster_2 = matrix(rnorm(2 * 3), 2, 3)
  )
  runs <- list(list(scan_name = "run1", type = "full", data = full_data_list))
  meta_df <- data.frame(cluster_id = c(1L, 2L), desc = c("A", "B"))

  write_parcellated_experiment_h5(tf, mask, clus, runs,
    cluster_metadata = meta_df,
    overwrite = TRUE, verbose = FALSE)

  h5f <- H5File$new(tf, mode = "r+")
  meta_grp <- h5f[["/clusters/cluster_meta"]]
  meta_list <- lapply(names(meta_grp), function(nm) meta_grp[[nm]][])
  names(meta_list) <- names(meta_grp)
  df <- as.data.frame(meta_list)
  h5f$link_delete("/clusters/cluster_meta")
  h5f$create_dataset("/clusters/cluster_meta", df)
  h5f$close_all()

  exp <- H5ParcellatedMultiScan(tf)
  expect_s3_class(exp@cluster_metadata, "data.frame")
  expect_setequal(names(exp@cluster_metadata), names(meta_df))
  expect_equal(nrow(exp@cluster_metadata), nrow(meta_df))
})
