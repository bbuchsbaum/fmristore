library(testthat)
library(hdf5r)
library(neuroim2)
library(fmristore)

# Test .dataset_path method for H5ParcellatedScan

test_that(".dataset_path validates inputs and constructs path", {
  mask <- fmristore:::create_minimal_LogicalNeuroVol(dims = c(2, 2, 1))
  clus <- fmristore:::create_minimal_ClusteredNeuroVol(mask_vol = mask, num_clusters = 2L)
  h5f <- H5File$new(tempfile(fileext = ".h5"), mode = "w")
  on.exit(
    {
      fn <- h5f$get_filename()
      h5f$close_all()
      unlink(fn)
    },
    add = TRUE)

  obj <- h5f
  run <- new("H5ParcellatedScan",
    obj = obj,
    mask = mask,
    clusters = clus,
    n_voxels = sum(mask),
    scan_name = "runA",
    n_time = 5L,
    compress = FALSE)

  expect_equal(fmristore:::.dataset_path(run, 1L), "/scans/runA/clusters/cluster_1")
  expect_error(fmristore:::.dataset_path(run, 0L), "single positive integer")

  run_bad <- new("H5ParcellatedScan",
    obj = obj,
    mask = mask,
    clusters = clus,
    n_voxels = sum(mask),
    scan_name = "",
    n_time = 5L,
    compress = FALSE)
  expect_error(fmristore:::.dataset_path(run_bad, 1L), "Invalid 'scan_name'")
})

# Test internal .get_cluster_timeseries_by_mask_index helper

test_that(".get_cluster_timeseries_by_mask_index extracts data and checks bounds", {
  # Create mask with specific TRUE voxels to match cluster count
  mask <- fmristore:::create_minimal_LogicalNeuroVol(dims = c(2, 2, 1),
    true_voxels = list(c(1L, 1L, 1L), c(2L, 1L, 1L),
      c(1L, 2L, 1L), c(2L, 2L, 1L)))
  clus_vals <- c(1L, 2L, 1L, 2L)
  clus <- ClusteredNeuroVol(mask, clusters = clus_vals)
  n_time <- 3L

  tmp <- tempfile(fileext = ".h5")
  h5f <- H5File$new(tmp, mode = "w")
  on.exit(
    {
      h5f$close_all()
      unlink(tmp)
    },
    add = TRUE)
  # Create nested group structure
  scans_grp <- h5f$create_group("scans")
  runA_grp <- scans_grp$create_group("runA")
  grp <- runA_grp$create_group("clusters")
  grp[["cluster_1"]] <- matrix(1:6, nrow = 2, ncol = 3)
  grp[["cluster_2"]] <- matrix(11:16, nrow = 2, ncol = 3)

  run <- new("H5ParcellatedScan",
    obj = h5f,
    mask = mask,
    clusters = clus,
    n_voxels = sum(mask),
    scan_name = "runA",
    n_time = n_time,
    compress = FALSE)

  res <- fmristore:::.get_cluster_timeseries_by_mask_index(run, c(1, 3), c(1, 3), n_time)
  # mask index 1 -> cluster 1, position 1: times 1,3 -> [1,5]
  # mask index 3 -> cluster 1, position 2: times 1,3 -> [2,6]
  expect_equal(res, matrix(c(1, 2, 5, 6), nrow = 2, ncol = 2))

  res2 <- fmristore:::.get_cluster_timeseries_by_mask_index(run, c(2, 4), NULL, n_time)
  # mask index 2 -> cluster 2, position 1: all times -> [11,13,15]
  # mask index 4 -> cluster 2, position 2: all times -> [12,14,16]
  expect_equal(res2, matrix(c(11, 12, 13, 14, 15, 16), nrow = 2, ncol = 3))

  expect_error(fmristore:::.get_cluster_timeseries_by_mask_index(run, 5L, NULL, n_time),
    "out of the valid mask range")
  expect_error(fmristore:::.get_cluster_timeseries_by_mask_index(run, 1L, 4L, n_time),
    "valid time range")
})
