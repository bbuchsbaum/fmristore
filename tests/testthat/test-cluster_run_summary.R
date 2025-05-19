# test-cluster_run_summary.R

library(testthat)
library(hdf5r)
library(neuroim2)
# Assuming fmristore classes/methods are loaded via NAMESPACE or devtools::load_all()

# Helper function to create a dummy HDF5 file with SUMMARY data for testing
create_dummy_clustered_summary_h5 <- function(filepath,
                                              dims = c(4, 4, 2), # x,y,z
                                              n_time = 10,
                                              scan_name = "scan_summary1",
                                              cluster_ids = 1:3,
                                              summary_dset_name = "summary_data") {
  
  n_clusters <- length(cluster_ids)
  
  # Create basic mask (needed for context, but data isn't voxel-wise)
  mask_vol <- LogicalNeuroVol(array(TRUE, dim = dims), NeuroSpace(dims))

  # Create cluster map (clusters arg can be NULL for make_run_summary, 
  # but let's create one for potential consistency checks if needed later)
  # This cluster map structure might not perfectly align with the summary data,
  # as the summary data is just [nTime, nClusters].
  clus_map_arr <- array(rep(cluster_ids, length.out = prod(dims)), dim = dims)
  clusters_vol <- ClusteredNeuroVol(mask_vol, clusters = clus_map_arr[mask_vol])

  # Create predictable summary data [nTime, nClusters]
  summary_mat <- matrix(0, nrow = n_time, ncol = n_clusters)
  colnames(summary_mat) <- paste0("Cluster", cluster_ids)
  for (cl_idx in 1:n_clusters) {
      cid <- cluster_ids[cl_idx]
      # Data: cluster_id * 100 + timepoint
      summary_mat[, cl_idx] <- cid * 100 + 1:n_time
  }

  # Create HDF5 file
  file <- H5File$new(filepath, mode = "w")
  

  # Create structure /scans/<scan_name>/clusters_summary/
  scans_top_grp <- file$create_group("scans")
  scan_grp <- scans_top_grp$create_group(scan_name)
  summary_grp <- scan_grp$create_group("clusters_summary")

  # Write the summary dataset
  summary_grp[[summary_dset_name]] <- summary_mat
  
  # Optionally add cluster_names/ids as attributes if format supports it
  # summary_grp$create_attr("cluster_names", dtype=H5T_STRING, space=H5S_SCALAR) $write(paste0("Cluster", cluster_ids))

  file$close_all()
  
  return(list(filepath = filepath, 
              mask = mask_vol, 
              clusters = clusters_vol, # May be NULL in actual use case
              dims = dims, 
              n_time = n_time, 
              scan_name = scan_name,
              summary_data = summary_mat, 
              cluster_ids = cluster_ids,
              cluster_names = colnames(summary_mat),
              summary_dset_name = summary_dset_name
              ))
}



test_that("H5ClusteredRunSummary constructor works with file path", {
  setup_info <- create_dummy_clustered_summary_h5(tempfile(fileext = ".h5"))
  on.exit(unlink(setup_info$filepath), add = TRUE)

  # Use new constructor
  run_summary <- H5ClusteredRunSummary(file = setup_info$filepath,
                                       scan_name = setup_info$scan_name,
                                       mask = setup_info$mask,
                                       clusters = setup_info$clusters,
                                       cluster_names = setup_info$cluster_names,
                                       cluster_ids = setup_info$cluster_ids,
                                       summary_dset = setup_info$summary_dset_name
                                       )

  expect_s4_class(run_summary, "H5ClusteredRunSummary")
  expect_equal(run_summary@scan_name, setup_info$scan_name)
  expect_equal(run_summary@n_time, setup_info$n_time)
  expect_equal(run_summary@cluster_names, setup_info$cluster_names)
  expect_equal(run_summary@cluster_ids, setup_info$cluster_ids)
  expect_equal(run_summary@summary_dset, setup_info$summary_dset_name)
  expect_true(run_summary@obj$is_valid)
  
  # Close handle managed by object
  h5file(run_summary)$close_all()
})

test_that("H5ClusteredRunSummary constructor works with open H5File handle (via make_run_summary)", {
  setup_info <- create_dummy_clustered_summary_h5(tempfile(fileext = ".h5"))
  on.exit(unlink(setup_info$filepath), add = TRUE)

  h5f <- H5File$new(setup_info$filepath, mode = "r")
  on.exit(if(h5f$is_valid) h5f$close_all(), add = TRUE)

  # Keep using make_run_summary for this test, expecting a deprecation warning
  expect_warning(
    run_summary_open <- make_run_summary(file_source = h5f,
                                       scan_name = setup_info$scan_name,
                                       mask = setup_info$mask,
                                       clusters = setup_info$clusters,
                                       cluster_names = setup_info$cluster_names,
                                       cluster_ids = setup_info$cluster_ids,
                                       summary_dset = setup_info$summary_dset_name
                                       ),
    "deprecated"
  )

  expect_s4_class(run_summary_open, "H5ClusteredRunSummary")
  expect_true(run_summary_open@obj$is_valid)
  expect_equal(h5file(run_summary_open)$get_filename(), h5f$get_filename())
})

test_that("H5ClusteredRunSummary constructor errors work", {
  setup_info <- create_dummy_clustered_summary_h5(tempfile(fileext = ".h5"))
  on.exit(unlink(setup_info$filepath), add = TRUE)
  
  # Need an open handle to test internal errors
  h5f_err <- H5File$new(setup_info$filepath, mode = "r")
  on.exit(if(h5f_err$is_valid) h5f_err$close_all(), add = TRUE)

  # Use new constructor for error checks
  expect_error(H5ClusteredRunSummary(file = "nonexistent.h5", scan_name = "s1", mask = setup_info$mask, 
  clusters = setup_info$clusters), regexp="does not exist")
  # Test invalid scan_name within an existing file
  expect_error(H5ClusteredRunSummary(file = setup_info$filepath, 
                                     scan_name = "wrong_scan_name", mask = setup_info$mask, 
                                     clusters = setup_info$clusters), 
                                     "Summary dataset not found")
  # Test invalid summary_dset name
  expect_error(H5ClusteredRunSummary(file = setup_info$filepath, scan_name = setup_info$scan_name, mask = setup_info$mask, clusters = setup_info$clusters, summary_dset="wrong_name"), "Summary dataset not found")
})

test_that("H5ClusteredRunSummary cluster name/ID reconciliation works", {
  setup_info <- create_dummy_clustered_summary_h5(tempfile(fileext = ".h5"))
  on.exit(unlink(setup_info$filepath), add = TRUE)
  
  h5f <- H5File$new(setup_info$filepath, mode = "r")
  on.exit(if(h5f$is_valid) h5f$close_all(), add = TRUE)
  
  # Case 1: Provided names/IDs match dataset columns
  # Use new constructor (need file path)
  run_summary <- H5ClusteredRunSummary(file = setup_info$filepath,
                                       scan_name = setup_info$scan_name,
                                       mask = setup_info$mask,
                                       clusters = setup_info$clusters,
                                       cluster_names = paste0("clus_", 1:3), # Explicitly correct
                                       cluster_ids = 1:3, 
                                       summary_dset = setup_info$summary_dset_name)
  expect_equal(run_summary@cluster_names, paste0("clus_", 1:3))
  expect_equal(run_summary@cluster_ids, 1:3)
  h5file(run_summary)$close_all() # Close object handle

  # Case 2: Mismatched provided names (warning, reset to Col_X)
  # Use new constructor
  expect_warning(
      run_summary_badnames <- H5ClusteredRunSummary(file = setup_info$filepath,
                                                scan_name = setup_info$scan_name,
                                                mask = setup_info$mask,
                                                clusters = setup_info$clusters,
                                                cluster_names = c("A", "B"), # Mismatch cols (3)
                                                cluster_ids = 1:3, # Match cols
                                                summary_dset = setup_info$summary_dset_name),
      "Final number of cluster names.*does not match dataset columns"
  )
  expect_warning( # Second warning about resetting
      H5ClusteredRunSummary(file = setup_info$filepath,
                         scan_name = setup_info$scan_name,
                         mask = setup_info$mask,
                         clusters = setup_info$clusters,
                         cluster_names = c("A", "B"), # Mismatch cols (3)
                         cluster_ids = 1:3, # Match cols
                         summary_dset = setup_info$summary_dset_name),
      "Resetting names/IDs to Col_X/sequential"
  )
  expect_equal(run_summary_badnames@cluster_names, paste0("Col_", 1:3))
  expect_equal(run_summary_badnames@cluster_ids, 1:3)
  h5file(run_summary_badnames)$close_all() # Close object handle
  
  # Case 3: No names/IDs provided, derive from clusters object
  # Use new constructor
  run_summary_derive <- H5ClusteredRunSummary(file = setup_info$filepath,
                                            scan_name = setup_info$scan_name,
                                            mask = setup_info$mask,
                                            clusters = setup_info$clusters,
                                            cluster_names = character(), # Explicitly empty
                                            cluster_ids = integer(),   # Explicitly empty
                                            summary_dset = setup_info$summary_dset_name)
  expect_equal(run_summary_derive@cluster_names, paste0("Clus_", 1:3)) # Derived from unique(clusters@clusters)
  expect_equal(run_summary_derive@cluster_ids, 1:3)
  h5file(run_summary_derive)$close_all()

  # Case 4: No names/IDs/clusters provided, derive from dataset cols (Col_X)
  # Use new constructor
  expect_warning(
      run_summary_defaults <- H5ClusteredRunSummary(file = setup_info$filepath,
                                                scan_name = setup_info$scan_name,
                                                mask = setup_info$mask,
                                                clusters = NULL, # No clusters object
                                                cluster_names = character(), 
                                                cluster_ids = integer(),   
                                                summary_dset = setup_info$summary_dset_name),
      "generated default column names \\(Col_X\\)"
  )
  expect_equal(run_summary_defaults@cluster_names, paste0("Col_", 1:3))
  expect_equal(run_summary_defaults@cluster_ids, 1:3)
  h5file(run_summary_defaults)$close_all()
  
})

test_that("as.matrix method for H5ClusteredRunSummary works", {
  setup_info <- create_dummy_clustered_summary_h5(tempfile(fileext = ".h5"))
  on.exit(unlink(setup_info$filepath), add = TRUE)
  
  h5f <- H5File$new(setup_info$filepath, mode = "r")
  on.exit(h5f$close_all(), add = TRUE)
  
  run_summary <- make_run_summary(file_source = h5f,
                                  scan_name = setup_info$scan_name,
                                  mask = setup_info$mask,
                                  clusters = setup_info$clusters,
                                  cluster_names = setup_info$cluster_names,
                                  cluster_ids = setup_info$cluster_ids,
                                  summary_dset = setup_info$summary_dset_name
                                  )
                                  
  # Call as.matrix
  mat <- as.matrix(run_summary)
  
  # Check dimensions
  expect_equal(nrow(mat), setup_info$n_time)
  expect_equal(ncol(mat), length(setup_info$cluster_ids))
  
  # Check column names
  expect_equal(colnames(mat), setup_info$cluster_names)
  
  # Check content
  expect_equal(mat, setup_info$summary_data)
  
  # Test case where cluster_names don't match columns (should warn)
   run_summary_badnames <- make_run_summary(file_source = h5f,
                                  scan_name = setup_info$scan_name,
                                  mask = setup_info$mask,
                                  clusters = setup_info$clusters,
                                  cluster_names = c("Wrong", "Names"), # Incorrect number
                                  cluster_ids = setup_info$cluster_ids
                                  )
  expect_warning(mat_badnames <- as.matrix(run_summary_badnames), "Length of cluster_names")
  expect_null(colnames(mat_badnames)) # Names should not be set
  expect_equal(mat_badnames, setup_info$summary_data, check.attributes = FALSE) # Content still matches
  
})

test_that("make_run_summary generates default names/ids and as.data.frame works", {
  setup_info <- create_dummy_clustered_summary_h5(tempfile(fileext = ".h5"))
  on.exit(unlink(setup_info$filepath), add = TRUE)
  
  h5f <- H5File$new(setup_info$filepath, mode = "r")
  on.exit(h5f$close_all(), add = TRUE)
  
  # Call constructor without names or ids
  run_summary_defaults <- NULL
  
    expect_warning(
       run_summary_defaults <- make_run_summary(file_source = h5f,
                                     scan_name = setup_info$scan_name,
                                     mask = setup_info$mask,
                                     clusters = NULL # Pass NULL to test default generation
                                     # Omit cluster_names and cluster_ids
                                     ),
      regexp="Generated default" # Warning for names
    )
  expect_s4_class(run_summary_defaults, "H5ClusteredRunSummary")
  n_clusters_expected <- ncol(setup_info$summary_data)
  expected_names <- paste0("Col_", seq_len(n_clusters_expected))
  expected_ids <- seq_len(n_clusters_expected)
  
  # Check if defaults were set in the object
  expect_equal(run_summary_defaults@cluster_names, expected_names)
  expect_equal(run_summary_defaults@cluster_ids, expected_ids)
  
  # Test as.data.frame
  df <- as.data.frame(run_summary_defaults)
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), setup_info$n_time)
  expect_equal(ncol(df), n_clusters_expected)
  expect_equal(names(df), expected_names) # Column names should be the generated defaults
  
  # Check content consistency with matrix form
  expect_equal(as.matrix(df), setup_info$summary_data, check.attributes = FALSE)
})


