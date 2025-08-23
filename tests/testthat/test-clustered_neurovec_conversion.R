test_that("ClusteredNeuroVec conversion defaults to single scan", {
  skip_if_not_installed("neuroim2")
  
  library(neuroim2)
  library(fmristore)
  
  # Create a simple test ClusteredNeuroVec
  # First create a mask
  mask_array <- array(FALSE, dim = c(10, 10, 10))
  mask_array[3:8, 3:8, 3:8] <- TRUE
  mask <- neuroim2::LogicalNeuroVol(mask_array, neuroim2::NeuroSpace(c(10, 10, 10)))
  
  # Create cluster assignments (must match number of TRUE voxels in mask)
  n_voxels <- sum(mask_array)
  n_clusters <- 5
  cluster_ids <- sample(1:n_clusters, n_voxels, replace = TRUE)
  
  # Create ClusteredNeuroVol using the constructor function
  cvol <- neuroim2::ClusteredNeuroVol(mask = mask, clusters = cluster_ids)
  
  # Create time series data (T x K matrix)
  n_time <- 100
  ts_data <- matrix(rnorm(n_time * n_clusters), nrow = n_time, ncol = n_clusters)
  
  # Create ClusteredNeuroVec - pass matrix as first arg, cvol as second
  cnvec <- neuroim2::ClusteredNeuroVec(x = ts_data, cvol = cvol)
  
  # Test write_dataset method with default (single scan)
  temp_file <- tempfile(fileext = ".h5")
  on.exit(unlink(temp_file), add = TRUE)
  
  result <- write_dataset(cnvec, file = temp_file, scan_name = "test_scan")
  
  # Check that result is H5ParcellatedScanSummary (single scan)
  expect_s4_class(result, "H5ParcellatedScanSummary")
  
  # Check file exists
  expect_true(file.exists(temp_file))
  
  # Check scan name
  expect_equal(result@scan_name, "test_scan")
  
  # Check time series dimensions
  expect_equal(result@n_time, n_time)
  
  # Check mask dimensions match
  expect_equal(dim(result@mask), dim(mask))
  
  # Check clusters match
  expect_equal(dim(result@clusters), dim(cvol))
  
  # Clean up
  close(result)
})

test_that("ClusteredNeuroVec conversion to multi-scan works", {
  skip_if_not_installed("neuroim2")
  
  library(neuroim2)
  library(fmristore)
  
  # Create a simple test ClusteredNeuroVec
  mask_array <- array(FALSE, dim = c(8, 8, 8))
  mask_array[3:6, 3:6, 3:6] <- TRUE
  mask <- neuroim2::LogicalNeuroVol(mask_array, neuroim2::NeuroSpace(c(8, 8, 8)))
  
  n_voxels <- sum(mask_array)
  n_clusters <- 3
  cluster_ids <- sample(1:n_clusters, n_voxels, replace = TRUE)
  
  cvol <- neuroim2::ClusteredNeuroVol(mask = mask, clusters = cluster_ids)
  
  n_time <- 50
  ts_data <- matrix(rnorm(n_time * n_clusters), nrow = n_time, ncol = n_clusters)
  
  cnvec <- neuroim2::ClusteredNeuroVec(x = ts_data, cvol = cvol)
  
  # Test write_dataset with as_multiscan=TRUE
  temp_file <- tempfile(fileext = ".h5")
  on.exit(unlink(temp_file), add = TRUE)
  
  result <- write_dataset(cnvec, file = temp_file, scan_name = "multi_test", as_multiscan = TRUE)
  
  # Check that result is H5ParcellatedMultiScan
  expect_s4_class(result, "H5ParcellatedMultiScan")
  
  # Check file exists
  expect_true(file.exists(temp_file))
  
  # Check that we can read it back
  h5_obj <- H5ParcellatedMultiScan(temp_file)
  expect_s4_class(h5_obj, "H5ParcellatedMultiScan")
  
  # Check number of runs
  runs <- h5_obj@runs
  expect_length(runs, 1)
  expect_equal(runs[[1]]@scan_name, "multi_test")
  
  # Check it's a summary scan
  expect_s4_class(runs[[1]], "H5ParcellatedScanSummary")
  
  # Check time series dimensions
  expect_equal(runs[[1]]@n_time, n_time)
  
  # Clean up
  close(h5_obj)
})

test_that("as_h5 method for ClusteredNeuroVec defaults to single scan", {
  skip_if_not_installed("neuroim2")
  
  library(neuroim2)
  library(fmristore)
  
  # Create a simple test ClusteredNeuroVec
  mask_array <- array(FALSE, dim = c(5, 5, 5))
  mask_array[2:4, 2:4, 2:4] <- TRUE
  mask <- neuroim2::LogicalNeuroVol(mask_array, neuroim2::NeuroSpace(c(5, 5, 5)))
  
  n_voxels <- sum(mask_array)
  n_clusters <- 3
  cluster_ids <- sample(1:n_clusters, n_voxels, replace = TRUE)
  
  cvol <- neuroim2::ClusteredNeuroVol(mask = mask, clusters = cluster_ids)
  
  n_time <- 50
  ts_data <- matrix(rnorm(n_time * n_clusters), nrow = n_time, ncol = n_clusters)
  
  cnvec <- neuroim2::ClusteredNeuroVec(x = ts_data, cvol = cvol)
  
  # Test as_h5 with explicit file (defaults to single scan)
  temp_file <- tempfile(fileext = ".h5")
  on.exit(unlink(temp_file), add = TRUE)
  
  result <- as_h5(cnvec, file = temp_file)
  
  expect_s4_class(result, "H5ParcellatedScanSummary")  # Single scan by default
  expect_true(file.exists(temp_file))
  
  # Test as_h5 without file (should use temp file)
  result2 <- as_h5(cnvec, file = NULL)
  expect_s4_class(result2, "H5ParcellatedScanSummary")  # Single scan by default
  
  # Test as_h5 with as_multiscan=TRUE
  temp_file3 <- tempfile(fileext = ".h5")
  result3 <- as_h5(cnvec, file = temp_file3, as_multiscan = TRUE)
  expect_s4_class(result3, "H5ParcellatedMultiScan")  # Multi-scan when requested
  
  # Clean up
  close(result)
  close(result2)
  close(result3)
  unlink(h5file(result2))
  unlink(temp_file3)
})

test_that("ClusteredNeuroVec conversion handles edge cases", {
  skip_if_not_installed("neuroim2")
  
  library(neuroim2)
  library(fmristore)
  
  # Test with single cluster
  mask_array <- array(FALSE, dim = c(5, 5, 5))
  mask_array[2:4, 2:4, 2:4] <- TRUE
  mask <- neuroim2::LogicalNeuroVol(mask_array, neuroim2::NeuroSpace(c(5, 5, 5)))
  
  n_voxels <- sum(mask_array)
  cluster_ids <- rep(1, n_voxels)  # All voxels in one cluster
  
  cvol <- neuroim2::ClusteredNeuroVol(mask = mask, clusters = cluster_ids)
  
  n_time <- 20
  ts_data <- matrix(rnorm(n_time * 1), nrow = n_time, ncol = 1)
  
  cnvec <- neuroim2::ClusteredNeuroVec(x = ts_data, cvol = cvol)
  
  temp_file <- tempfile(fileext = ".h5")
  on.exit(unlink(temp_file), add = TRUE)
  
  result <- write_dataset(cnvec, file = temp_file)
  expect_s4_class(result, "H5ParcellatedScanSummary")  # Now defaults to single scan
  
  # Check single cluster was preserved
  expect_equal(result@n_time, n_time)
  
  close(result)
})

test_that("ClusteredNeuroVec conversion validates input", {
  skip_if_not_installed("neuroim2")
  
  library(neuroim2)
  library(fmristore)
  
  # Create valid ClusteredNeuroVec
  mask_array <- array(FALSE, dim = c(5, 5, 5))
  mask_array[2:4, 2:4, 2:4] <- TRUE
  mask <- neuroim2::LogicalNeuroVol(mask_array, neuroim2::NeuroSpace(c(5, 5, 5)))
  
  n_voxels <- sum(mask_array)
  n_clusters <- 3
  cluster_ids <- sample(1:n_clusters, n_voxels, replace = TRUE)
  
  cvol <- neuroim2::ClusteredNeuroVol(mask = mask, clusters = cluster_ids)
  
  n_time <- 50
  ts_data <- matrix(rnorm(n_time * n_clusters), nrow = n_time, ncol = n_clusters)
  
  cnvec <- neuroim2::ClusteredNeuroVec(x = ts_data, cvol = cvol)
  
  # Test missing file argument
  expect_error(write_dataset(cnvec), "'file' argument is required")
  
  # Test with custom scan name
  temp_file <- tempfile(fileext = ".h5")
  on.exit(unlink(temp_file), add = TRUE)
  
  result <- write_dataset(cnvec, file = temp_file, scan_name = "custom_scan")
  expect_s4_class(result, "H5ParcellatedScanSummary")  # Now defaults to single scan
  expect_equal(result@scan_name, "custom_scan")
  
  close(result)
})

test_that("ClusteredNeuroVec conversion preserves data integrity", {
  skip_if_not_installed("neuroim2")
  
  library(neuroim2)
  library(fmristore)
  
  # Create ClusteredNeuroVec with known data
  mask_array <- array(FALSE, dim = c(4, 4, 4))
  mask_array[2:3, 2:3, 2:3] <- TRUE
  mask <- neuroim2::LogicalNeuroVol(mask_array, neuroim2::NeuroSpace(c(4, 4, 4)))
  
  n_voxels <- sum(mask_array)
  n_clusters <- 2
  # Ensure we have both clusters represented
  cluster_ids <- c(rep(1, floor(n_voxels/2)), rep(2, ceiling(n_voxels/2)))
  
  cvol <- neuroim2::ClusteredNeuroVol(mask = mask, clusters = cluster_ids)
  
  n_time <- 10
  # Create known time series data
  ts_data <- matrix(seq_len(n_time * n_clusters), nrow = n_time, ncol = n_clusters)
  
  cnvec <- neuroim2::ClusteredNeuroVec(x = ts_data, cvol = cvol)
  
  temp_file <- tempfile(fileext = ".h5")
  on.exit(unlink(temp_file), add = TRUE)
  
  result <- write_dataset(cnvec, file = temp_file)
  
  # Check that result is single scan summary
  expect_s4_class(result, "H5ParcellatedScanSummary")
  expect_equal(result@n_time, n_time)
  
  # Check cluster IDs match
  expect_equal(length(result@cluster_ids), n_clusters)
  expect_true(all(result@cluster_ids %in% 1:n_clusters))
  
  close(result)
})