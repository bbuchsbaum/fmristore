# Load required packages
library(fmristore)
library(neuroim2)
library(hdf5r)

test_that("detect_h5_type correctly identifies H5 file types", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("neuroim2")
  
  # Test H5ParcellatedMultiScan detection
  temp_file <- tempfile(fileext = ".h5")
  on.exit(unlink(temp_file))
  
  # Create a simple parcellated experiment
  dims <- c(10, 10, 5)
  mask_data <- array(sample(c(TRUE, FALSE), prod(dims), replace = TRUE), dim = dims)
  mask <- neuroim2::LogicalNeuroVol(mask_data, neuroim2::NeuroSpace(dims))
  
  cluster_vec <- sample(0:3, sum(mask_data), replace = TRUE)
  clusters <- neuroim2::ClusteredNeuroVol(
    mask = mask,
    clusters = cluster_vec
  )
  
  # Create simple runs data
  runs_data <- list(
    list(
      scan_name = "run1",
      type = "summary",
      data = matrix(rnorm(20 * 4), nrow = 20, ncol = 4)
    )
  )
  
  # Write using new function
  write_parcellated_experiment_h5(
    filepath = temp_file,
    mask = mask,
    clusters = clusters,
    runs_data = runs_data,
    verbose = FALSE
  )
  
  # Test detection
  detected_type <- detect_h5_type(temp_file)
  expect_equal(detected_type, "H5ParcellatedMultiScan")
  
  # Test with H5File object
  h5 <- hdf5r::H5File$new(temp_file, "r")
  detected_type2 <- detect_h5_type(h5)
  h5$close_all()
  expect_equal(detected_type2, "H5ParcellatedMultiScan")
})

test_that("read_dataset auto-detects and reads H5ParcellatedMultiScan", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("neuroim2")
  
  temp_file <- tempfile(fileext = ".h5")
  on.exit(unlink(temp_file))
  
  # Create test data
  dims <- c(10, 10, 5)
  mask_data <- array(sample(c(TRUE, FALSE), prod(dims), replace = TRUE), dim = dims)
  mask <- neuroim2::LogicalNeuroVol(mask_data, neuroim2::NeuroSpace(dims))
  
  cluster_vec <- sample(0:3, sum(mask_data), replace = TRUE)
  clusters <- neuroim2::ClusteredNeuroVol(
    mask = mask,
    clusters = cluster_vec
  )
  
  runs_data <- list(
    list(
      scan_name = "run1",
      type = "summary",
      data = matrix(rnorm(20 * 4), nrow = 20, ncol = 4)
    )
  )
  
  write_parcellated_experiment_h5(
    filepath = temp_file,
    mask = mask,
    clusters = clusters,
    runs_data = runs_data,
    verbose = FALSE
  )
  
  # Test auto-detection
  obj <- read_dataset(temp_file)
  expect_s4_class(obj, "H5ParcellatedMultiScan")
  
  # Test explicit type
  obj2 <- read_dataset(temp_file, type = "H5ParcellatedMultiScan")
  expect_s4_class(obj2, "H5ParcellatedMultiScan")
  
  # Clean up
  close(obj)
  close(obj2)
})

test_that("write_dataset methods work for NeuroVec objects", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("neuroim2")
  
  temp_file <- tempfile(fileext = ".h5")
  on.exit(unlink(temp_file))
  
  # Create a NeuroVec
  dims <- c(10, 10, 5, 20)
  data <- array(rnorm(prod(dims)), dim = dims)
  space <- neuroim2::NeuroSpace(dims)
  nvec <- neuroim2::NeuroVec(data, space)
  
  # Create clusters
  mask_data <- array(TRUE, dim = dims[1:3])
  mask <- neuroim2::LogicalNeuroVol(mask_data, neuroim2::NeuroSpace(dims[1:3]))
  
  cluster_vec <- sample(1:3, prod(dims[1:3]), replace = TRUE)
  clusters <- neuroim2::ClusteredNeuroVol(
    mask = mask,
    clusters = cluster_vec
  )
  
  # Test write_dataset with summary
  write_dataset(nvec, file = temp_file, clusters = clusters, summary = TRUE)
  expect_true(file.exists(temp_file))
  
  # Read back and check
  obj <- read_dataset(temp_file)
  expect_s4_class(obj, "H5ParcellatedMultiScan")
  close(obj)
})

test_that("write_dataset methods work for list of NeuroVecs", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("neuroim2")
  
  temp_file <- tempfile(fileext = ".h5")
  on.exit(unlink(temp_file))
  
  # Create multiple NeuroVecs
  dims <- c(10, 10, 5, 20)
  nvecs <- lapply(1:3, function(i) {
    data <- array(rnorm(prod(dims)), dim = dims)
    space <- neuroim2::NeuroSpace(dims)
    neuroim2::NeuroVec(data, space)
  })
  
  # Create clusters
  mask_data <- array(TRUE, dim = dims[1:3])
  mask <- neuroim2::LogicalNeuroVol(mask_data, neuroim2::NeuroSpace(dims[1:3]))
  
  cluster_vec <- sample(1:3, prod(dims[1:3]), replace = TRUE)
  clusters <- neuroim2::ClusteredNeuroVol(
    mask = mask,
    clusters = cluster_vec
  )
  
  # Test write_dataset for list
  write_dataset(nvecs, file = temp_file, clusters = clusters, 
                mask = mask,  # Need to provide mask explicitly
                scan_names = c("scan1", "scan2", "scan3"),
                summary = FALSE)
  expect_true(file.exists(temp_file))
  
  # Read back and check
  obj <- read_dataset(temp_file)
  expect_s4_class(obj, "H5ParcellatedMultiScan")
  expect_equal(n_scans(obj), 3)
  close(obj)
})

test_that("write_dataset handles H5-backed objects correctly", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("neuroim2")
  
  temp_file1 <- tempfile(fileext = ".h5")
  temp_file2 <- tempfile(fileext = ".h5")
  on.exit({
    unlink(temp_file1)
    unlink(temp_file2)
  })
  
  # Create initial H5 file
  dims <- c(10, 10, 5)
  mask_data <- array(sample(c(TRUE, FALSE), prod(dims), replace = TRUE), dim = dims)
  mask <- neuroim2::LogicalNeuroVol(mask_data, neuroim2::NeuroSpace(dims))
  
  cluster_vec <- sample(0:3, sum(mask_data), replace = TRUE)
  clusters <- neuroim2::ClusteredNeuroVol(
    mask = mask,
    clusters = cluster_vec
  )
  
  runs_data <- list(
    list(
      scan_name = "run1",
      type = "summary",
      data = matrix(rnorm(20 * 4), nrow = 20, ncol = 4)
    )
  )
  
  write_parcellated_experiment_h5(
    filepath = temp_file1,
    mask = mask,
    clusters = clusters,
    runs_data = runs_data,
    verbose = FALSE
  )
  
  # Read as H5-backed object
  h5_obj <- H5ParcellatedMultiScan(temp_file1)
  
  # Test writing to same file (should return early with message)
  expect_message(
    result <- write_dataset(h5_obj, file = temp_file1),
    "already backed by target file"
  )
  expect_equal(result, temp_file1)
  
  # Test writing to different file (should copy)
  expect_message(
    write_dataset(h5_obj, file = temp_file2),
    "Copying H5ParcellatedMultiScan"
  )
  expect_true(file.exists(temp_file2))
  
  # Verify copy worked
  h5_obj2 <- H5ParcellatedMultiScan(temp_file2)
  expect_equal(n_scans(h5_obj), n_scans(h5_obj2))
  
  # Clean up
  close(h5_obj)
  close(h5_obj2)
})

test_that("detect_h5_type returns 'unknown' for unrecognized files", {
  skip_if_not_installed("hdf5r")
  
  temp_file <- tempfile(fileext = ".h5")
  on.exit(unlink(temp_file))
  
  # Create an H5 file with unrecognized structure
  h5 <- hdf5r::H5File$new(temp_file, "w")
  h5[["random_data"]] <- matrix(1:10, nrow = 2)
  h5$close_all()
  
  detected_type <- detect_h5_type(temp_file)
  expect_equal(detected_type, "unknown")
})

test_that("read_dataset fails gracefully for unknown types", {
  skip_if_not_installed("hdf5r")
  
  temp_file <- tempfile(fileext = ".h5")
  on.exit(unlink(temp_file))
  
  # Create an H5 file with unrecognized structure
  h5 <- hdf5r::H5File$new(temp_file, "w")
  h5[["random_data"]] <- matrix(1:10, nrow = 2)
  h5$close_all()
  
  expect_error(
    read_dataset(temp_file),
    "Cannot determine dataset type"
  )
})

test_that("deprecated write_clustered_experiment_h5 shows warning", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("neuroim2")
  
  temp_file <- tempfile(fileext = ".h5")
  on.exit(unlink(temp_file))
  
  # Create minimal test data
  dims <- c(10, 10, 5)
  mask_data <- array(sample(c(TRUE, FALSE), prod(dims), replace = TRUE), dim = dims)
  mask <- neuroim2::LogicalNeuroVol(mask_data, neuroim2::NeuroSpace(dims))
  
  cluster_vec <- sample(0:3, sum(mask_data), replace = TRUE)
  clusters <- neuroim2::ClusteredNeuroVol(
    mask = mask,
    clusters = cluster_vec
  )
  
  runs_data <- list(
    list(
      scan_name = "run1",
      type = "summary",
      data = matrix(rnorm(20 * 4), nrow = 20, ncol = 4)
    )
  )
  
  # Test that deprecated function shows warning
  expect_warning(
    write_clustered_experiment_h5(
      filepath = temp_file,
      mask = mask,
      clusters = clusters,
      runs_data = runs_data,
      verbose = FALSE
    ),
    "deprecated"
  )
  
  # But it should still work
  expect_true(file.exists(temp_file))
})

test_that("class attributes are written correctly", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("neuroim2")
  
  temp_file <- tempfile(fileext = ".h5")
  on.exit(unlink(temp_file))
  
  # Create test data
  dims <- c(10, 10, 5)
  mask_data <- array(sample(c(TRUE, FALSE), prod(dims), replace = TRUE), dim = dims)
  mask <- neuroim2::LogicalNeuroVol(mask_data, neuroim2::NeuroSpace(dims))
  
  cluster_vec <- sample(0:3, sum(mask_data), replace = TRUE)
  clusters <- neuroim2::ClusteredNeuroVol(
    mask = mask,
    clusters = cluster_vec
  )
  
  runs_data <- list(
    list(
      scan_name = "run1",
      type = "summary",
      data = matrix(rnorm(20 * 4), nrow = 20, ncol = 4)
    )
  )
  
  write_parcellated_experiment_h5(
    filepath = temp_file,
    mask = mask,
    clusters = clusters,
    runs_data = runs_data,
    verbose = FALSE
  )
  
  # Check that class attribute was written
  h5 <- hdf5r::H5File$new(temp_file, "r")
  expect_true("fmristore_class" %in% hdf5r::h5attr_names(h5))
  expect_equal(hdf5r::h5attr(h5, "fmristore_class"), "H5ParcellatedMultiScan")
  expect_true("fmristore_version" %in% hdf5r::h5attr_names(h5))
  h5$close_all()
})