library(fmristore)
library(neuroim2)

test_that("write_cluster_result works with list input", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("neuroim2")
  
  # Create test data
  brain_dims <- c(10, 10, 5)
  n_time <- 20
  
  # Create mask
  mask_data <- array(FALSE, brain_dims)
  mask_data[3:7, 3:7, 2:4] <- TRUE
  mask <- neuroim2::LogicalNeuroVol(mask_data, neuroim2::NeuroSpace(brain_dims))
  
  # Create cluster assignments
  n_voxels <- sum(mask_data)
  cluster_assignments <- rep(1:3, length.out = n_voxels)
  
  # Create neurovec
  fmri_data <- array(rnorm(prod(brain_dims) * n_time), dim = c(brain_dims, n_time))
  nvec <- neuroim2::NeuroVec(fmri_data, neuroim2::NeuroSpace(c(brain_dims, n_time)))
  
  # Create cluster_result as a list (simulating supervoxels output)
  cluster_result <- list(
    clusters = cluster_assignments,
    mask = mask
  )
  
  # Test full type
  temp_file <- tempfile(fileext = ".h5")
  expect_silent(
    write_cluster_result(
      cluster_result = cluster_result,
      neurovec = nvec,
      filename = temp_file,
      type = "full",
      verbose = FALSE
    )
  )
  expect_true(file.exists(temp_file))
  
  # Clean up
  unlink(temp_file)
  
  # Test summary type
  temp_file <- tempfile(fileext = ".h5")
  expect_silent(
    write_cluster_result(
      cluster_result = cluster_result,
      neurovec = nvec,
      filename = temp_file,
      type = "summary",
      verbose = FALSE
    )
  )
  expect_true(file.exists(temp_file))
  
  # Clean up
  unlink(temp_file)
})

test_that("write_cluster_result works with ClusteredNeuroVol input", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("neuroim2")
  
  # Create test data
  brain_dims <- c(10, 10, 5)
  n_time <- 20
  
  # Create mask
  mask_data <- array(FALSE, brain_dims)
  mask_data[3:7, 3:7, 2:4] <- TRUE
  mask <- neuroim2::LogicalNeuroVol(mask_data, neuroim2::NeuroSpace(brain_dims))
  
  # Create ClusteredNeuroVol directly
  n_voxels <- sum(mask_data)
  cluster_assignments <- rep(1:3, length.out = n_voxels)
  clusters <- neuroim2::ClusteredNeuroVol(mask, cluster_assignments)
  
  # Create neurovec
  fmri_data <- array(rnorm(prod(brain_dims) * n_time), dim = c(brain_dims, n_time))
  nvec <- neuroim2::NeuroVec(fmri_data, neuroim2::NeuroSpace(c(brain_dims, n_time)))
  
  # Test with ClusteredNeuroVol
  temp_file <- tempfile(fileext = ".h5")
  expect_silent(
    write_cluster_result(
      cluster_result = clusters,
      neurovec = nvec,
      filename = temp_file,
      type = "full",
      verbose = FALSE
    )
  )
  expect_true(file.exists(temp_file))
  
  # Clean up
  unlink(temp_file)
})

test_that("write_cluster_result handles compression correctly", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("neuroim2")
  
  # Create minimal test data
  brain_dims <- c(5, 5, 2)
  n_time <- 10
  
  mask_data <- array(TRUE, brain_dims)
  mask <- neuroim2::LogicalNeuroVol(mask_data, neuroim2::NeuroSpace(brain_dims))
  
  n_voxels <- sum(mask_data)
  clusters <- neuroim2::ClusteredNeuroVol(mask, rep(1:2, length.out = n_voxels))
  
  fmri_data <- array(rnorm(prod(brain_dims) * n_time), dim = c(brain_dims, n_time))
  nvec <- neuroim2::NeuroVec(fmri_data, neuroim2::NeuroSpace(c(brain_dims, n_time)))
  
  # Test logical TRUE
  temp_file <- tempfile(fileext = ".h5")
  expect_silent(
    write_cluster_result(
      cluster_result = clusters,
      neurovec = nvec,
      filename = temp_file,
      compress = TRUE,
      verbose = FALSE
    )
  )
  unlink(temp_file)
  
  # Test logical FALSE
  temp_file <- tempfile(fileext = ".h5")
  expect_silent(
    write_cluster_result(
      cluster_result = clusters,
      neurovec = nvec,
      filename = temp_file,
      compress = FALSE,
      verbose = FALSE
    )
  )
  unlink(temp_file)
  
  # Test numeric compression level
  temp_file <- tempfile(fileext = ".h5")
  expect_silent(
    write_cluster_result(
      cluster_result = clusters,
      neurovec = nvec,
      filename = temp_file,
      compress = 7,
      verbose = FALSE
    )
  )
  unlink(temp_file)
  
  # Test invalid compression level
  expect_error(
    write_cluster_result(
      cluster_result = clusters,
      neurovec = nvec,
      filename = tempfile(fileext = ".h5"),
      compress = 10,
      verbose = FALSE
    ),
    "Compression level must be between 0 and 9"
  )
})

test_that("write_cluster_result prevents overwriting by default", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("neuroim2")
  
  # Create minimal test data
  brain_dims <- c(5, 5, 2)
  mask_data <- array(TRUE, brain_dims)
  mask <- neuroim2::LogicalNeuroVol(mask_data, neuroim2::NeuroSpace(brain_dims))
  clusters <- neuroim2::ClusteredNeuroVol(mask, rep(1, sum(mask_data)))
  
  fmri_data <- array(1, dim = c(brain_dims, 5))
  nvec <- neuroim2::NeuroVec(fmri_data, neuroim2::NeuroSpace(c(brain_dims, 5)))
  
  # Create a file
  temp_file <- tempfile(fileext = ".h5")
  write_cluster_result(
    cluster_result = clusters,
    neurovec = nvec,
    filename = temp_file,
    verbose = FALSE
  )
  
  # Try to overwrite without flag
  expect_error(
    write_cluster_result(
      cluster_result = clusters,
      neurovec = nvec,
      filename = temp_file,
      verbose = FALSE
    ),
    "File already exists"
  )
  
  # Should work with overwrite=TRUE
  expect_silent(
    write_cluster_result(
      cluster_result = clusters,
      neurovec = nvec,
      filename = temp_file,
      overwrite = TRUE,
      verbose = FALSE
    )
  )
  
  # Clean up
  unlink(temp_file)
})

test_that("write_cluster_result validates input correctly", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("neuroim2")
  
  # Create valid components
  brain_dims <- c(5, 5, 2)
  mask_data <- array(TRUE, brain_dims)
  mask <- neuroim2::LogicalNeuroVol(mask_data, neuroim2::NeuroSpace(brain_dims))
  
  fmri_data <- array(1, dim = c(brain_dims, 5))
  nvec <- neuroim2::NeuroVec(fmri_data, neuroim2::NeuroSpace(c(brain_dims, 5)))
  
  # Test missing clusters element
  expect_error(
    write_cluster_result(
      cluster_result = list(mask = mask),
      neurovec = nvec,
      filename = tempfile(fileext = ".h5"),
      verbose = FALSE
    ),
    "must contain 'clusters' and 'mask' elements"
  )
  
  # Test mismatched cluster vector length
  expect_error(
    write_cluster_result(
      cluster_result = list(
        clusters = 1:10,  # Wrong length
        mask = mask
      ),
      neurovec = nvec,
      filename = tempfile(fileext = ".h5"),
      verbose = FALSE
    ),
    "doesn't match number of TRUE voxels"
  )
  
  # Test invalid neurovec
  expect_error(
    write_cluster_result(
      cluster_result = list(clusters = rep(1, sum(mask_data)), mask = mask),
      neurovec = "not_a_neurovec",
      filename = tempfile(fileext = ".h5"),
      verbose = FALSE
    ),
    "neurovec must be a NeuroVec"
  )
})