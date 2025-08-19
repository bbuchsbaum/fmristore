# Comprehensive tests for NeuroVecSeq to HDF5 functionality
library(testthat)
library(fmristore)
library(neuroim2)
library(hdf5r)

test_that("neurovecseq_to_h5 creates correct HDF5 structure", {
  skip_if_not_installed("neuroim2")
  skip_if_not_installed("hdf5r")

  # Create test data with known values
  dims1 <- c(10, 10, 5, 20)
  dims2 <- c(10, 10, 5, 30)
  dims3 <- c(10, 10, 5, 25)

  # Create data with specific patterns to verify later
  data1 <- array(seq_len(prod(dims1)), dim = dims1)
  data2 <- array(seq_len(prod(dims2)) * 2, dim = dims2)
  data3 <- array(seq_len(prod(dims3)) * 3, dim = dims3)

  # Create NeuroVec objects using helper or with proper NeuroSpace
  if (exists("create_minimal_DenseNeuroVec", where = asNamespace("fmristore"))) {
    vec1 <- fmristore:::create_minimal_DenseNeuroVec(dims1)
    vec1@.Data <- data1
    vec2 <- fmristore:::create_minimal_DenseNeuroVec(dims2)
    vec2@.Data <- data2
    vec3 <- fmristore:::create_minimal_DenseNeuroVec(dims3)
    vec3@.Data <- data3
  } else {
    vec1 <- NeuroVec(data1, NeuroSpace(dims1))
    vec2 <- NeuroVec(data2, NeuroSpace(dims2))
    vec3 <- NeuroVec(data3, NeuroSpace(dims3))
  }

  nvs <- NeuroVecSeq(vec1, vec2, vec3)

  # Convert to HDF5
  temp_file <- tempfile(fileext = ".h5")
  on.exit(unlink(temp_file), add = TRUE)

  result <- neurovecseq_to_h5(
    nvs,
    file = temp_file,
    scan_names = c("scan_A", "scan_B", "scan_C"),
    data_type = "FLOAT",
    chunk_dim = c(5, 5, 5, 10),
    compression = 4
  )

  expect_equal(result, temp_file)
  expect_true(file.exists(temp_file))

  # Verify file structure in detail
  h5f <- H5File$new(temp_file, mode = "r")
  on.exit(try(h5f$close_all(), silent = TRUE), add = TRUE)

  # Check root attributes
  expect_equal(h5attr(h5f, "rtype"), "NeuroVecSeq")
  expect_equal(h5attr(h5f, "n_scans"), 3)

  # Check space group
  expect_true(h5f$exists("/space"))
  expect_true(h5f$exists("/space/dim"))
  expect_true(h5f$exists("/space/origin"))
  expect_true(h5f$exists("/space/trans"))

  # Verify spatial dimensions
  space_dims <- h5f[["space/dim"]][]
  expect_equal(space_dims, c(10, 10, 5))

  # Check scans structure
  expect_true(h5f$exists("/scans"))
  expect_true(h5f$exists("/scans/scan_A"))
  expect_true(h5f$exists("/scans/scan_B"))
  expect_true(h5f$exists("/scans/scan_C"))

  # Verify scan attributes
  expect_equal(h5attr(h5f[["scans/scan_A"]], "n_time"), 20)
  expect_equal(h5attr(h5f[["scans/scan_B"]], "n_time"), 30)
  expect_equal(h5attr(h5f[["scans/scan_C"]], "n_time"), 25)

  # Check data arrays exist and have correct dimensions
  expect_true(h5f$exists("/scans/scan_A/data"))
  expect_equal(h5f[["scans/scan_A/data"]]$dims, dims1)
  expect_equal(h5f[["scans/scan_B/data"]]$dims, dims2)
  expect_equal(h5f[["scans/scan_C/data"]]$dims, dims3)

  # Verify actual data values for a sample
  scan_A_sample <- h5f[["scans/scan_A/data"]][1:2, 1:2, 1, 1]
  expected_sample <- data1[1:2, 1:2, 1, 1]
  expect_equal(scan_A_sample, expected_sample, tolerance = 1e-6)
})

test_that("neurovecseq_to_h5 handles metadata correctly", {
  skip_if_not_installed("neuroim2")
  skip_if_not_installed("hdf5r")

  # Create simple test data
  dims <- c(5, 5, 3, 10)
  vec1 <- NeuroVec(array(1, dim = dims), NeuroSpace(dims))
  vec2 <- NeuroVec(array(2, dim = dims), NeuroSpace(dims))
  vec3 <- NeuroVec(array(3, dim = dims), NeuroSpace(dims))

  nvs <- NeuroVecSeq(vec1, vec2, vec3)

  # Create comprehensive metadata
  scan_metadata <- list(
    run1 = list(
      TR = 2.0,
      TE = 30.0,
      flip_angle = 90,
      task = "rest",
      subject_id = "sub-01",
      session = "ses-01",
      acquisition_date = "2024-01-15",
      scanner = "Siemens Prisma",
      field_strength = 3.0,
      n_volumes = 10L,
      slice_order = c(1L, 3L, 5L, 2L, 4L)
    ),
    run2 = list(
      TR = 2.5,
      TE = 35.0,
      flip_angle = 85,
      task = "motor",
      subject_id = "sub-01",
      session = "ses-01",
      acquisition_date = "2024-01-15",
      n_trials = 40L,
      stimulus_type = "visual"
    ),
    run3 = list(
      TR = 2.0,
      task = "visual",
      notes = "Subject reported slight movement"
    )
  )

  temp_file <- tempfile(fileext = ".h5")
  on.exit(unlink(temp_file), add = TRUE)

  result <- neurovecseq_to_h5(
    nvs,
    file = temp_file,
    scan_names = c("run1", "run2", "run3"),
    scan_metadata = scan_metadata
  )

  # Verify metadata was written correctly
  h5f <- H5File$new(temp_file, mode = "r")
  on.exit(try(h5f$close_all(), silent = TRUE), add = TRUE)

  # Check run1 metadata (comprehensive)
  expect_true(h5f$exists("/scans/run1/metadata"))
  expect_equal(h5f[["scans/run1/metadata/TR"]][], 2.0)
  expect_equal(h5f[["scans/run1/metadata/TE"]][], 30.0)
  expect_equal(h5f[["scans/run1/metadata/flip_angle"]][], 90)
  expect_equal(h5f[["scans/run1/metadata/task"]][], "rest")
  expect_equal(h5f[["scans/run1/metadata/subject_id"]][], "sub-01")
  expect_equal(h5f[["scans/run1/metadata/field_strength"]][], 3.0)
  expect_equal(h5f[["scans/run1/metadata/n_volumes"]][], 10L)
  expect_equal(h5f[["scans/run1/metadata/slice_order"]][], c(1L, 3L, 5L, 2L, 4L))

  # Check run2 metadata (partial)
  expect_equal(h5f[["scans/run2/metadata/TR"]][], 2.5)
  expect_equal(h5f[["scans/run2/metadata/task"]][], "motor")
  expect_equal(h5f[["scans/run2/metadata/n_trials"]][], 40L)

  # Check run3 metadata (minimal)
  expect_equal(h5f[["scans/run3/metadata/task"]][], "visual")
  expect_equal(h5f[["scans/run3/metadata/notes"]][], "Subject reported slight movement")
})

test_that("neurovecseq_to_h5 chunking and compression work correctly", {
  skip_if_not_installed("neuroim2")
  skip_if_not_installed("hdf5r")

  dims <- c(20, 20, 10, 50)
  vec <- NeuroVec(array(rnorm(prod(dims)), dim = dims), NeuroSpace(dims))
  nvs <- NeuroVecSeq(vec)

  # Test 1: Custom chunk dimensions
  temp_file1 <- tempfile(fileext = ".h5")
  on.exit(unlink(temp_file1), add = TRUE)

  neurovecseq_to_h5(
    nvs,
    file = temp_file1,
    chunk_dim = c(10, 10, 5, 25),
    compression = 6
  )

  h5f1 <- H5File$new(temp_file1, mode = "r")
  on.exit(try(h5f1$close_all(), silent = TRUE), add = TRUE)

  # Check chunk dimensions
  dset1 <- h5f1[["scans/scan_1/data"]]
  chunk_info <- dset1$chunk_dims
  expect_equal(chunk_info, c(10, 10, 5, 25))

  # Test 2: Default (time-optimized) chunking
  temp_file2 <- tempfile(fileext = ".h5")
  on.exit(unlink(temp_file2), add = TRUE)

  neurovecseq_to_h5(
    nvs,
    file = temp_file2,
    chunk_dim = NULL, # Use default
    compression = 4
  )

  h5f2 <- H5File$new(temp_file2, mode = "r")
  on.exit(try(h5f2$close_all(), silent = TRUE), add = TRUE)

  dset2 <- h5f2[["scans/scan_1/data"]]
  chunk_info2 <- dset2$chunk_dims
  expect_equal(chunk_info2, c(10, 10, 10, 50)) # Default time-optimized

  # Test 3: File size comparison (compressed vs uncompressed)
  temp_file3 <- tempfile(fileext = ".h5")
  on.exit(unlink(temp_file3), add = TRUE)

  neurovecseq_to_h5(
    nvs,
    file = temp_file3,
    compression = 0 # No compression
  )

  size_compressed <- file.info(temp_file1)$size
  size_uncompressed <- file.info(temp_file3)$size

  # Compressed file should be smaller
  expect_true(size_compressed < size_uncompressed)

  h5f1$close_all()
  h5f2$close_all()
})

test_that("neurovecseq_to_h5 data integrity is maintained", {
  skip_if_not_installed("neuroim2")
  skip_if_not_installed("hdf5r")

  # Create data with specific patterns
  dims <- c(8, 8, 4, 15)

  # Different patterns for each scan
  data1 <- array(0, dim = dims)
  data1[4, 4, 2, ] <- 1:15 # Time series at specific voxel

  data2 <- array(0, dim = dims)
  data2[, , , 8] <- 1 # All voxels at time point 8

  data3 <- array(0, dim = dims)
  data3[1:4, 1:4, 1, 1] <- matrix(1:16, 4, 4) # Spatial pattern at t=1

  vec1 <- NeuroVec(data1, NeuroSpace(dims))
  vec2 <- NeuroVec(data2, NeuroSpace(dims))
  vec3 <- NeuroVec(data3, NeuroSpace(dims))

  nvs <- NeuroVecSeq(vec1, vec2, vec3)

  temp_file <- tempfile(fileext = ".h5")
  on.exit(unlink(temp_file), add = TRUE)

  neurovecseq_to_h5(
    nvs,
    file = temp_file,
    scan_names = c("pattern1", "pattern2", "pattern3"),
    compression = 4
  )

  # Read back and verify
  h5f <- H5File$new(temp_file, mode = "r")
  on.exit(try(h5f$close_all(), silent = TRUE), add = TRUE)

  # Test 1: Time series extraction
  ts_extracted <- h5f[["scans/pattern1/data"]][4, 4, 2, ]
  expect_equal(ts_extracted, 1:15)

  # Test 2: Full volume extraction
  vol_extracted <- h5f[["scans/pattern2/data"]][, , , 8]
  expect_true(all(vol_extracted == 1))

  # Test 3: Spatial pattern
  spatial_extracted <- h5f[["scans/pattern3/data"]][1:4, 1:4, 1, 1]
  expect_equal(spatial_extracted, matrix(1:16, 4, 4))

  # Test 4: Verify zeros elsewhere
  other_voxel <- h5f[["scans/pattern1/data"]][1, 1, 1, ]
  expect_true(all(other_voxel == 0))
})

test_that("neurovecseq_to_h5 handles edge cases", {
  skip_if_not_installed("neuroim2")
  skip_if_not_installed("hdf5r")

  # Test 1: Single scan NeuroVecSeq
  dims <- c(5, 5, 3, 10)
  vec <- NeuroVec(array(1, dim = dims), NeuroSpace(dims))
  nvs_single <- NeuroVecSeq(vec)

  temp_file1 <- tempfile(fileext = ".h5")
  on.exit(unlink(temp_file1), add = TRUE)

  result1 <- neurovecseq_to_h5(nvs_single, file = temp_file1)
  expect_true(file.exists(temp_file1))

  h5f1 <- H5File$new(temp_file1, mode = "r")
  on.exit(try(h5f1$close_all(), silent = TRUE), add = TRUE)
  expect_equal(h5attr(h5f1, "n_scans"), 1)
  expect_true(h5f1$exists("/scans/scan_1"))
  h5f1$close_all()

  # Test 2: Large number of scans
  vec_list <- lapply(1:10, function(i) {
    NeuroVec(array(i, dim = c(3, 3, 2, 5)), NeuroSpace(c(3, 3, 2, 5)))
  })
  nvs_many <- do.call(NeuroVecSeq, vec_list)

  temp_file2 <- tempfile(fileext = ".h5")
  on.exit(unlink(temp_file2), add = TRUE)

  result2 <- neurovecseq_to_h5(
    nvs_many,
    file = temp_file2,
    scan_names = paste0("scan_", letters[1:10])
  )

  h5f2 <- H5File$new(temp_file2, mode = "r")
  on.exit(try(h5f2$close_all(), silent = TRUE), add = TRUE)
  expect_equal(h5attr(h5f2, "n_scans"), 10)
  expect_true(h5f2$exists("/scans/scan_j"))
  h5f2$close_all()

  # Test 3: Different data types
  temp_file3 <- tempfile(fileext = ".h5")
  on.exit(unlink(temp_file3), add = TRUE)

  result3 <- neurovecseq_to_h5(
    nvs_single,
    file = temp_file3,
    data_type = "DOUBLE"
  )

  h5f3 <- H5File$new(temp_file3, mode = "r")
  on.exit(try(h5f3$close_all(), silent = TRUE), add = TRUE)

  # Check data type
  dset <- h5f3[["scans/scan_1/data"]]
  dtype_obj <- dset$get_type()
  dtype_size <- dtype_obj$get_size()
  expect_equal(dtype_size, 8) # DOUBLE is 8 bytes

  h5f3$close_all()
})

test_that("neurovecseq_to_h5 error handling works correctly", {
  skip_if_not_installed("neuroim2")
  skip_if_not_installed("hdf5r")

  # Test 1: Empty NeuroVecSeq - can't create empty NeuroVecSeq, so test with NULL input
  expect_error(
    neurovecseq_to_h5(NULL),
    "Input must be a NeuroVecSeq object"
  )

  # Test 2: Invalid input type
  expect_error(
    neurovecseq_to_h5(list(1, 2, 3)),
    "Input must be a NeuroVecSeq object"
  )

  # Test 3: Mismatched scan names length
  vec1 <- NeuroVec(array(1, dim = c(3, 3, 2, 5)), NeuroSpace(c(3, 3, 2, 5)))
  vec2 <- NeuroVec(array(2, dim = c(3, 3, 2, 5)), NeuroSpace(c(3, 3, 2, 5)))
  nvs <- NeuroVecSeq(vec1, vec2)

  expect_error(
    neurovecseq_to_h5(nvs, scan_names = c("only_one")),
    "Length of scan_names must match number of NeuroVec objects"
  )

  # Test 4: Invalid data type
  expect_error(
    neurovecseq_to_h5(nvs, data_type = "INVALID"),
    "Unsupported data_type"
  )

  # Test 5: Invalid compression level
  expect_error(
    neurovecseq_to_h5(nvs, compression = 10),
    "compression.*must be.*0-9"
  )

  # Test 6: Invalid chunk dimensions
  expect_error(
    neurovecseq_to_h5(nvs, chunk_dim = c(100, 100, 100, 100)), # Larger than data
    "chunk_dims.*must be.*<= dims"
  )
})

test_that("neurovecseq_to_h5 performance is acceptable", {
  skip_if_not_installed("neuroim2")
  skip_if_not_installed("hdf5r")
  skip_on_cran() # Skip performance tests on CRAN

  # Create reasonably sized data
  dims <- c(64, 64, 30, 100) # ~12MB per scan
  vec1 <- NeuroVec(array(rnorm(prod(dims)), dim = dims), NeuroSpace(dims))
  vec2 <- NeuroVec(array(rnorm(prod(dims)), dim = dims), NeuroSpace(dims))
  vec3 <- NeuroVec(array(rnorm(prod(dims)), dim = dims), NeuroSpace(dims))

  nvs <- NeuroVecSeq(vec1, vec2, vec3)

  temp_file <- tempfile(fileext = ".h5")
  on.exit(unlink(temp_file), add = TRUE)

  # Time the operation
  write_time <- system.time({
    neurovecseq_to_h5(nvs, file = temp_file, compression = 4)
  })["elapsed"]

  # Should complete in reasonable time (adjust based on system)
  expect_true(write_time < 30,
    info = sprintf("Write took %.1f seconds", write_time)
  )

  # Test read performance
  h5f <- H5File$new(temp_file, mode = "r")
  on.exit(try(h5f$close_all(), silent = TRUE), add = TRUE)

  # Time series access should be fast
  ts_time <- system.time({
    for (i in 1:10) {
      ts <- h5f[["scans/scan_1/data"]][32, 32, 15, ]
    }
  })["elapsed"]

  expect_true(ts_time < 1.0,
    info = sprintf("10 time series extractions took %.3f seconds", ts_time)
  )

  # Volume access should also be reasonable
  vol_time <- system.time({
    vol <- h5f[["scans/scan_2/data"]][, , , 50]
  })["elapsed"]

  expect_true(vol_time < 2.0,
    info = sprintf("Volume extraction took %.3f seconds", vol_time)
  )
})
