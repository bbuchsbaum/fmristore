# Test neurovecseq_to_h5 function for NeuroVecSeq
library(testthat)
library(fmristore)
library(neuroim2)

test_that("neurovecseq_to_h5 works for NeuroVecSeq", {
  skip_if_not_installed("neuroim2")
  skip_if_not_installed("hdf5r")

  # Create test data - 3 NeuroVec objects with different time dimensions
  dims1 <- c(10, 10, 5, 20)  # 20 time points
  dims2 <- c(10, 10, 5, 30)  # 30 time points
  dims3 <- c(10, 10, 5, 25)  # 25 time points

  # Create NeuroVec objects using a helper function if available
  if (exists("create_minimal_DenseNeuroVec", where = asNamespace("fmristore"))) {
    vec1 <- fmristore:::create_minimal_DenseNeuroVec(dims1)
    vec2 <- fmristore:::create_minimal_DenseNeuroVec(dims2)
    vec3 <- fmristore:::create_minimal_DenseNeuroVec(dims3)
  } else {
    # Fallback: create manually
    data1 <- array(rnorm(prod(dims1)), dim = dims1)
    data2 <- array(rnorm(prod(dims2)), dim = dims2)
    data3 <- array(rnorm(prod(dims3)), dim = dims3)

    # Create 4D NeuroSpace objects
    space1 <- NeuroSpace(dims1)
    space2 <- NeuroSpace(dims2)
    space3 <- NeuroSpace(dims3)

    vec1 <- NeuroVec(data1, space1)
    vec2 <- NeuroVec(data2, space2)
    vec3 <- NeuroVec(data3, space3)
  }

  # Create NeuroVecSeq
  nvs <- NeuroVecSeq(vec1, vec2, vec3)

  # Check that we have a NeuroVecSeq
  expect_s4_class(nvs, "NeuroVecSeq")

  # Test basic conversion
  temp_file <- tempfile(fileext = ".h5")
  on.exit(unlink(temp_file), add = TRUE)

  # Debug: Check method dispatch
  # message("Class of nvs: ", class(nvs))
  # message("Methods for as_h5: ", paste(showMethods("as_h5", printTo = FALSE), collapse = "\n"))

  result <- neurovecseq_to_h5(nvs, file = temp_file,
    scan_names = c("run1", "run2", "run3"),
    compression = 4)

  expect_equal(result, temp_file)
  expect_true(file.exists(temp_file))

  # Verify file structure
  h5f <- hdf5r::H5File$new(temp_file, mode = "r")
  on.exit(try(h5f$close_all(), silent = TRUE), add = TRUE)

  # Check global attributes
  expect_equal(h5attr(h5f, "rtype"), "NeuroVecSeq")
  expect_equal(h5attr(h5f, "n_scans"), 3)

  # Check space information
  expect_true(h5f$exists("/space"))
  expect_equal(h5f[["space/dim"]][], dims1[1:3])

  # Check scans structure
  expect_true(h5f$exists("/scans"))
  expect_true(h5f$exists("/scans/run1"))
  expect_true(h5f$exists("/scans/run2"))
  expect_true(h5f$exists("/scans/run3"))

  # Check time dimensions
  expect_equal(h5attr(h5f[["scans/run1"]], "n_time"), 20)
  expect_equal(h5attr(h5f[["scans/run2"]], "n_time"), 30)
  expect_equal(h5attr(h5f[["scans/run3"]], "n_time"), 25)

  # Check data dimensions
  expect_equal(h5f[["scans/run1/data"]]$dims, dims1)
  expect_equal(h5f[["scans/run2/data"]]$dims, dims2)
  expect_equal(h5f[["scans/run3/data"]]$dims, dims3)
})

test_that("neurovecseq_to_h5 handles metadata", {
  skip_if_not_installed("neuroim2")
  skip_if_not_installed("hdf5r")

  # Create simple test data
  dims <- c(5, 5, 3, 10)
  if (exists("create_minimal_DenseNeuroVec", where = asNamespace("fmristore"))) {
    vec1 <- fmristore:::create_minimal_DenseNeuroVec(dims)
    vec1@.Data[] <- 1  # Fill with 1s
    vec2 <- fmristore:::create_minimal_DenseNeuroVec(dims)
    vec2@.Data[] <- 2  # Fill with 2s
  } else {
    vec1 <- NeuroVec(array(1, dim = dims), NeuroSpace(dims))
    vec2 <- NeuroVec(array(2, dim = dims), NeuroSpace(dims))
  }

  nvs <- NeuroVecSeq(vec1, vec2)

  # Create metadata
  scan_metadata <- list(
    scan_1 = list(TR = 2.0, task = "rest", subject = "sub01"),
    scan_2 = list(TR = 2.5, task = "motor", subject = "sub01")
  )

  temp_file <- tempfile(fileext = ".h5")
  on.exit(unlink(temp_file), add = TRUE)

  result <- neurovecseq_to_h5(nvs, file = temp_file, scan_metadata = scan_metadata)

  # Verify metadata was written
  h5f <- hdf5r::H5File$new(temp_file, mode = "r")
  on.exit(try(h5f$close_all(), silent = TRUE), add = TRUE)

  expect_true(h5f$exists("/scans/scan_1/metadata"))
  expect_equal(h5f[["scans/scan_1/metadata/TR"]][], 2.0)
  expect_equal(h5f[["scans/scan_1/metadata/task"]][], "rest")
  expect_equal(h5f[["scans/scan_2/metadata/TR"]][], 2.5)
  expect_equal(h5f[["scans/scan_2/metadata/task"]][], "motor")
})

test_that("neurovecseq_to_h5 validates inputs", {
  skip_if_not_installed("neuroim2")

  # Test with mismatched spatial dimensions
  if (exists("create_minimal_DenseNeuroVec", where = asNamespace("fmristore"))) {
    vec1 <- fmristore:::create_minimal_DenseNeuroVec(c(10, 10, 5, 20))
    vec2 <- fmristore:::create_minimal_DenseNeuroVec(c(8, 8, 5, 20))  # Different spatial dims

    # NeuroVecSeq constructor already validates spatial dimensions match
    expect_error(NeuroVecSeq(vec1, vec2), "All NeuroVec objects must have the same spatial dimensions")

    # Test with wrong number of scan names
    vec3 <- fmristore:::create_minimal_DenseNeuroVec(c(10, 10, 5, 20))
    nvs2 <- NeuroVecSeq(vec1, vec3)
  } else {
    vec1 <- NeuroVec(array(1, dim = c(10, 10, 5, 20)), NeuroSpace(c(10, 10, 5, 20)))
    vec2 <- NeuroVec(array(1, dim = c(8, 8, 5, 20)), NeuroSpace(c(8, 8, 5, 20)))

    # NeuroVecSeq constructor already validates spatial dimensions match
    expect_error(NeuroVecSeq(vec1, vec2), "All NeuroVec objects must have the same spatial dimensions")

    # Test with wrong number of scan names
    vec3 <- NeuroVec(array(1, dim = c(10, 10, 5, 20)), NeuroSpace(c(10, 10, 5, 20)))
    nvs2 <- NeuroVecSeq(vec1, vec3)
  }

  expect_error(neurovecseq_to_h5(nvs2, scan_names = c("run1")),
    "Length of scan_names must match number of NeuroVec objects")
})
