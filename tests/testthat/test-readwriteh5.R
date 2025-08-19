library(neuroim2)
library(fmristore)

context("HDF5 read/write functions")

test_that("round trip write/read preserves data integrity", {
  # Create test data with explicit dimensions
  dims <- c(10, 10, 10)
  nvols <- 10
  space <- NeuroSpace(c(dims, nvols), spacing = c(2, 2, 2), origin = c(1, 1, 1))

  # Create volumes with known pattern
  vols <- array(0, dim = c(dims, nvols))
  for (i in 1:nvols) {
    vols[, , , i] <- array(i + seq_len(prod(dims)) / 1000, dim = dims)
  }

  # Create NeuroVec and labels
  vec <- NeuroVec(vols, space)
  labels <- paste0("vol_", 1:nvols)

  # Create mask where values > 0.5 in first volume
  mask_array <- array(as.logical(vols[, , , 1] > 0.5), dim = dims) # Ensure 3D logical array
  mask_space <- NeuroSpace(dims, spacing = spacing(space)[1:3], origin = origin(space))
  mask <- LogicalNeuroVol(mask_array, mask_space)
  expect_true(inherits(mask, "LogicalNeuroVol"))
  expect_equal(dim(mask), dims) # Verify mask dimensions

  # Write to temporary file
  tmp <- tempfile(fileext = ".h5")
  on.exit(unlink(tmp))

  write_labeled_vec(vec, mask = mask, labels = labels, file = tmp)

  # Read back
  result <- read_labeled_vec(tmp)

  # Test space properties

  expect_equal(spacing(space(vec)), spacing(space(result[[1]])))
  expect_equal(trans(vec), trans(result[[1]]))

  # Test labels
  expect_equal(labels, names(result))

  # Test mask preservation
  expect_true(!is.null(attr(result, "mask")))
  expect_equal(dim(attr(result, "mask")), dims) # Verify read mask dimensions
  expect_equal(as.array(mask@.Data), as.array(attr(result, "mask")@.Data))

  # Test volume data
  for (i in 1:nvols) {
    expect_equal(as.array(vec[[i]]@.Data), as.array(result[[i]]@.Data))
  }
})
