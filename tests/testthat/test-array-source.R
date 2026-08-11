test_that("h5_array_source implements the reconstructible source contract", {
  path <- tempfile(fileext = ".h5")
  matrix <- matrix(seq_len(30), nrow = 5L, ncol = 6L)
  h5 <- hdf5r::H5File$new(path, mode = "w")
  h5$create_group("assays")
  h5$create_dataset("assays/beta", robj = matrix, chunk_dims = c(2L, 3L))
  h5$close_all()
  on.exit(unlink(path), add = TRUE)

  source <- h5_array_source(path, "/assays/beta")

  expect_s3_class(source, "array_source")
  expect_identical(fmridataset::source_shape(source), c(5L, 6L))
  expect_identical(fmridataset::source_chunks(source), c(2L, 3L))
  expect_true(fmridataset::source_dtype(source) %in% c("int32", "int64"))
  expect_setequal(
    fmridataset::source_capabilities(source),
    c("row_slice", "column_slice", "block_slice", "serializable")
  )
  expect_false(any(vapply(unclass(source), is.function, logical(1))))
  expect_false(any(vapply(unclass(source), is.environment, logical(1))))
  expect_false(any(vapply(
    unclass(source),
    function(x) typeof(x) == "externalptr",
    logical(1)
  )))
  expect_equal(
    fmridataset::source_read(
      source,
      observations = c(5L, 2L, 2L),
      features = c(6L, 1L, 3L)
    ),
    matrix[c(5L, 2L, 2L), c(6L, 1L, 3L), drop = FALSE]
  )

  rds <- tempfile(fileext = ".rds")
  saveRDS(source, rds)
  restored <- readRDS(rds)
  expect_equal(fmridataset::source_read(restored), matrix)
})

test_that("h5_array_source normalizes transposed physical axes", {
  path <- tempfile(fileext = ".h5")
  matrix <- matrix(seq_len(20), nrow = 4L, ncol = 5L)
  h5 <- hdf5r::H5File$new(path, mode = "w")
  h5$create_dataset("beta", robj = t(matrix))
  h5$close_all()
  on.exit(unlink(path), add = TRUE)

  source <- h5_array_source(
    path,
    "beta",
    physical_axes = c("feature", "observation")
  )
  expect_identical(fmridataset::source_shape(source), c(4L, 5L))
  expect_equal(
    fmridataset::source_read(source, 4:2, c(5L, 1L)),
    matrix[4:2, c(5L, 1L), drop = FALSE]
  )
})

test_that("explicit HDF5 source sessions close their handles", {
  path <- tempfile(fileext = ".h5")
  h5 <- hdf5r::H5File$new(path, mode = "w")
  h5[["beta"]] <- matrix(runif(12), 3L, 4L)
  h5$close_all()
  on.exit(unlink(path), add = TRUE)

  source <- h5_array_source(path, "beta")
  handle <- fmridataset::source_open(source)
  expect_true(handle$file$is_valid)
  expect_equal(
    fmridataset::source_read(handle, 1L, 2L),
    fmridataset::source_read(source, 1L, 2L)
  )
  fmridataset::source_close(handle)
  expect_false(handle$file$is_valid)
  expect_error(fmridataset::source_read(handle), "closed")
})

test_that("invalid HDF5 source descriptors fail before use", {
  expect_error(h5_array_source("missing.h5", "beta"), "does not exist")

  path <- tempfile(fileext = ".h5")
  h5 <- hdf5r::H5File$new(path, mode = "w")
  h5[["cube"]] <- array(1, c(2L, 2L, 2L))
  h5$close_all()
  on.exit(unlink(path), add = TRUE)

  expect_error(h5_array_source(path, "missing"), "does not exist")
  expect_error(h5_array_source(path, "cube"), "two dimensional")
  expect_error(
    h5_array_source(path, "cube", physical_axes = c("feature", "feature")),
    "exactly once"
  )
})
