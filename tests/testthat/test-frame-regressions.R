.regression_frame <- function(n = 4L, assays = NULL) {
  observations <- fmridataset::axis_frame(
    data.frame(
      .obs_id = paste0("obs-", seq_len(n)),
      condition = factor(rep(c("A", "B"), length.out = n), levels = c("A", "B"))
    ),
    blocks = list(
      motion = fmridataset::axis_block(
        matrix(as.double(seq_len(2L * n)), nrow = n),
        components = data.frame(.component_id = c("x", "y")),
        role = "confound"
      )
    )
  )
  space <- fmridataset::volume_space(
    dim = c(2L, 2L, 2L), support = c(1L, 3L, 4L, 7L), template = "regression"
  )
  if (is.null(assays)) {
    beta <- matrix(as.double(seq_len(4L * n)), nrow = n)
    assays <- list(beta = beta, variance = beta / 10)
  }
  fmridataset::fmri_frame(
    assays = assays,
    observations = observations,
    features = fmridataset::feature_axis(
      data.frame(.feature_id = fmridataset::feature_ids(space)), space = space
    ),
    active_assay = names(assays)[[1L]],
    metadata = list(study = "regression")
  )
}

.store_hashes <- function(path) {
  files <- list.files(path, recursive = TRUE, all.files = TRUE, full.names = TRUE)
  stats::setNames(unname(tools::md5sum(files)), files)
}

test_that("recovery on a healthy frame file finds and removes nothing", {
  path <- tempfile(fileext = ".fds.h5")
  on.exit(unlink(path), add = TRUE)
  write_frame_h5(.regression_frame(), path)
  before <- unname(tools::md5sum(path))

  for (action in c("report", "rollback")) {
    plan <- recover_frame_store(path, action)
    expect_s3_class(plan, "data.frame")
    expect_identical(nrow(plan), 0L)
    expect_named(plan, c("path", "kind", "removable", "removed"))
    expect_type(plan$removed, "logical")
  }
  expect_identical(unname(tools::md5sum(path)), before)
  expect_identical(frame_store_info(path)$status, "complete")
})

test_that("recovery on a healthy sharded store finds and removes nothing", {
  frame <- .regression_frame()
  path <- tempfile("healthy-shards-")
  on.exit(unlink(path, recursive = TRUE), add = TRUE)
  write_sharded_frame_h5(frame[1:2, ], path, shard_id = "a")
  append_frame_shard_h5(frame[3:4, ], path, shard_id = "b")
  before <- .store_hashes(path)

  for (action in c("report", "rollback")) {
    plan <- recover_frame_store(path, action)
    expect_identical(nrow(plan), 0L)
    expect_named(plan, c("path", "kind", "removable", "removed"))
  }
  expect_identical(.store_hashes(path), before)
  expect_identical(frame_store_info(path)$status, "complete")
})

test_that("rollback never deletes look-alike files outside a valid store", {
  path <- tempfile("not-a-store-")
  dir.create(file.path(path, "shards"), recursive = TRUE)
  on.exit(unlink(path, recursive = TRUE), add = TRUE)
  notes <- file.path(path, ".manifest-user-notes.txt")
  config <- file.path(path, "shards", ".shard-keep.cfg")
  writeLines("keep me", notes)
  writeLines("keep me too", config)

  expect_identical(frame_store_info(path)$status, "invalid")
  plan <- recover_frame_store(path, "rollback")
  expect_identical(nrow(plan), 0L)
  expect_true(file.exists(notes))
  expect_true(file.exists(config))
})

test_that("h5_read_plan sizes very large sources without integer overflow", {
  path <- tempfile(fileext = ".h5")
  on.exit(unlink(path), add = TRUE)
  h5 <- hdf5r::H5File$new(path, mode = "w")
  # Chunked and never written, so the file stays tiny.
  dataset <- h5$create_dataset(
    "beta",
    dtype = hdf5r::h5types$H5T_IEEE_F32LE,
    dims = c(200000L, 50000L),
    chunk_dims = c(128L, 8192L)
  )
  dataset$close()
  h5$close_all()
  source <- h5_array_source(path, "beta")
  expect_identical(source$shape, c(200000L, 50000L))

  full <- h5_read_plan(source)
  expect_identical(full$requested_values, 1e10)
  expect_identical(full$physical_values, 1e10)
  expect_identical(full$requested_bytes, 4e10)
  expect_identical(full$amplification, 1)

  block <- h5_read_plan(source, observations = 1:50000, features = 1:50000)
  expect_identical(block$requested_values, 2.5e9)
  # 391 observation chunks of 128; the feature axis is covered exactly.
  expect_identical(block$physical_values, 50048 * 50000)
  expect_lt(block$amplification, 1.2)
})

test_that("a crash while appending a shard leaves only artifacts recovery removes", {
  frame <- .regression_frame()
  path <- tempfile("crash-shards-")
  saved <- tempfile("crash-artifacts-")
  dir.create(saved)
  on.exit(unlink(c(path, saved), recursive = TRUE), add = TRUE)
  write_sharded_frame_h5(frame[1:2, ], path, shard_id = "a")
  committed <- .store_hashes(path)
  shard_dir <- file.path(path, "shards")

  # A killed process never runs its cleanup handlers. Capture what is on disk
  # at the moment of the injected failure, before the writer unwinds.
  artifacts <- character()
  capture <- function(condition) {
    present <- list.files(shard_dir, all.files = TRUE, full.names = TRUE, no.. = TRUE)
    artifacts <<- setdiff(present, names(committed))
    file.copy(artifacts, saved)
  }
  withr::with_options(
    list(fmristore.frame_writer_fault_after_array = 1L),
    expect_error(
      withCallingHandlers(
        append_frame_shard_h5(frame[3:4, ], path, shard_id = "b"),
        error = capture
      ),
      "Injected"
    )
  )
  expect_gte(length(artifacts), 1L)
  expect_identical(.store_hashes(path), committed)
  file.copy(file.path(saved, basename(artifacts)), artifacts)

  info <- frame_store_info(path)
  reported <- info$recovery[info$recovery$kind == "shard_partial", , drop = FALSE]
  expect_setequal(reported$path, normalizePath(artifacts, winslash = "/"))
  expect_true(all(reported$removable))

  rolled_back <- recover_frame_store(path, "rollback")
  expect_true(all(rolled_back$removed))
  expect_false(any(file.exists(artifacts)))
  expect_identical(.store_hashes(path), committed)

  append_frame_shard_h5(frame[3:4, ], path, shard_id = "b")
  expect_identical(nrow(recover_frame_store(path, "report")), 0L)
  expect_equal(
    fmridataset::collect_assay(open_sharded_frame_h5(path)),
    fmridataset::collect_assay(frame)
  )
})

test_that("logical arrays round trip as logical with NA preserved", {
  n <- 4L
  mask <- matrix(rep(c(TRUE, FALSE, NA, TRUE), length.out = 4L * n), nrow = n)
  frame <- .regression_frame(n, assays = list(mask = mask))
  path <- tempfile(fileext = ".fds.h5")
  on.exit(unlink(path), add = TRUE)
  write_frame_h5(frame, path, chunk_dims = c(3L, 3L))

  reopened <- open_frame_h5(path)
  source <- fmridataset::assay(reopened, "mask")$source
  expect_identical(fmridataset::source_dtype(source), "logical")
  values <- fmridataset::collect_assay(reopened, "mask")
  expect_type(values, "logical")
  expect_identical(unname(values), unname(mask))
  expect_identical(sum(is.na(values)), sum(is.na(mask)))
})

test_that("h5_array_source reports logical for hdf5r logical enums", {
  path <- tempfile(fileext = ".h5")
  on.exit(unlink(path), add = TRUE)
  mask <- matrix(c(TRUE, FALSE, NA, TRUE), 2L, 2L)
  h5 <- hdf5r::H5File$new(path, mode = "w")
  h5[["mask"]] <- mask
  dataset <- h5$create_dataset(
    "mask_na", dtype = hdf5r::h5types$H5T_LOGICAL_NA, dims = c(2L, 2L)
  )
  dataset[, ] <- mask
  dataset$close()
  h5$close_all()

  for (name in c("mask", "mask_na")) {
    source <- h5_array_source(path, name)
    expect_identical(fmridataset::source_dtype(source), "logical")
    expect_identical(fmridataset::source_read(source), mask)
    expect_identical(h5_read_plan(source)$requested_bytes, 4)
  }
})

test_that("narrow integer dtypes reject values they cannot store", {
  path <- tempfile(fileext = ".h5")
  on.exit(unlink(path), add = TRUE)
  h5 <- hdf5r::H5File$new(path, mode = "w")
  on.exit(if (h5$is_valid) h5$close_all(), add = TRUE, after = FALSE)
  write_case <- function(name, dtype, values) {
    fmristore:::.write_frame_array(
      h5, values, list(key = name, shape = dim(values), dtype = dtype),
      name, memory_budget = 1e6
    )
  }

  expect_error(write_case("u8_na", "uint8", matrix(c(1, NA, 3, 4), 2L)), "NA")
  expect_error(write_case("u8_high", "uint8", matrix(c(1, 2, 300, 4), 2L)), "uint8")
  expect_error(write_case("u8_low", "uint8", matrix(c(1, 2, -1, 4), 2L)), "uint8")
  expect_error(write_case("u8_frac", "uint8", matrix(c(1, 2.5, 3, 4), 2L)), "uint8")
  expect_error(write_case("i8_low", "int8", matrix(c(1L, -129L, 3L, 4L), 2L)), "int8")
  expect_error(write_case("u16_na", "uint16", matrix(c(1L, NA, 3L, 4L), 2L)), "NA")
  expect_error(write_case("i16_high", "int16", matrix(c(1L, NA, 40000L, -5L), 2L)), "int16")
  expect_error(write_case("u32_low", "uint32", matrix(c(1, -1, 3, 4), 2L)), "uint32")
  expect_error(write_case("u32_na", "uint32", matrix(c(1, NA, 3, 4), 2L)), "NA")

  valid <- list(
    uint8 = c(0, 255, 1, 2), int8 = c(-128, 127, 1, 2),
    uint16 = c(0, 65535, 1, 2), int16 = c(-32768, 32767, 1, 2),
    uint32 = c(0, 2^32 - 1, 1, 2), int32 = c(1L, NA, -5L, .Machine$integer.max)
  )
  for (dtype in names(valid)) {
    write_case(paste0("ok_", dtype), dtype, matrix(valid[[dtype]], 2L))
  }
  h5$close_all()
  for (dtype in names(valid)) {
    expect_equal(
      fmridataset::source_read(h5_array_source(path, paste0("ok_", dtype))),
      matrix(valid[[dtype]], 2L)
    )
  }
})

test_that("frame writes reject unrepresentable integer data before commit", {
  values <- matrix(c(1, 2, NA, 4), 4L, 4L)
  source <- fmridataset::memory_source(values, dtype = "uint8")
  frame <- .regression_frame(4L, assays = list(counts = source))
  parent <- tempfile("narrow-write-")
  dir.create(parent)
  on.exit(unlink(parent, recursive = TRUE), add = TRUE)
  path <- file.path(parent, "counts.fds.h5")

  expect_error(write_frame_h5(frame, path), "uint8")
  expect_false(file.exists(path))
  expect_length(list.files(parent, all.files = TRUE, no.. = TRUE), 0L)

  write_frame_h5(.regression_frame(), path)
  before <- unname(tools::md5sum(path))
  expect_error(write_frame_h5(frame, path, overwrite = TRUE), "uint8")
  expect_identical(unname(tools::md5sum(path)), before)
  expect_identical(list.files(parent, all.files = TRUE, no.. = TRUE), basename(path))
})
