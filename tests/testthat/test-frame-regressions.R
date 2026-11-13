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
