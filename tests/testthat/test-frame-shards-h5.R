.example_sharded_frame <- function() {
  observations <- fmridataset::axis_frame(
    data.frame(
      .obs_id = paste0("obs-", 1:6),
      subject_id = factor(rep(c("sub-01", "sub-02", "sub-03"), each = 2L)),
      condition = factor(rep(c("A", "B"), 3L), levels = c("A", "B"))
    ),
    blocks = list(
      motion = fmridataset::axis_block(
        matrix(seq_len(12), nrow = 6L),
        components = data.frame(.component_id = c("x", "y")),
        role = "confound"
      )
    )
  )
  space <- fmridataset::volume_space(
    dim = c(2L, 2L, 2L), support = c(1L, 3L, 4L, 7L), template = "shard-test"
  )
  beta <- matrix(seq_len(24), nrow = 6L)
  fmridataset::fmri_frame(
    assays = list(beta = beta, variance = beta / 10),
    observations = observations,
    features = fmridataset::feature_axis(
      data.frame(.feature_id = fmridataset::feature_ids(space)), space = space
    ),
    active_assay = "beta",
    metadata = list(study = "append-test")
  )
}

test_that("frame shards append without rewriting prior HDF5 data", {
  frame <- .example_sharded_frame()
  first <- frame[1:2, ]
  second <- frame[3:4, ]
  path <- tempfile("frame-shards-")
  on.exit(unlink(path, recursive = TRUE), add = TRUE)

  expect_identical(
    write_sharded_frame_h5(first, path, shard_id = "sub-01"),
    normalizePath(path)
  )
  before <- frame_shard_manifest(path)
  first_file <- file.path(path, before$.file[[1L]])
  first_hash <- unname(tools::md5sum(first_file))

  expect_identical(
    append_frame_shard_h5(second, path, shard_id = "sub-02"),
    normalizePath(path)
  )
  after <- frame_shard_manifest(path)
  expect_identical(after$.shard_id, c("sub-01", "sub-02"))
  expect_identical(after$.n_observation, c(2L, 2L))
  expect_identical(unname(tools::md5sum(first_file)), first_hash)

  reopened <- open_sharded_frame_h5(path)
  expect_s3_class(fmridataset::assay(reopened, "beta")$source, "row_sharded_source")
  expect_identical(
    fmridataset::shard_manifest(fmridataset::assay(reopened, "beta")$source)$.shard_id,
    c("sub-01", "sub-02")
  )
  expected <- fmridataset::bind_observations(first, second)
  expect_identical(fmridataset::observation_ids(reopened), fmridataset::observation_ids(expected))
  expect_identical(fmridataset::observations(reopened), fmridataset::observations(expected))
  expect_equal(
    fmridataset::collect_assay(reopened, "beta"),
    fmridataset::collect_assay(expected, "beta")
  )
  expect_equal(
    fmridataset::source_read(
      fmridataset::axis_block_data(fmridataset::obs_blocks(reopened)$motion)
    ),
    fmridataset::axis_block_data(fmridataset::obs_blocks(expected)$motion)
  )
})

test_that("failed shard append leaves the committed store byte-identical", {
  frame <- .example_sharded_frame()
  path <- tempfile("frame-shards-")
  on.exit(unlink(path, recursive = TRUE), add = TRUE)
  write_sharded_frame_h5(frame[1:2, ], path, shard_id = "sub-01")
  files_before <- list.files(path, recursive = TRUE, full.names = TRUE)
  hashes_before <- unname(tools::md5sum(files_before))

  withr::local_options(fmristore.shard_writer_fault_after_file = TRUE)
  expect_error(
    append_frame_shard_h5(frame[3:4, ], path, shard_id = "sub-02"),
    "Injected"
  )
  files_after <- list.files(path, recursive = TRUE, full.names = TRUE)
  expect_identical(files_after, files_before)
  expect_identical(unname(tools::md5sum(files_after)), hashes_before)
  expect_identical(frame_shard_manifest(path)$.shard_id, "sub-01")
  expect_equal(
    fmridataset::collect_assay(open_sharded_frame_h5(path)),
    fmridataset::collect_assay(frame[1:2, ])
  )
})

test_that("sharded frames reject incompatible or ambiguous appends before writing", {
  frame <- .example_sharded_frame()
  path <- tempfile("frame-shards-")
  on.exit(unlink(path, recursive = TRUE), add = TRUE)
  write_sharded_frame_h5(frame[1:2, ], path, shard_id = "sub-01")
  files_before <- list.files(path, recursive = TRUE)

  expect_error(
    append_frame_shard_h5(frame[3:4, 1:3], path, shard_id = "sub-02"),
    "incompatible"
  )
  expect_error(
    append_frame_shard_h5(frame[1:2, ], path, shard_id = "sub-02"),
    "collide"
  )
  expect_error(
    append_frame_shard_h5(frame[3:4, ], path, shard_id = "sub-01"),
    "already contains"
  )
  expect_error(
    append_frame_shard_h5(frame[3:4, ], path, shard_id = "bad/id"),
    "shard_id"
  )
  expect_identical(list.files(path, recursive = TRUE), files_before)
})
