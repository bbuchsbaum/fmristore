.recovery_frame <- function() {
  observations <- fmridataset::axis_frame(
    data.frame(
      .obs_id = paste0("obs-", 1:4),
      group = factor(c("A", "A", "B", "B"), levels = c("A", "B"))
    ),
    blocks = list(
      motion = fmridataset::axis_block(
        matrix(seq_len(8), nrow = 4L),
        components = data.frame(.component_id = c("x", "y"))
      )
    )
  )
  space <- fmridataset::volume_space(
    c(2L, 2L, 2L), support = c(1L, 2L, 4L, 8L), template = "migration-test"
  )
  beta <- matrix(seq_len(16), nrow = 4L)
  fmridataset::fmri_frame(
    assays = list(beta = beta, variance = beta / 10),
    observations = observations,
    features = fmridataset::feature_axis(
      data.frame(.feature_id = fmridataset::feature_ids(space)), space = space
    ),
    active_assay = "beta",
    metadata = list(label = "legacy-fixture"),
    provenance = fmridataset::as_provenance_graph(list(step = "fixture"))
  )
}

.write_provisional_frame <- function(x, path) {
  assays <- lapply(fmridataset::assays(x), function(value) {
    list(
      name = value$name,
      dataset = paste0("/assays/", value$name),
      dtype = value$dtype,
      role = value$role,
      units = value$units,
      metadata = value$metadata
    )
  })
  manifest <- list(
    schema_version = "0.1-provisional",
    observations = x$observations,
    features = x$features,
    entities = x$entities,
    relations = x$relations,
    tables = x$tables,
    active_assay = fmridataset::active_assay(x),
    metadata = x$metadata,
    provenance = x$provenance,
    assays = assays
  )
  h5 <- hdf5r::H5File$new(path, mode = "w")
  on.exit(h5$close_all(), add = TRUE)
  hdf5r::h5attr(h5, "fds_schema_version") <- "0.1-provisional"
  hdf5r::h5attr(h5, "fds_object_type") <- "fmri_frame"
  h5$create_group("assays")
  h5$create_group("manifest")
  for (name in names(assays)) {
    h5$create_dataset(
      paste0("assays/", name),
      robj = fmridataset::collect_assay(x, name),
      chunk_dims = c(2L, 2L)
    )
  }
  h5$create_dataset(
    "manifest/rds",
    robj = as.integer(serialize(manifest, NULL, version = 3L))
  )
  h5$flush()
  invisible(path)
}

test_that("frame store inspection distinguishes current, legacy, and incomplete files", {
  frame <- .recovery_frame()
  current <- tempfile(fileext = ".h5")
  legacy <- tempfile(fileext = ".h5")
  incomplete <- tempfile(fileext = ".h5")
  corrupt <- tempfile(fileext = ".h5")
  on.exit(unlink(c(current, legacy, incomplete, corrupt)), add = TRUE)
  write_frame_h5(frame, current)
  .write_provisional_frame(frame, legacy)
  file.copy(current, incomplete)
  writeLines("not hdf5", corrupt)
  h5 <- hdf5r::H5File$new(incomplete, mode = "r+")
  hdf5r::h5attr(h5, "fds_commit_state") <- "writing"
  h5$close_all()

  expect_identical(frame_store_info(current)$status, "complete")
  legacy_info <- frame_store_info(legacy)
  expect_identical(legacy_info$status, "legacy")
  expect_true(legacy_info$migratable)
  expect_identical(frame_store_info(incomplete)$status, "incomplete")
  expect_identical(frame_store_info(corrupt)$status, "corrupt")
  expect_identical(frame_schema_migrations()$from_schema, "0.1-provisional")
})

test_that("provisional frames migrate lazily and atomically to FDS v1", {
  frame <- .recovery_frame()
  legacy <- tempfile(fileext = ".h5")
  migrated <- tempfile(fileext = ".h5")
  on.exit(unlink(c(legacy, migrated)), add = TRUE)
  .write_provisional_frame(frame, legacy)

  expect_identical(
    migrate_frame_h5(legacy, migrated, memory_budget = 8L),
    normalizePath(migrated)
  )
  expect_identical(frame_store_info(migrated)$status, "complete")
  reopened <- open_frame_h5(migrated)
  expect_identical(fmridataset::observations(reopened), fmridataset::observations(frame))
  expect_identical(
    fmridataset::space_digest(fmridataset::space(reopened)),
    fmridataset::space_digest(fmridataset::space(frame))
  )
  expect_equal(fmridataset::collect_assay(reopened, "beta"), fmridataset::collect_assay(frame, "beta"))
  expect_equal(
    fmridataset::source_read(
      fmridataset::axis_block_data(fmridataset::obs_blocks(reopened)$motion)
    ),
    fmridataset::axis_block_data(fmridataset::obs_blocks(frame)$motion)
  )

  in_place <- tempfile(fileext = ".h5")
  on.exit(unlink(in_place), add = TRUE)
  .write_provisional_frame(frame, in_place)
  before <- unname(tools::md5sum(in_place))
  withr::local_options(fmristore.frame_writer_fault_after_array = 1L)
  expect_error(
    migrate_frame_h5(in_place, in_place, overwrite = TRUE, memory_budget = 8L),
    "Injected"
  )
  expect_identical(unname(tools::md5sum(in_place)), before)
  expect_identical(frame_store_info(in_place)$status, "legacy")
})

test_that("recovery removes only proven uncommitted shard artifacts", {
  frame <- .recovery_frame()
  path <- tempfile("recovery-shards-")
  on.exit(unlink(path, recursive = TRUE), add = TRUE)
  write_sharded_frame_h5(frame[1:2, ], path, shard_id = "sub-01")
  committed <- c(
    file.path(path, "manifest.rds"),
    file.path(path, "shards", "s000001.fds.h5")
  )
  committed_hashes <- unname(tools::md5sum(committed))

  orphan <- file.path(path, "shards", "s000002.fds.h5")
  expect_true(file.copy(committed[[2L]], orphan))
  shard_partial <- file.path(path, "shards", ".shard-interrupted.h5")
  manifest_partial <- file.path(path, ".manifest-interrupted.rds")
  file.create(shard_partial, manifest_partial)
  sibling <- file.path(dirname(path), paste0(".", basename(path), ".partial-interrupted"))
  dir.create(sibling)
  file.create(file.path(sibling, "payload"))
  on.exit(unlink(sibling, recursive = TRUE), add = TRUE)

  plan <- recover_frame_store(path, "report")
  expect_setequal(
    plan$kind,
    c("sibling_partial", "manifest_partial", "shard_partial", "orphan_shard")
  )
  expect_true(all(file.exists(plan$path) | dir.exists(plan$path)))
  expect_false(any(plan$removed))
  expect_identical(unname(tools::md5sum(committed)), committed_hashes)

  rolled_back <- recover_frame_store(path, "rollback")
  expect_true(all(rolled_back$removed))
  expect_false(any(file.exists(rolled_back$path) | dir.exists(rolled_back$path)))
  expect_identical(unname(tools::md5sum(committed)), committed_hashes)
  expect_identical(frame_store_info(path)$status, "complete")
  expect_equal(
    fmridataset::collect_assay(open_sharded_frame_h5(path)),
    fmridataset::collect_assay(frame[1:2, ])
  )
})

test_that("recovery reports missing committed shards without deleting evidence", {
  frame <- .recovery_frame()
  path <- tempfile("recovery-shards-")
  on.exit(unlink(path, recursive = TRUE), add = TRUE)
  write_sharded_frame_h5(frame[1:2, ], path, shard_id = "sub-01")
  shard <- file.path(path, "shards", "s000001.fds.h5")
  moved <- paste0(shard, ".held")
  file.rename(shard, moved)
  on.exit(if (file.exists(moved)) file.rename(moved, shard), add = TRUE)

  info <- frame_store_info(path)
  expect_identical(info$status, "incomplete")
  expect_identical(info$recovery$kind, "missing_committed_shard")
  result <- recover_frame_store(path, "rollback")
  expect_false(result$removable)
  expect_false(result$removed)
  expect_true(file.exists(moved))
})
