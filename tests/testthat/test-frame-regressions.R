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
