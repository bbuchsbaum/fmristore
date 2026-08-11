.example_h5_frame <- function() {
  observations <- data.frame(
    .obs_id = sprintf("obs-%02d", 1:6),
    subject_id = factor(rep(c("s01", "s02", "s03"), each = 2L)),
    condition = factor(c("A", "B", "A", "A", "B", "A"), levels = c("A", "B")),
    stringsAsFactors = FALSE
  )
  observations <- fmridataset::axis_frame(
    observations,
    blocks = list(
      embedding = fmridataset::axis_block(
        matrix(seq_len(18), nrow = 6L),
        components = data.frame(
          .component_id = paste0("pc", 1:3),
          name = paste0("PC", 1:3)
        )
      )
    )
  )
  space <- fmridataset::volume_space(
    dim = c(2L, 2L, 2L),
    affine = diag(4),
    support = c(1L, 3L, 4L, 7L),
    template = "test-space"
  )
  features <- fmridataset::feature_axis(
    data.frame(.feature_id = fmridataset::feature_ids(space)),
    space = space
  )
  beta <- matrix(seq_len(24), nrow = 6L, ncol = 4L)
  fmridataset::fmri_frame(
    assays = list(beta = beta, variance = beta / 10 + 0.1),
    observations = observations,
    features = features,
    entities = list(subject = data.frame(subject_id = c("s01", "s02", "s03"))),
    relations = list(observation_subject = observations$data$subject_id),
    tables = list(contrasts = data.frame(name = c("A", "B"))),
    active_assay = "beta",
    metadata = list(label = "roundtrip"),
    provenance = list(step = "fixture")
  )
}

test_that("provisional HDF5 frames round trip without realizing on open", {
  frame <- .example_h5_frame()
  path <- tempfile(fileext = ".fds.h5")
  on.exit(unlink(path), add = TRUE)

  expect_identical(write_frame_h5(frame, path), normalizePath(path, mustWork = FALSE))
  reopened <- open_frame_h5(path)

  expect_s3_class(reopened, "fmri_frame")
  expect_identical(fmridataset::observation_ids(reopened), fmridataset::observation_ids(frame))
  expect_identical(fmridataset::feature_ids(reopened), fmridataset::feature_ids(frame))
  expect_identical(
    fmridataset::space_digest(fmridataset::space(reopened)),
    fmridataset::space_digest(fmridataset::space(frame))
  )
  expect_identical(fmridataset::observations(reopened), fmridataset::observations(frame))
  expect_identical(
    fmridataset::block_components(fmridataset::obs_blocks(reopened)$embedding),
    fmridataset::block_components(fmridataset::obs_blocks(frame)$embedding)
  )
  expect_s3_class(fmridataset::assay(reopened, "beta")$source, "h5_array_source")
  expect_equal(fmridataset::collect_assay(reopened, "beta"), fmridataset::collect_assay(frame, "beta"))
  expect_equal(
    fmridataset::collect_assay(reopened, "variance"),
    fmridataset::collect_assay(frame, "variance")
  )
  expect_equal(
    neuroim2::values(fmridataset::spatial_map(reopened, "obs-01")),
    neuroim2::values(fmridataset::spatial_map(frame, "obs-01"))
  )
})

test_that("frame writes are bounded and atomically fail", {
  frame <- .example_h5_frame()
  parent <- tempfile("frame-write-")
  dir.create(parent)
  on.exit(unlink(parent, recursive = TRUE), add = TRUE)
  path <- file.path(parent, "study.fds.h5")

  expect_error(write_frame_h5(frame, path, memory_budget = 8L), "memory_budget")
  expect_false(file.exists(path))
  expect_length(list.files(parent, pattern = "partial", all.files = TRUE), 0L)

  withr::local_options(fmristore.frame_writer_fault_after_assay = 1L)
  expect_error(write_frame_h5(frame, path), "Injected")
  expect_false(file.exists(path))
  expect_length(list.files(parent, pattern = "partial", all.files = TRUE), 0L)
})

test_that("frame writers protect existing destinations and schema", {
  frame <- .example_h5_frame()
  path <- tempfile(fileext = ".fds.h5")
  on.exit(unlink(path), add = TRUE)
  write_frame_h5(frame, path)

  expect_error(write_frame_h5(frame, path), "already exists")
  expect_silent(write_frame_h5(frame, path, overwrite = TRUE))

  invalid <- tempfile(fileext = ".h5")
  on.exit(unlink(invalid), add = TRUE)
  h5 <- hdf5r::H5File$new(invalid, mode = "w")
  h5[["data"]] <- 1
  h5$close_all()
  expect_error(open_frame_h5(invalid), "schema")
})
