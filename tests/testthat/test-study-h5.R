.study_h5_frame <- function(prefix, subjects, space = NULL) {
  if (is.null(space)) {
    space <- fmridataset::index_space(
      3L, ids = paste0(prefix, "-feature-", 1:3), namespace = prefix
    )
  }
  observations <- tibble::tibble(
    .obs_id = paste0(prefix, "-", seq_along(subjects)),
    subject_id = subjects
  )
  fmridataset::fmri_frame(
    assays = list(
      beta = matrix(seq_len(length(subjects) * fmridataset::n_features(space)),
                    nrow = length(subjects)),
      variance = matrix(seq_len(length(subjects) * fmridataset::n_features(space)) / 10,
                        nrow = length(subjects))
    ),
    observations = observations,
    space = space,
    entities = list(subject = fmridataset::entity_frame(
      tibble::tibble(subject_id = unique(subjects)), key = "subject_id"
    )),
    relations = list(observation_subject = fmridataset::key_relation(
      "subject_id", target = "subject"
    )),
    active_assay = "beta"
  )
}

.study_h5_fixture <- function() {
  phenotype <- fmridataset::counting_source(fmridataset::memory_source(
    matrix(c(.2, .8, .4, .6), nrow = 2L)
  ))
  entities <- list(subject = fmridataset::entity_frame(
    tibble::tibble(subject_id = c("sub-01", "sub-02"), age = c(63, 71)),
    key = "subject_id",
    blocks = list(phenotype = fmridataset::axis_block(
      phenotype,
      components = tibble::tibble(.component_id = c("memory", "attention"))
    ))
  ))
  bold <- .study_h5_frame("bold", c("sub-01", "sub-02"))
  beta <- .study_h5_frame("beta", c("sub-01", "sub-02"))
  native <- fmridataset::fmri_collection(list(
    sub_01 = .study_h5_frame(
      "native-01", "sub-01",
      fmridataset::volume_space(
        c(2L, 2L, 2L), support = 1:3, template = "sub-01"
      )
    ),
    sub_02 = .study_h5_frame(
      "native-02", "sub-02",
      fmridataset::volume_space(
        c(3L, 2L, 2L), support = c(1L, 2L, 4L, 5L), template = "sub-02"
      )
    )
  ), metadata = list(domain = "native"))
  study <- fmridataset::fmri_study(
    frames = list(bold = bold, beta = beta, native = native),
    entities = entities,
    links = list(beta_from_bold = fmridataset::frame_link(
      "beta", "bold", "derived_from",
      map = tibble::tibble(
        .from_id = fmridataset::observation_ids(beta),
        .to_id = fmridataset::observation_ids(bold)
      )
    )),
    tables = list(events = fmridataset::event_table(tibble::tibble(
      event_id = c("event-1", "event-2"),
      subject_id = c("sub-01", "sub-02"),
      onset = c(0, 2), duration = c(1, 1)
    ))),
    metadata = list(title = "atomic study"),
    provenance = fmridataset::as_provenance_graph(list(step = "fixture"))
  )
  list(study = study, phenotype = phenotype)
}

test_that("HDF5 studies round trip frames collections entities links and tables lazily", {
  fixture <- .study_h5_fixture()
  path <- tempfile("study-store-")
  on.exit(unlink(path, recursive = TRUE), add = TRUE)

  expect_identical(
    write_study_h5(fixture$study, path, memory_budget = 16L),
    normalizePath(path, mustWork = FALSE)
  )
  after_write <- fmridataset::source_counts(fixture$phenotype)
  reopened <- open_study_h5(path)

  expect_s3_class(reopened, "fmri_study")
  expect_identical(
    fmridataset::fds_study_manifest(reopened),
    fmridataset::fds_study_manifest(fixture$study)
  )
  expect_identical(fmridataset::study_ids(reopened), c("bold", "beta", "native"))
  expect_s3_class(fmridataset::study_frame(reopened, "native"), "fmri_collection")
  expect_identical(fmridataset::study_links(reopened), fmridataset::study_links(fixture$study))
  expect_identical(fmridataset::study_tables(reopened), fmridataset::study_tables(fixture$study))
  expect_identical(reopened$metadata, fixture$study$metadata)
  expect_identical(reopened$provenance, fixture$study$provenance)
  expect_identical(fmridataset::source_counts(fixture$phenotype), after_write)
  expect_s3_class(unserialize(serialize(reopened, NULL)), "fmri_study")

  for (name in c("bold", "beta")) {
    expect_equal(
      fmridataset::collect_assay(fmridataset::study_frame(reopened, name), "beta"),
      fmridataset::collect_assay(fmridataset::study_frame(fixture$study, name), "beta")
    )
  }
  reopened_block <- fmridataset::entity_blocks(
    fmridataset::entity(reopened, "subject")
  )$phenotype
  expect_s3_class(fmridataset::axis_block_data(reopened_block), "h5_array_source")
  expect_equal(
    fmridataset::source_read(fmridataset::as_array_source(
      fmridataset::axis_block_data(reopened_block)
    )),
    matrix(c(.2, .8, .4, .6), nrow = 2L)
  )
})

test_that("study writes fail before visibility and clean all temporary artifacts", {
  fixture <- .study_h5_fixture()
  parent <- tempfile("study-atomic-")
  dir.create(parent)
  on.exit(unlink(parent, recursive = TRUE), add = TRUE)
  path <- file.path(parent, "study.fds")

  withr::local_options(fmristore.study_writer_fault_after_representation = 1L)
  expect_error(write_study_h5(fixture$study, path), "Injected")
  expect_false(file.exists(path))
  expect_length(list.files(parent, pattern = "partial", all.files = TRUE), 0L)

  withr::local_options(fmristore.study_writer_fault_after_representation = Inf)
  withr::local_options(fmristore.study_writer_fault_after_array = 1L)
  expect_error(write_study_h5(fixture$study, path), "Injected")
  expect_false(file.exists(path))
  expect_length(list.files(parent, pattern = "partial", all.files = TRUE), 0L)
})

test_that("failed study overwrite restores the previous complete destination", {
  fixture <- .study_h5_fixture()
  parent <- tempfile("study-overwrite-")
  dir.create(parent)
  on.exit(unlink(parent, recursive = TRUE), add = TRUE)
  path <- file.path(parent, "study.fds")
  write_study_h5(fixture$study, path)
  before <- unname(tools::md5sum(file.path(path, "manifest.rds")))

  replacement <- fixture$study
  replacement$metadata$title <- "replacement"
  withr::local_options(fmristore.study_writer_fault_after_backup = TRUE)
  expect_error(write_study_h5(replacement, path, overwrite = TRUE), "Injected")

  expect_identical(unname(tools::md5sum(file.path(path, "manifest.rds"))), before)
  expect_identical(open_study_h5(path)$metadata$title, "atomic study")
  expect_length(list.files(parent, pattern = "partial|backup", all.files = TRUE), 0L)

  withr::local_options(fmristore.study_writer_fault_after_backup = FALSE)
  expect_silent(write_study_h5(replacement, path, overwrite = TRUE))
  expect_identical(open_study_h5(path)$metadata$title, "replacement")
  expect_length(list.files(parent, pattern = "partial|backup", all.files = TRUE), 0L)
})

test_that("study stores reject incomplete roots and semantic substitution", {
  fixture <- .study_h5_fixture()
  parent <- tempfile("study-invalid-")
  dir.create(parent)
  on.exit(unlink(parent, recursive = TRUE), add = TRUE)

  expect_error(open_study_h5(parent), "manifest")

  path <- file.path(parent, "valid.fds")
  write_study_h5(fixture$study, path)
  root <- readRDS(file.path(path, "manifest.rds"))
  root$semantic$metadata$title <- "tampered"
  saveRDS(root, file.path(path, "manifest.rds"))
  expect_error(open_study_h5(path), "digest")
})

test_that("filtered study views persist as compact self-contained studies", {
  fixture <- .study_h5_fixture()
  selected <- fmridataset::filter_entities(fixture$study, subject, age >= 65)
  path <- tempfile("study-filtered-")
  on.exit(unlink(path, recursive = TRUE), add = TRUE)

  write_study_h5(selected, path)
  reopened <- open_study_h5(path)

  expect_identical(
    fmridataset::entity_ids(fmridataset::entity(reopened, "subject")),
    "sub-02"
  )
  expect_identical(nrow(fmridataset::study_frame(reopened, "bold")), 1L)
  expect_identical(nrow(fmridataset::event_data(fmridataset::events(reopened))), 1L)
})

test_that("studies without shared blocks omit the study array file", {
  fixture <- .study_h5_fixture()
  subjects <- fmridataset::entity_data(fmridataset::entity(fixture$study, "subject"))
  study <- fmridataset::fmri_study(
    frames = fmridataset::study_frames(fixture$study, contextual = FALSE),
    entities = list(subject = fmridataset::entity_frame(subjects, key = "subject_id")),
    links = fmridataset::study_links(fixture$study),
    tables = fmridataset::study_tables(fixture$study),
    metadata = fixture$study$metadata,
    provenance = fixture$study$provenance
  )
  path <- tempfile("study-no-shared-arrays-")
  on.exit(unlink(path, recursive = TRUE), add = TRUE)

  write_study_h5(study, path)
  expect_false(file.exists(file.path(path, "shared-arrays.fds.h5")))
  reopened <- open_study_h5(path)
  expect_identical(fmridataset::fds_study_manifest(reopened), fmridataset::fds_study_manifest(study))
})
