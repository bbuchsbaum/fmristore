.frame_shard_schema_id <- "org.fmridataset.fds-shards/v1"
.frame_shard_schema_version <- 1L

.validate_frame_shard_id <- function(shard_id) {
  if (!is.character(shard_id) || length(shard_id) != 1L || is.na(shard_id) ||
    !grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", shard_id)) {
    stop(
      "`shard_id` must start with a letter or digit and contain only letters, digits, dot, underscore, or hyphen.",
      call. = FALSE
    )
  }
  shard_id
}

.frame_shard_contract <- function(x) {
  manifest <- fmridataset::fds_frame_manifest(x)
  observation_columns <- lapply(manifest$axes$observation$data, function(column) {
    list(
      class = class(column),
      typeof = typeof(column),
      levels = if (is.factor(column)) levels(column) else NULL,
      ordered = is.ordered(column)
    )
  })
  arrays <- lapply(manifest$arrays, function(array) {
    if (identical(array$axes[[1L]], "observation")) array$shape[[1L]] <- NA_integer_
    array
  })
  assays <- lapply(manifest$assays, function(assay) {
    assay$shape[[1L]] <- NA_integer_
    assay$observation_digest <- NULL
    assay
  })
  list(
    feature = manifest$axes$feature,
    observation_columns = observation_columns,
    observation_blocks = manifest$axes$observation$blocks,
    observation_metadata = manifest$axes$observation$metadata,
    arrays = arrays,
    assays = assays,
    entities = manifest$entities,
    relations = manifest$relations,
    tables = manifest$tables,
    active_assay = manifest$active_assay,
    metadata = manifest$metadata,
    provenance = manifest$provenance,
    extensions = manifest$extensions
  )
}

.frame_shard_root <- function(path) file.path(path, "manifest.rds")

.read_frame_shard_root <- function(path) {
  if (!is.character(path) || length(path) != 1L || !dir.exists(path)) {
    stop("Sharded HDF5 frame directory does not exist.", call. = FALSE)
  }
  root <- .frame_shard_root(path)
  if (!file.exists(root)) stop("Sharded HDF5 frame manifest is missing.", call. = FALSE)
  manifest <- readRDS(root)
  required <- c("schema_id", "schema_version", "object_type", "contract", "shards")
  if (!is.list(manifest) || !identical(names(manifest), required) ||
    !identical(manifest$schema_id, .frame_shard_schema_id) ||
    !identical(manifest$schema_version, .frame_shard_schema_version) ||
    !identical(manifest$object_type, "fmri_frame_shards") ||
    !is.data.frame(manifest$shards) || !nrow(manifest$shards)) {
    stop("Unsupported or invalid sharded HDF5 frame manifest.", call. = FALSE)
  }
  expected <- c(".shard_id", ".file", ".n_observation", ".semantic_digest")
  expected_files <- file.path(
    "shards",
    paste0("s", sprintf("%06d", seq_len(nrow(manifest$shards))), ".fds.h5")
  )
  if (!identical(names(manifest$shards), expected) ||
    !is.character(manifest$shards$.shard_id) ||
    anyNA(manifest$shards$.shard_id) ||
    any(!grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", manifest$shards$.shard_id)) ||
    anyDuplicated(manifest$shards$.shard_id) ||
    !identical(manifest$shards$.file, expected_files) ||
    !is.integer(manifest$shards$.n_observation) ||
    anyNA(manifest$shards$.n_observation) ||
    any(manifest$shards$.n_observation < 1L) ||
    !is.character(manifest$shards$.semantic_digest) ||
    anyNA(manifest$shards$.semantic_digest) ||
    any(!nzchar(manifest$shards$.semantic_digest)) ||
    !is.list(manifest$contract)) {
    stop("Invalid sharded HDF5 frame entries.", call. = FALSE)
  }
  manifest
}

.write_frame_shard_root <- function(path, manifest) {
  temporary <- tempfile(pattern = ".manifest-", tmpdir = path, fileext = ".rds")
  committed <- FALSE
  on.exit(if (!committed && file.exists(temporary)) unlink(temporary), add = TRUE)
  saveRDS(manifest, temporary, version = 3L)
  if (!file.rename(temporary, .frame_shard_root(path))) {
    stop("Could not atomically commit the sharded frame manifest.", call. = FALSE)
  }
  committed <- TRUE
  invisible(manifest)
}

.frame_shard_entry <- function(x, shard_id, file) {
  semantic <- fmridataset::fds_frame_manifest(x)
  data.frame(
    .shard_id = shard_id,
    .file = file,
    .n_observation = as.integer(nrow(x)),
    .semantic_digest = fmridataset::fds_manifest_digest(semantic),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

#' Write, append, and open row-sharded HDF5 frames
#'
#' A sharded frame is a directory containing immutable FDS v1 HDF5 frame
#' files and one small root manifest. Appending writes a new shard and then
#' atomically replaces only the root manifest; previously committed HDF5 files
#' are never opened for writing or rewritten.
#'
#' All shards must have identical feature identity, assay declarations,
#' observation-block definitions, and non-observation registries. Observation
#' IDs must be globally unique. Numerical arrays remain lazy when reopened.
#'
#' @param x An `fmri_frame` containing one new observation shard.
#' @param path Sharded frame directory.
#' @param shard_id Stable unique shard identifier.
#' @param chunk_dims,layout,target_chunk_bytes,memory_budget Passed to
#'   [write_frame_h5()].
#' @return Writers invisibly return the normalized directory path.
#'   `open_sharded_frame_h5()` returns one lazy `fmri_frame`, and
#'   `frame_shard_manifest()` returns root metadata without opening a shard.
#' @export
write_sharded_frame_h5 <- function(
    x,
    path,
    shard_id = "shard-000001",
    chunk_dims = NULL,
    layout = c("balanced", "imagewise", "featurewise"),
    target_chunk_bytes = getOption("fmristore.target_chunk_bytes", 1024^2),
    memory_budget = getOption("fmridataset.block_budget", 512 * 1024^2)) {
  if (!inherits(x, "fmri_frame")) stop("`x` must be an fmri_frame.", call. = FALSE)
  shard_id <- .validate_frame_shard_id(shard_id)
  layout <- match.arg(layout)
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  if (file.exists(path) || dir.exists(path)) {
    stop("Sharded frame destination already exists.", call. = FALSE)
  }
  parent <- dirname(path)
  if (!dir.exists(parent)) stop("Destination directory does not exist.", call. = FALSE)
  temporary <- tempfile(pattern = paste0(".", basename(path), ".partial-"), tmpdir = parent)
  dir.create(file.path(temporary, "shards"), recursive = TRUE)
  committed <- FALSE
  on.exit(if (!committed && dir.exists(temporary)) unlink(temporary, recursive = TRUE), add = TRUE)

  relative <- file.path("shards", "s000001.fds.h5")
  write_frame_h5(
    x, file.path(temporary, relative), chunk_dims = chunk_dims, layout = layout,
    target_chunk_bytes = target_chunk_bytes, memory_budget = memory_budget
  )
  manifest <- list(
    schema_id = .frame_shard_schema_id,
    schema_version = .frame_shard_schema_version,
    object_type = "fmri_frame_shards",
    contract = .frame_shard_contract(x),
    shards = .frame_shard_entry(x, shard_id, relative)
  )
  saveRDS(manifest, .frame_shard_root(temporary), version = 3L)
  if (!file.rename(temporary, path)) {
    stop("Could not atomically commit sharded HDF5 frame directory.", call. = FALSE)
  }
  committed <- TRUE
  invisible(normalizePath(path, winslash = "/", mustWork = TRUE))
}

#' @rdname write_sharded_frame_h5
#' @export
append_frame_shard_h5 <- function(
    x,
    path,
    shard_id,
    chunk_dims = NULL,
    layout = c("balanced", "imagewise", "featurewise"),
    target_chunk_bytes = getOption("fmristore.target_chunk_bytes", 1024^2),
    memory_budget = getOption("fmridataset.block_budget", 512 * 1024^2)) {
  if (!inherits(x, "fmri_frame")) stop("`x` must be an fmri_frame.", call. = FALSE)
  shard_id <- .validate_frame_shard_id(shard_id)
  layout <- match.arg(layout)
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  manifest <- .read_frame_shard_root(path)
  if (shard_id %in% manifest$shards$.shard_id) {
    stop("Sharded frame already contains `shard_id`: ", shard_id, call. = FALSE)
  }
  if (!identical(.frame_shard_contract(x), manifest$contract)) {
    stop("Appended frame is incompatible with the committed shard contract.", call. = FALSE)
  }
  existing <- open_sharded_frame_h5(path)
  if (any(fmridataset::observation_ids(x) %in% fmridataset::observation_ids(existing))) {
    stop("Observation IDs collide with the committed sharded frame.", call. = FALSE)
  }

  index <- nrow(manifest$shards) + 1L
  relative <- file.path("shards", paste0("s", sprintf("%06d", index), ".fds.h5"))
  final <- file.path(path, relative)
  if (file.exists(final)) stop("Target shard file already exists.", call. = FALSE)
  temporary <- tempfile(pattern = ".shard-", tmpdir = file.path(path, "shards"), fileext = ".h5")
  committed <- FALSE
  on.exit({
    if (!committed && file.exists(temporary)) unlink(temporary)
    if (!committed && file.exists(final)) unlink(final)
  }, add = TRUE)
  write_frame_h5(
    x, temporary, chunk_dims = chunk_dims, layout = layout,
    target_chunk_bytes = target_chunk_bytes, memory_budget = memory_budget
  )
  if (!file.rename(temporary, final)) stop("Could not commit new frame shard.", call. = FALSE)
  if (isTRUE(getOption("fmristore.shard_writer_fault_after_file", FALSE))) {
    stop("Injected sharded frame append failure.", call. = FALSE)
  }

  manifest$shards <- rbind(
    manifest$shards,
    .frame_shard_entry(x, shard_id, relative)
  )
  row.names(manifest$shards) <- seq_len(nrow(manifest$shards))
  .write_frame_shard_root(path, manifest)
  committed <- TRUE
  invisible(path)
}

#' @rdname write_sharded_frame_h5
#' @export
frame_shard_manifest <- function(path) .read_frame_shard_root(path)$shards

#' @rdname write_sharded_frame_h5
#' @export
open_sharded_frame_h5 <- function(path) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  manifest <- .read_frame_shard_root(path)
  frames <- lapply(seq_len(nrow(manifest$shards)), function(index) {
    entry <- manifest$shards[index, , drop = FALSE]
    shard_path <- file.path(path, entry$.file)
    if (!file.exists(shard_path)) stop("Committed frame shard is missing: ", entry$.file, call. = FALSE)
    frame <- open_frame_h5(shard_path)
    semantic <- fmridataset::fds_frame_manifest(frame)
    if (!identical(fmridataset::fds_manifest_digest(semantic), entry$.semantic_digest) ||
      !identical(.frame_shard_contract(frame), manifest$contract) ||
      nrow(frame) != entry$.n_observation) {
      stop("Committed frame shard does not match the root manifest: ", entry$.shard_id, call. = FALSE)
    }
    frame
  })
  combined <- do.call(fmridataset::bind_observations, frames)
  shard_ids <- manifest$shards$.shard_id
  for (name in names(fmridataset::assays(combined))) {
    combined$assays[[name]]$source <- fmridataset::row_sharded_source(
      lapply(frames, function(frame) fmridataset::assay(frame, name)$source),
      shard_ids = shard_ids
    )
  }
  for (name in names(fmridataset::obs_blocks(combined))) {
    combined$observations$blocks[[name]]$data <- fmridataset::row_sharded_source(
      lapply(frames, function(frame) {
        fmridataset::axis_block_data(fmridataset::obs_blocks(frame)[[name]])
      }),
      shard_ids = shard_ids
    )
  }
  combined
}
