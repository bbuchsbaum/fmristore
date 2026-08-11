.study_h5_schema_id <- "org.fmridataset.fds-study-store/v1"
.study_h5_schema_version <- 1L

.study_store_root <- function(path) file.path(path, "manifest.rds")

.study_path_exists <- function(path) file.exists(path) || dir.exists(path)

.study_representation_file <- function(index) {
  file.path("representations", paste0("r", sprintf("%06d", index), ".fds.h5"))
}

.study_collection_member_file <- function(representation_index, member_index) {
  file.path(
    "representations",
    paste0("r", sprintf("%06d", representation_index)),
    paste0("m", sprintf("%06d", member_index), ".fds.h5")
  )
}

.study_frame_storage_entry <- function(frame, relative) {
  list(
    file = relative,
    semantic_digest = fmridataset::fds_manifest_digest(
      fmridataset::fds_frame_manifest(frame)
    )
  )
}

.write_study_representations <- function(
    representations,
    semantic,
    path,
    chunk_dims,
    layout,
    target_chunk_bytes,
    memory_budget) {
  dir.create(file.path(path, "representations"), recursive = TRUE)
  storage <- vector("list", length(representations))
  names(storage) <- names(representations)
  fail_after <- getOption("fmristore.study_writer_fault_after_representation", Inf)
  for (index in seq_along(representations)) {
    name <- names(representations)[[index]]
    value <- representations[[name]]
    descriptor <- semantic$representations[[name]]
    if (identical(descriptor$type, "fmri_frame")) {
      relative <- .study_representation_file(index)
      write_frame_h5(
        value, file.path(path, relative),
        chunk_dims = chunk_dims, layout = layout,
        target_chunk_bytes = target_chunk_bytes,
        memory_budget = memory_budget
      )
      storage[[name]] <- c(
        list(type = "fmri_frame"),
        .study_frame_storage_entry(value, relative)
      )
    } else {
      member_values <- fmridataset::collection_frames(value)
      member_storage <- vector("list", length(member_values))
      names(member_storage) <- names(member_values)
      for (member_index in seq_along(member_values)) {
        member_name <- names(member_values)[[member_index]]
        relative <- .study_collection_member_file(index, member_index)
        dir.create(dirname(file.path(path, relative)), recursive = TRUE, showWarnings = FALSE)
        write_frame_h5(
          member_values[[member_name]], file.path(path, relative),
          chunk_dims = chunk_dims, layout = layout,
          target_chunk_bytes = target_chunk_bytes,
          memory_budget = memory_budget
        )
        member_storage[[member_name]] <- .study_frame_storage_entry(
          member_values[[member_name]], relative
        )
      }
      storage[[name]] <- list(
        type = "fmri_collection",
        members = member_storage
      )
    }
    if (is.finite(fail_after) && index >= fail_after) {
      stop("Injected HDF5 study writer failure after representation.", call. = FALSE)
    }
  }
  storage
}

.write_study_arrays <- function(
    semantic,
    bindings,
    path,
    chunk_dims,
    layout,
    target_chunk_bytes,
    memory_budget) {
  if (!length(semantic$arrays)) return(list(file = NULL, arrays = list()))
  relative <- "shared-arrays.fds.h5"
  h5 <- hdf5r::H5File$new(file.path(path, relative), mode = "w")
  on.exit(close_h5_safely(h5), add = TRUE)
  hdf5r::h5attr(h5, "fds_schema_id") <- .study_h5_schema_id
  hdf5r::h5attr(h5, "fds_schema_version") <- .study_h5_schema_version
  hdf5r::h5attr(h5, "fds_object_type") <- "fmri_study_arrays"
  hdf5r::h5attr(h5, "fds_commit_state") <- "writing"
  semantic_digest <- fmridataset::fds_study_manifest_digest(semantic)
  hdf5r::h5attr(h5, "fds_semantic_digest") <- semantic_digest
  h5$create_group("arrays")
  storage <- vector("list", length(semantic$arrays))
  names(storage) <- names(semantic$arrays)
  fail_after <- getOption("fmristore.study_writer_fault_after_array", Inf)
  for (index in seq_along(semantic$arrays)) {
    key <- names(semantic$arrays)[[index]]
    storage[[key]] <- .write_frame_array(
      h5,
      binding = bindings[[key]],
      declaration = semantic$arrays[[key]],
      dataset_path = paste0("arrays/a", sprintf("%06d", index)),
      memory_budget = memory_budget,
      chunk_dims = chunk_dims,
      layout = layout,
      target_chunk_bytes = target_chunk_bytes
    )
    if (is.finite(fail_after) && index >= fail_after) {
      stop("Injected HDF5 study writer failure after shared array.", call. = FALSE)
    }
  }
  hdf5r::h5attr(h5, "fds_commit_state") <- "complete"
  h5$flush()
  close_h5_safely(h5)
  list(file = relative, arrays = storage)
}

.commit_study_directory <- function(temporary, path, overwrite) {
  if (!.study_path_exists(path)) {
    if (!file.rename(temporary, path)) {
      stop("Could not atomically commit HDF5 study: ", path, call. = FALSE)
    }
    return(invisible(TRUE))
  }
  if (!isTRUE(overwrite)) {
    stop("Destination already exists; set `overwrite = TRUE`.", call. = FALSE)
  }
  backup <- tempfile(
    pattern = paste0(".", basename(path), ".backup-"),
    tmpdir = dirname(path)
  )
  old_moved <- FALSE
  new_committed <- FALSE
  on.exit({
    if (old_moved && !new_committed && !.study_path_exists(path) &&
        .study_path_exists(backup)) {
      file.rename(backup, path)
    }
    if (new_committed && .study_path_exists(backup)) {
      unlink(backup, recursive = TRUE)
    }
  }, add = TRUE)
  if (!file.rename(path, backup)) {
    stop("Could not stage the existing HDF5 study for replacement.", call. = FALSE)
  }
  old_moved <- TRUE
  if (isTRUE(getOption("fmristore.study_writer_fault_after_backup", FALSE))) {
    stop("Injected HDF5 study writer failure after backup.", call. = FALSE)
  }
  if (!file.rename(temporary, path)) {
    stop("Could not atomically install the replacement HDF5 study.", call. = FALSE)
  }
  new_committed <- TRUE
  invisible(TRUE)
}

#' Write and open atomic FDS HDF5 studies
#'
#' `write_study_h5()` persists a canonical `fmri_study` as a directory of
#' immutable FDS HDF5 frames plus one source-free semantic manifest. Shared
#' entity blocks are written once to a study-level HDF5 array file. The root
#' directory becomes visible only after every representation and manifest is
#' complete. Existing destinations are restored if replacement fails.
#'
#' `open_study_h5()` validates semantic and physical digests while opening only
#' metadata. Assays and shared blocks remain lazy `h5_array_source` objects.
#'
#' @param x An `fmri_study` or filtered study view.
#' @param path Destination or source study directory.
#' @param overwrite Atomically replace an existing complete destination.
#' @param chunk_dims Optional two-element HDF5 chunk shape.
#' @param layout Workload profile passed to [write_frame_h5()].
#' @param target_chunk_bytes Desired physical chunk size.
#' @param memory_budget Hard maximum bytes for one source-read block.
#' @return `write_study_h5()` invisibly returns the normalized committed path;
#'   `open_study_h5()` returns a lazy `fmri_study`.
#' @export
write_study_h5 <- function(
    x,
    path,
    overwrite = FALSE,
    chunk_dims = NULL,
    layout = c("balanced", "imagewise", "featurewise"),
    target_chunk_bytes = getOption("fmristore.target_chunk_bytes", 1024^2),
    memory_budget = getOption("fmridataset.block_budget", 512 * 1024^2)) {
  if (!inherits(x, "fmri_study")) {
    stop("`x` must be an fmri_study or fmri_study_view.", call. = FALSE)
  }
  fmridataset::validate_fmri_study(x)
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("`path` must be one non-empty directory path.", call. = FALSE)
  }
  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite)) {
    stop("`overwrite` must be TRUE or FALSE.", call. = FALSE)
  }
  memory_budget <- .validate_frame_budget(memory_budget)
  layout <- match.arg(layout)
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  parent <- dirname(path)
  if (!dir.exists(parent)) {
    stop("Destination directory does not exist: ", parent, call. = FALSE)
  }
  if (.study_path_exists(path) && !isTRUE(overwrite)) {
    stop("Destination already exists; set `overwrite = TRUE`.", call. = FALSE)
  }

  semantic <- fmridataset::fds_study_manifest(x)
  representations <- fmridataset::fds_study_representations(x)
  bindings <- fmridataset::fds_study_bindings(x)
  temporary <- tempfile(
    pattern = paste0(".", basename(path), ".partial-"),
    tmpdir = parent
  )
  dir.create(temporary)
  committed <- FALSE
  on.exit({
    if (!committed && dir.exists(temporary)) unlink(temporary, recursive = TRUE)
  }, add = TRUE)

  representation_storage <- .write_study_representations(
    representations, semantic, temporary, chunk_dims, layout,
    target_chunk_bytes, memory_budget
  )
  shared_storage <- .write_study_arrays(
    semantic, bindings, temporary, chunk_dims, layout,
    target_chunk_bytes, memory_budget
  )
  root <- list(
    schema_id = .study_h5_schema_id,
    schema_version = .study_h5_schema_version,
    object_type = "fmri_study",
    commit_state = "complete",
    semantic_digest = fmridataset::fds_study_manifest_digest(semantic),
    semantic = semantic,
    representations = representation_storage,
    shared_arrays = shared_storage
  )
  saveRDS(root, .study_store_root(temporary), version = 3L)
  .commit_study_directory(temporary, path, overwrite)
  committed <- TRUE
  invisible(normalizePath(path, winslash = "/", mustWork = TRUE))
}

.read_study_root <- function(path) {
  if (!is.character(path) || length(path) != 1L || !dir.exists(path)) {
    stop("HDF5 study directory does not exist.", call. = FALSE)
  }
  root_path <- .study_store_root(path)
  if (!file.exists(root_path)) stop("HDF5 study manifest is missing.", call. = FALSE)
  root <- readRDS(root_path)
  required <- c(
    "schema_id", "schema_version", "object_type", "commit_state",
    "semantic_digest", "semantic", "representations", "shared_arrays"
  )
  if (!is.list(root) || !identical(names(root), required) ||
      !identical(root$schema_id, .study_h5_schema_id) ||
      !identical(root$schema_version, .study_h5_schema_version) ||
      !identical(root$object_type, "fmri_study") ||
      !identical(root$commit_state, "complete")) {
    stop("Unsupported or incomplete HDF5 study manifest.", call. = FALSE)
  }
  fmridataset::validate_fds_study_manifest(root$semantic)
  digest <- fmridataset::fds_study_manifest_digest(root$semantic)
  if (!identical(root$semantic_digest, digest)) {
    stop("HDF5 study semantic digest does not match its manifest.", call. = FALSE)
  }
  if (!is.list(root$representations) ||
      !identical(names(root$representations), names(root$semantic$representations)) ||
      !is.list(root$shared_arrays)) {
    stop("HDF5 study physical manifest does not match its semantic manifest.", call. = FALSE)
  }
  root
}

.validate_study_relative_file <- function(actual, expected, field) {
  if (!is.character(actual) || length(actual) != 1L || is.na(actual) ||
      !identical(actual, expected)) {
    stop("Invalid HDF5 study physical path: ", field, call. = FALSE)
  }
  actual
}

.open_study_frame <- function(path, storage, semantic, expected_file, field) {
  required <- c("file", "semantic_digest")
  if (!is.list(storage) || !identical(names(storage), required)) {
    stop("Invalid HDF5 study frame descriptor: ", field, call. = FALSE)
  }
  relative <- .validate_study_relative_file(storage$file, expected_file, field)
  frame_path <- file.path(path, relative)
  if (!file.exists(frame_path)) stop("Committed HDF5 study frame is missing: ", relative, call. = FALSE)
  frame <- open_frame_h5(frame_path)
  digest <- fmridataset::fds_manifest_digest(fmridataset::fds_frame_manifest(frame))
  expected_digest <- fmridataset::fds_manifest_digest(semantic)
  if (!identical(storage$semantic_digest, expected_digest) ||
      !identical(digest, expected_digest)) {
    stop("HDF5 study frame digest mismatch: ", field, call. = FALSE)
  }
  frame
}

.open_study_representations <- function(path, root) {
  out <- vector("list", length(root$semantic$representations))
  names(out) <- names(root$semantic$representations)
  for (index in seq_along(out)) {
    name <- names(out)[[index]]
    semantic <- root$semantic$representations[[name]]
    storage <- root$representations[[name]]
    if (!is.list(storage) || !identical(storage$type, semantic$type)) {
      stop("HDF5 study representation type mismatch: ", name, call. = FALSE)
    }
    if (identical(semantic$type, "fmri_frame")) {
      if (!identical(names(storage), c("type", "file", "semantic_digest"))) {
        stop("Invalid HDF5 study frame descriptor: ", name, call. = FALSE)
      }
      out[[name]] <- .open_study_frame(
        path,
        storage[c("file", "semantic_digest")],
        semantic$manifest,
        .study_representation_file(index),
        name
      )
      next
    }
    if (!identical(names(storage), c("type", "members")) ||
        !is.list(storage$members) ||
        !identical(names(storage$members), names(semantic$members))) {
      stop("Invalid HDF5 study collection descriptor: ", name, call. = FALSE)
    }
    members <- vector("list", length(semantic$members))
    names(members) <- names(semantic$members)
    for (member_index in seq_along(members)) {
      member_name <- names(members)[[member_index]]
      members[[member_name]] <- .open_study_frame(
        path,
        storage$members[[member_name]],
        semantic$members[[member_name]],
        .study_collection_member_file(index, member_index),
        paste0(name, ".", member_name)
      )
    }
    out[[name]] <- fmridataset::fmri_collection(
      members,
      metadata = semantic$metadata,
      provenance = semantic$provenance
    )
  }
  out
}

.open_study_bindings <- function(path, root) {
  semantic_arrays <- root$semantic$arrays
  storage <- root$shared_arrays
  if (!length(semantic_arrays)) {
    if (!identical(storage, list(file = NULL, arrays = list()))) {
      stop("HDF5 study declares unexpected shared arrays.", call. = FALSE)
    }
    return(list())
  }
  if (!is.list(storage) || !identical(names(storage), c("file", "arrays")) ||
      !identical(storage$file, "shared-arrays.fds.h5") ||
      !is.list(storage$arrays) ||
      !identical(names(storage$arrays), names(semantic_arrays))) {
    stop("Invalid HDF5 study shared-array manifest.", call. = FALSE)
  }
  array_path <- file.path(path, storage$file)
  if (!file.exists(array_path)) stop("Committed HDF5 study shared arrays are missing.", call. = FALSE)
  h5 <- hdf5r::H5File$new(array_path, mode = "r")
  on.exit(close_h5_safely(h5), add = TRUE)
  attribute <- function(name) {
    tryCatch(hdf5r::h5attr(h5, name), error = function(error) NULL)
  }
  if (!identical(attribute("fds_schema_id"), .study_h5_schema_id) ||
      !identical(attribute("fds_schema_version"), .study_h5_schema_version) ||
      !identical(attribute("fds_object_type"), "fmri_study_arrays") ||
      !identical(attribute("fds_commit_state"), "complete") ||
      !identical(attribute("fds_semantic_digest"), root$semantic_digest)) {
    stop("Unsupported, incomplete, or mismatched HDF5 study shared arrays.", call. = FALSE)
  }
  out <- lapply(names(semantic_arrays), function(key) {
    info <- storage$arrays[[key]]
    descriptor <- semantic_arrays[[key]]
    if (!is.list(info) || !identical(info$key, key) ||
        !identical(info$shape, descriptor$shape) ||
        !identical(info$dtype, descriptor$dtype)) {
      stop("HDF5 study shared-array metadata mismatch: ", key, call. = FALSE)
    }
    h5_array_source(array_path, info$dataset)
  })
  names(out) <- names(semantic_arrays)
  close_h5_safely(h5)
  out
}

#' @rdname write_study_h5
#' @export
open_study_h5 <- function(path) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  root <- .read_study_root(path)
  representations <- .open_study_representations(path, root)
  bindings <- .open_study_bindings(path, root)
  fmridataset::study_from_fds_manifest(root$semantic, representations, bindings)
}
