.frame_schema_id <- "org.fmridataset.fds/v1"
.frame_schema_version <- 1L

.validate_frame_assay_names <- function(x) {
  assay_names <- names(fmridataset::assays(x))
  valid <- grepl("^[A-Za-z][A-Za-z0-9_.-]*$", assay_names)
  if (!all(valid)) {
    stop(
      "HDF5 frame assay names must start with a letter and contain only ",
      "letters, digits, dot, underscore, or hyphen.",
      call. = FALSE
    )
  }
  invisible(assay_names)
}

.frame_h5_dtype <- function(dtype) {
  switch(
    dtype,
    logical = hdf5r::h5types$H5T_NATIVE_HBOOL,
    uint8 = hdf5r::h5types$H5T_STD_U8LE,
    int8 = hdf5r::h5types$H5T_STD_I8LE,
    uint16 = hdf5r::h5types$H5T_STD_U16LE,
    int16 = hdf5r::h5types$H5T_STD_I16LE,
    uint32 = hdf5r::h5types$H5T_STD_U32LE,
    int32 = hdf5r::h5types$H5T_STD_I32LE,
    uint64 = hdf5r::h5types$H5T_STD_U64LE,
    int64 = hdf5r::h5types$H5T_STD_I64LE,
    float16 = hdf5r::h5types$H5T_IEEE_F32LE,
    bfloat16 = hdf5r::h5types$H5T_IEEE_F32LE,
    float32 = hdf5r::h5types$H5T_IEEE_F32LE,
    float64 = hdf5r::h5types$H5T_IEEE_F64LE,
    stop("The HDF5 frame codec does not support dtype: ", dtype, call. = FALSE)
  )
}

.validate_frame_budget <- function(memory_budget) {
  if (!is.numeric(memory_budget) || length(memory_budget) != 1L ||
    is.na(memory_budget) || !is.finite(memory_budget) || memory_budget <= 0) {
    stop("`memory_budget` must be one positive finite number of bytes.", call. = FALSE)
  }
  as.double(memory_budget)
}

.frame_chunk_dims <- function(
    shape,
    dtype,
    memory_budget,
    requested = NULL,
    layout = c("balanced", "imagewise", "featurewise"),
    target_chunk_bytes = 1024^2) {
  layout <- match.arg(layout)
  shape <- as.integer(shape)
  if (length(shape) != 2L || anyNA(shape) || any(shape < 0L)) {
    stop("Frame array shape must contain two non-negative integers.", call. = FALSE)
  }
  chunk_shape <- pmax(shape, 1L)
  if (!is.numeric(target_chunk_bytes) || length(target_chunk_bytes) != 1L ||
    is.na(target_chunk_bytes) || !is.finite(target_chunk_bytes) ||
    target_chunk_bytes <= 0) {
    stop("`target_chunk_bytes` must be one positive finite number.", call. = FALSE)
  }
  bytes <- .h5_dtype_bytes(dtype)
  capacity <- floor(memory_budget / bytes)
  if (capacity < 1L) {
    stop("`memory_budget` cannot hold one source value.", call. = FALSE)
  }
  if (!is.null(requested)) {
    requested <- as.integer(requested)
    if (length(requested) != 2L || anyNA(requested) || any(requested < 1L)) {
      stop("`chunk_dims` must contain two positive integers.", call. = FALSE)
    }
    chunks <- pmin(requested, chunk_shape)
  } else {
    capacity <- min(capacity, floor(target_chunk_bytes / bytes))
    if (capacity < 1L) {
      stop("`target_chunk_bytes` cannot hold one source value.", call. = FALSE)
    }
    if (layout == "imagewise") {
      feature <- min(chunk_shape[[2L]], capacity)
      observation <- min(chunk_shape[[1L]], max(1L, floor(capacity / feature)))
      chunks <- as.integer(c(observation, feature))
    } else if (layout == "featurewise") {
      observation <- min(chunk_shape[[1L]], capacity)
      feature <- min(chunk_shape[[2L]], max(1L, floor(capacity / observation)))
      chunks <- as.integer(c(observation, feature))
    } else {
      base <- pmin(chunk_shape, c(64L, 4096L))
      scale <- sqrt(capacity / prod(as.double(base)))
      observation <- min(chunk_shape[[1L]], max(1L, floor(base[[1L]] * scale)))
      feature <- min(chunk_shape[[2L]], max(1L, floor(capacity / observation)))
      chunks <- as.integer(c(observation, feature))
    }
  }
  if (prod(as.double(chunks)) > capacity) {
    observation <- floor(sqrt(capacity * chunks[[1L]] / chunks[[2L]]))
    observation <- min(chunk_shape[[1L]], max(1L, observation))
    feature <- min(chunk_shape[[2L]], max(1L, floor(capacity / observation)))
    observation <- min(chunk_shape[[1L]], max(1L, floor(capacity / feature)))
    chunks <- as.integer(c(observation, feature))
  }
  chunks
}

.frame_rds_write <- function(h5, path, value) {
  raw <- as.integer(serialize(value, connection = NULL, version = 3L))
  h5$create_dataset(path, robj = raw)
}

.frame_rds_read <- function(h5, path) {
  if (!h5$exists(path)) stop("HDF5 frame manifest is missing: ", path, call. = FALSE)
  dataset <- h5[[path]]
  on.exit(close_h5_safely(dataset), add = TRUE)
  unserialize(as.raw(dataset$read()))
}

.frame_binding_read <- function(binding, observations, features) {
  if (inherits(binding, "array_source")) {
    fmridataset::source_read(
      binding,
      observations = observations,
      features = features
    )
  } else {
    binding[observations, features, drop = FALSE]
  }
}

.write_frame_array <- function(
    h5,
    binding,
    declaration,
    dataset_path,
    memory_budget,
    chunk_dims = NULL,
    layout = "balanced",
    target_chunk_bytes = 1024^2) {
  shape <- as.integer(declaration$shape)
  if (length(shape) != 2L) {
    stop(
      "The current HDF5 codec supports two-dimensional bindings only: ",
      declaration$key,
      call. = FALSE
    )
  }
  chunks <- .frame_chunk_dims(
    shape, declaration$dtype, memory_budget, chunk_dims,
    layout = layout, target_chunk_bytes = target_chunk_bytes
  )
  dataset <- h5$create_dataset(
    dataset_path,
    dtype = .frame_h5_dtype(declaration$dtype),
    dims = shape,
    chunk_dims = chunks
  )
  on.exit(close_h5_safely(dataset), add = TRUE)
  hdf5r::h5attr(dataset, "fds_dtype") <- declaration$dtype
  if (any(shape == 0L)) {
    return(list(
      key = declaration$key,
      dataset = paste0("/", dataset_path),
      shape = shape,
      dtype = declaration$dtype,
      chunks = chunks,
      layout = layout,
      target_chunk_bytes = as.double(target_chunk_bytes)
    ))
  }
  observation_starts <- seq.int(1L, shape[[1L]], by = chunks[[1L]])
  feature_starts <- seq.int(1L, shape[[2L]], by = chunks[[2L]])
  for (observation_start in observation_starts) {
    observations <- observation_start:min(
      shape[[1L]], observation_start + chunks[[1L]] - 1L
    )
    for (feature_start in feature_starts) {
      features <- feature_start:min(
        shape[[2L]], feature_start + chunks[[2L]] - 1L
      )
      values <- .frame_binding_read(binding, observations, features)
      dataset[observations, features] <- values
    }
  }
  list(
    key = declaration$key,
    dataset = paste0("/", dataset_path),
    shape = shape,
    dtype = declaration$dtype,
    chunks = chunks,
    layout = layout,
    target_chunk_bytes = as.double(target_chunk_bytes)
  )
}

#' Write and open FDS v1 HDF5 fmri frames
#'
#' `write_frame_h5()` validates the source-free FDS v1 semantic manifest, then
#' streams every separately bound array into a sibling temporary HDF5 file.
#' The destination becomes visible only after arrays, manifests, digests, and a
#' complete commit marker have been flushed. `open_frame_h5()` reads manifests
#' and dataset metadata only; numerical arrays remain lazy.
#'
#' @param x An `fmri_frame`.
#' @param path Destination or source HDF5 path.
#' @param overwrite Atomically replace an existing destination.
#' @param chunk_dims Optional two-element HDF5 chunk shape applied to arrays.
#' @param layout Workload profile used when `chunk_dims` is absent: balanced,
#'   imagewise, or featurewise.
#' @param target_chunk_bytes Desired physical chunk size. The hard
#'   `memory_budget` always takes precedence.
#' @param memory_budget Hard maximum bytes for one source-read block.
#' @return `write_frame_h5()` invisibly returns the normalized committed path;
#'   `open_frame_h5()` returns a lazy `fmri_frame`.
#' @export
write_frame_h5 <- function(
    x,
    path,
    overwrite = FALSE,
    chunk_dims = NULL,
    layout = c("balanced", "imagewise", "featurewise"),
    target_chunk_bytes = getOption("fmristore.target_chunk_bytes", 1024^2),
    memory_budget = getOption("fmridataset.block_budget", 512 * 1024^2)) {
  if (!inherits(x, "fmri_frame")) {
    stop("`x` must be an fmri_frame.", call. = FALSE)
  }
  if (!is.character(path) || length(path) != 1L || !nzchar(path)) {
    stop("`path` must be one non-empty file path.", call. = FALSE)
  }
  .validate_frame_assay_names(x)
  memory_budget <- .validate_frame_budget(memory_budget)
  layout <- match.arg(layout)
  semantic <- fmridataset::fds_frame_manifest(x)
  bindings <- fmridataset::fds_frame_bindings(x)
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  parent <- dirname(path)
  if (!dir.exists(parent)) {
    stop("Destination directory does not exist: ", parent, call. = FALSE)
  }
  if (file.exists(path) && !isTRUE(overwrite)) {
    stop("Destination already exists; set `overwrite = TRUE`.", call. = FALSE)
  }

  temporary <- tempfile(
    pattern = paste0(".", basename(path), ".partial-"),
    tmpdir = parent,
    fileext = ".h5"
  )
  committed <- FALSE
  on.exit({
    if (!committed && file.exists(temporary)) unlink(temporary)
  }, add = TRUE)

  h5 <- hdf5r::H5File$new(temporary, mode = "w")
  on.exit(close_h5_safely(h5), add = TRUE)
  hdf5r::h5attr(h5, "fds_schema_id") <- .frame_schema_id
  hdf5r::h5attr(h5, "fds_schema_version") <- .frame_schema_version
  hdf5r::h5attr(h5, "fds_object_type") <- "fmri_frame"
  hdf5r::h5attr(h5, "fds_commit_state") <- "writing"
  h5$create_group("arrays")
  h5$create_group("manifest")

  storage_arrays <- vector("list", length(semantic$arrays))
  names(storage_arrays) <- names(semantic$arrays)
  fail_after <- getOption(
    "fmristore.frame_writer_fault_after_array",
    getOption("fmristore.frame_writer_fault_after_assay", Inf)
  )
  for (index in seq_along(semantic$arrays)) {
    key <- names(semantic$arrays)[[index]]
    storage_arrays[[key]] <- .write_frame_array(
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
      stop("Injected HDF5 frame writer failure.", call. = FALSE)
    }
  }
  storage <- list(
    schema_id = .frame_schema_id,
    schema_version = .frame_schema_version,
    object_type = "fmri_frame",
    semantic_digest = fmridataset::fds_manifest_digest(semantic),
    arrays = storage_arrays
  )
  .frame_rds_write(h5, "manifest/semantic_rds", semantic)
  .frame_rds_write(h5, "manifest/storage_rds", storage)
  hdf5r::h5attr(h5, "fds_semantic_digest") <- storage$semantic_digest
  hdf5r::h5attr(h5, "fds_commit_state") <- "complete"
  h5$flush()
  close_h5_safely(h5)

  if (!file.rename(temporary, path)) {
    stop("Could not atomically commit HDF5 frame: ", path, call. = FALSE)
  }
  committed <- TRUE
  invisible(normalizePath(path, winslash = "/", mustWork = TRUE))
}

#' @rdname write_frame_h5
#' @export
open_frame_h5 <- function(path) {
  if (!is.character(path) || length(path) != 1L || !file.exists(path)) {
    stop("HDF5 frame path does not exist.", call. = FALSE)
  }
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  h5 <- hdf5r::H5File$new(path, mode = "r")
  on.exit(close_h5_safely(h5), add = TRUE)
  attribute <- function(name) {
    tryCatch(hdf5r::h5attr(h5, name), error = function(error) NULL)
  }
  if (!identical(attribute("fds_schema_id"), .frame_schema_id) ||
    !identical(attribute("fds_schema_version"), .frame_schema_version) ||
    !identical(attribute("fds_object_type"), "fmri_frame") ||
    !identical(attribute("fds_commit_state"), "complete")) {
    stop("Unsupported, incomplete, or missing FDS v1 frame schema.", call. = FALSE)
  }
  semantic <- .frame_rds_read(h5, "manifest/semantic_rds")
  storage <- .frame_rds_read(h5, "manifest/storage_rds")
  fmridataset::validate_fds_manifest(semantic)
  digest <- fmridataset::fds_manifest_digest(semantic)
  if (!identical(storage$semantic_digest, digest) ||
    !identical(attribute("fds_semantic_digest"), digest) ||
    !setequal(names(storage$arrays), names(semantic$arrays))) {
    stop("FDS frame storage manifest does not match its semantic manifest.", call. = FALSE)
  }
  bindings <- lapply(names(semantic$arrays), function(key) {
    info <- storage$arrays[[key]]
    if (!identical(info$key, key) || !identical(info$shape, semantic$arrays[[key]]$shape) ||
      !identical(info$dtype, semantic$arrays[[key]]$dtype)) {
      stop("FDS physical array metadata mismatch: ", key, call. = FALSE)
    }
    h5_array_source(path, info$dataset)
  })
  names(bindings) <- names(semantic$arrays)
  close_h5_safely(h5)
  fmridataset::frame_from_fds_manifest(semantic, bindings)
}
