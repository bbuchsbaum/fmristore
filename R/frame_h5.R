.frame_schema_version <- "0.1-provisional"

.frame_manifest <- function(x) {
  assay_manifest <- lapply(fmridataset::assays(x), function(value) {
    list(
      name = value$name,
      dataset = paste0("/assays/", value$name),
      dtype = value$dtype,
      role = value$role,
      units = value$units,
      metadata = value$metadata
    )
  })
  list(
    schema_version = .frame_schema_version,
    observations = x$observations,
    features = x$features,
    entities = x$entities,
    relations = x$relations,
    tables = x$tables,
    active_assay = fmridataset::active_assay(x),
    metadata = x$metadata,
    provenance = x$provenance,
    assays = assay_manifest
  )
}

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

.frame_chunk_dims <- function(shape, requested = NULL) {
  if (!is.null(requested)) {
    requested <- as.integer(requested)
    if (length(requested) != 2L || anyNA(requested) || any(requested < 1L)) {
      stop("`chunk_dims` must contain two positive integers.", call. = FALSE)
    }
    return(pmin(requested, shape))
  }
  pmin(as.integer(shape), c(64L, 4096L))
}

#' Write and open provisional HDF5 fmri frames
#'
#' `write_frame_h5()` persists aligned assays plus a serializable semantic
#' manifest. `open_frame_h5()` reads only the manifest and HDF5 dataset
#' metadata; assay values remain lazy through [h5_array_source()].
#'
#' The provisional schema is intentionally versioned separately from the
#' frozen FDS 1.0 schema. Files are written through a sibling temporary file
#' and renamed only after every assay and the manifest have been committed.
#'
#' @param x An `fmri_frame`.
#' @param path Destination or source HDF5 path.
#' @param overwrite Replace an existing destination.
#' @param chunk_dims Optional two-element HDF5 chunk shape.
#' @param memory_budget Maximum bytes allowed for one realized assay.
#' @return `write_frame_h5()` invisibly returns the normalized path;
#'   `open_frame_h5()` returns an `fmri_frame` backed by HDF5 sources.
#' @export
write_frame_h5 <- function(
  x,
  path,
  overwrite = FALSE,
  chunk_dims = NULL,
  memory_budget = getOption("fmridataset.collect_budget", 2 * 1024^3)
) {
  if (!inherits(x, "fmri_frame")) {
    stop("`x` must be an fmri_frame.", call. = FALSE)
  }
  if (!is.character(path) || length(path) != 1L || !nzchar(path)) {
    stop("`path` must be one non-empty file path.", call. = FALSE)
  }
  .validate_frame_assay_names(x)
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
  on.exit(
    {
      if (!committed && file.exists(temporary)) unlink(temporary)
    },
    add = TRUE
  )

  h5 <- hdf5r::H5File$new(temporary, mode = "w")
  on.exit(close_h5_safely(h5), add = TRUE)
  hdf5r::h5attr(h5, "fds_schema_version") <- .frame_schema_version
  hdf5r::h5attr(h5, "fds_object_type") <- "fmri_frame"
  h5$create_group("assays")
  h5$create_group("manifest")

  assay_names <- names(fmridataset::assays(x))
  fail_after <- getOption("fmristore.frame_writer_fault_after_assay", Inf)
  for (index in seq_along(assay_names)) {
    assay_name <- assay_names[[index]]
    value <- fmridataset::collect_assay(
      x,
      assay = assay_name,
      memory_budget = memory_budget
    )
    chunks <- .frame_chunk_dims(dim(value), chunk_dims)
    h5$create_dataset(
      paste0("assays/", assay_name),
      robj = value,
      chunk_dims = chunks
    )
    if (is.finite(fail_after) && index >= fail_after) {
      stop("Injected HDF5 frame writer failure.", call. = FALSE)
    }
  }

  manifest <- as.integer(serialize(
    .frame_manifest(x),
    connection = NULL,
    version = 3L
  ))
  h5$create_dataset("manifest/rds", robj = manifest)
  h5$flush()
  close_h5_safely(h5)

  if (file.exists(path) && unlink(path) != 0L) {
    stop("Could not replace existing destination: ", path, call. = FALSE)
  }
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
  schema <- tryCatch(
    hdf5r::h5attr(h5, "fds_schema_version"),
    error = function(e) NULL
  )
  object_type <- tryCatch(
    hdf5r::h5attr(h5, "fds_object_type"),
    error = function(e) NULL
  )
  if (!identical(schema, .frame_schema_version) ||
    !identical(object_type, "fmri_frame")) {
    stop("Unsupported or missing provisional FDS frame schema.", call. = FALSE)
  }
  if (!h5$exists("manifest/rds")) {
    stop("HDF5 frame manifest is missing.", call. = FALSE)
  }
  dataset <- h5[["manifest/rds"]]
  on.exit(close_h5_safely(dataset), add = TRUE)
  manifest <- unserialize(as.raw(dataset$read()))
  close_h5_safely(dataset)
  close_h5_safely(h5)

  assay_values <- lapply(manifest$assays, function(info) {
    source <- h5_array_source(path, info$dataset)
    structure(
      list(
        name = info$name,
        source = source,
        dtype = fmridataset::source_dtype(source),
        observation_digest = NULL,
        feature_digest = NULL,
        role = info$role,
        units = info$units,
        metadata = info$metadata
      ),
      class = "aligned_assay"
    )
  })
  names(assay_values) <- vapply(manifest$assays, `[[`, character(1), "name")

  fmridataset::fmri_frame(
    assays = assay_values,
    observations = manifest$observations,
    features = manifest$features,
    entities = manifest$entities,
    relations = manifest$relations,
    tables = manifest$tables,
    active_assay = manifest$active_assay,
    metadata = manifest$metadata,
    provenance = manifest$provenance
  )
}
