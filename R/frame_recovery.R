.provisional_frame_schema_version <- "0.1-provisional"

.frame_h5_identity <- function(path) {
  tryCatch({
    h5 <- hdf5r::H5File$new(path, mode = "r")
    on.exit(close_h5_safely(h5), add = TRUE)
    attribute <- function(name) {
      tryCatch(hdf5r::h5attr(h5, name), error = function(error) NULL)
    }
    list(
      schema_id = attribute("fds_schema_id"),
      schema_version = attribute("fds_schema_version"),
      object_type = attribute("fds_object_type"),
      commit_state = attribute("fds_commit_state")
    )
  }, error = function(error) structure(list(message = conditionMessage(error)), class = "frame_h5_error"))
}

.frame_sibling_partials <- function(path) {
  parent <- dirname(path)
  if (!dir.exists(parent)) return(character())
  candidates <- list.files(parent, all.files = TRUE, full.names = TRUE, no.. = TRUE)
  prefix <- paste0(".", basename(path), ".partial-")
  candidates[startsWith(basename(candidates), prefix)]
}

.frame_recovery_plan <- function(path) {
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  targets <- .frame_sibling_partials(path)
  kinds <- rep("sibling_partial", length(targets))
  removable <- rep(TRUE, length(targets))

  if (dir.exists(path)) {
    root_partials <- list.files(
      path, pattern = "^\\.manifest-", all.files = TRUE, full.names = TRUE
    )
    targets <- c(targets, root_partials)
    kinds <- c(kinds, rep("manifest_partial", length(root_partials)))
    removable <- c(removable, rep(TRUE, length(root_partials)))

    shard_dir <- file.path(path, "shards")
    shard_files <- if (dir.exists(shard_dir)) {
      list.files(shard_dir, all.files = TRUE, full.names = TRUE, no.. = TRUE)
    } else {
      character()
    }
    shard_partials <- shard_files[startsWith(basename(shard_files), ".shard-")]
    targets <- c(targets, shard_partials)
    kinds <- c(kinds, rep("shard_partial", length(shard_partials)))
    removable <- c(removable, rep(TRUE, length(shard_partials)))

    manifest <- tryCatch(.read_frame_shard_root(path), error = function(error) NULL)
    if (!is.null(manifest)) {
      committed <- normalizePath(
        file.path(path, manifest$shards$.file), winslash = "/", mustWork = FALSE
      )
      stable <- shard_files[grepl("^s[0-9]{6}\\.fds\\.h5$", basename(shard_files))]
      orphans <- setdiff(normalizePath(stable, winslash = "/", mustWork = FALSE), committed)
      missing <- committed[!file.exists(committed)]
      targets <- c(targets, orphans, missing)
      kinds <- c(
        kinds,
        rep("orphan_shard", length(orphans)),
        rep("missing_committed_shard", length(missing))
      )
      removable <- c(removable, rep(TRUE, length(orphans)), rep(FALSE, length(missing)))
    }
  }

  data.frame(
    path = unname(targets),
    kind = unname(kinds),
    removable = unname(removable),
    stringsAsFactors = FALSE
  )
}

#' Inspect and recover FDS frame stores
#'
#' `frame_store_info()` classifies current FDS v1 files, the historical
#' `0.1-provisional` HDF5 layout, sharded frame directories, incomplete files,
#' and corrupt inputs without reading numerical arrays.
#'
#' `recover_frame_store()` reports or rolls back only artifacts whose names and
#' manifest status prove they are uncommitted: sibling temporary stores,
#' manifest/shard temporaries, and stable shard files absent from a valid root
#' manifest. Missing committed shards are reported but never removed.
#'
#' @param path Frame file or sharded-frame directory.
#' @param action Either `"report"` (read-only) or `"rollback"`.
#' @return `frame_store_info()` returns a serializable status list.
#'   `recover_frame_store()` returns a data frame of recovery actions.
#' @export
frame_store_info <- function(path) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("`path` must be one non-empty path.", call. = FALSE)
  }
  normalized <- normalizePath(path, winslash = "/", mustWork = FALSE)
  recovery <- .frame_recovery_plan(normalized)
  if (dir.exists(normalized)) {
    manifest <- tryCatch(.read_frame_shard_root(normalized), error = function(error) error)
    if (inherits(manifest, "error")) {
      return(list(
        path = normalized, kind = "sharded_frame", schema_id = .frame_shard_schema_id,
        schema_version = NA_integer_, status = "invalid", migratable = FALSE,
        message = conditionMessage(manifest), recovery = recovery
      ))
    }
    missing <- any(recovery$kind == "missing_committed_shard")
    return(list(
      path = normalized, kind = "sharded_frame", schema_id = manifest$schema_id,
      schema_version = manifest$schema_version,
      status = if (missing) "incomplete" else "complete", migratable = FALSE,
      message = NULL, recovery = recovery
    ))
  }
  if (!file.exists(normalized)) {
    return(list(
      path = normalized, kind = "missing", schema_id = NULL, schema_version = NULL,
      status = if (nrow(recovery)) "interrupted" else "missing", migratable = FALSE,
      message = NULL, recovery = recovery
    ))
  }
  identity <- .frame_h5_identity(normalized)
  if (inherits(identity, "frame_h5_error")) {
    return(list(
      path = normalized, kind = "hdf5", schema_id = NULL, schema_version = NULL,
      status = "corrupt", migratable = FALSE, message = identity$message,
      recovery = recovery
    ))
  }
  current <- identical(identity$schema_id, .frame_schema_id) &&
    identical(identity$schema_version, .frame_schema_version) &&
    identical(identity$object_type, "fmri_frame")
  provisional <- identical(identity$schema_version, .provisional_frame_schema_version) &&
    identical(identity$object_type, "fmri_frame")
  status <- if (current && identical(identity$commit_state, "complete")) {
    "complete"
  } else if (current) {
    "incomplete"
  } else if (provisional) {
    "legacy"
  } else {
    "unsupported"
  }
  list(
    path = normalized,
    kind = "hdf5_frame",
    schema_id = identity$schema_id,
    schema_version = identity$schema_version,
    status = status,
    migratable = provisional,
    message = NULL,
    recovery = recovery
  )
}

#' @rdname frame_store_info
#' @export
recover_frame_store <- function(path, action = c("report", "rollback")) {
  action <- match.arg(action)
  plan <- .frame_recovery_plan(path)
  plan$removed <- FALSE
  if (action == "rollback" && nrow(plan)) {
    for (index in which(plan$removable)) {
      target <- plan$path[[index]]
      existed <- file.exists(target) || dir.exists(target)
      if (existed) unlink(target, recursive = dir.exists(target))
      plan$removed[[index]] <- existed && !file.exists(target) && !dir.exists(target)
    }
  }
  plan
}

.read_provisional_frame_manifest <- function(path) {
  h5 <- hdf5r::H5File$new(path, mode = "r")
  on.exit(close_h5_safely(h5), add = TRUE)
  if (!h5$exists("manifest/rds")) stop("Provisional frame manifest is missing.", call. = FALSE)
  dataset <- h5[["manifest/rds"]]
  on.exit(close_h5_safely(dataset), add = TRUE)
  manifest <- unserialize(as.raw(dataset$read()))
  required <- c(
    "schema_version", "observations", "features", "entities", "relations",
    "tables", "active_assay", "metadata", "provenance", "assays"
  )
  if (!is.list(manifest) || !all(required %in% names(manifest)) ||
    !identical(manifest$schema_version, .provisional_frame_schema_version) ||
    !is.list(manifest$assays) || !length(manifest$assays)) {
    stop("Invalid provisional frame manifest.", call. = FALSE)
  }
  manifest
}

.open_provisional_frame_h5 <- function(path) {
  info <- frame_store_info(path)
  if (!isTRUE(info$migratable)) stop("HDF5 frame is not a supported legacy schema.", call. = FALSE)
  manifest <- .read_provisional_frame_manifest(info$path)
  assay_values <- lapply(manifest$assays, function(assay) {
    if (!is.list(assay) || !identical(
      assay$dataset, paste0("/assays/", assay$name)
    )) {
      stop("Invalid provisional assay declaration.", call. = FALSE)
    }
    source <- h5_array_source(info$path, assay$dataset)
    structure(
      list(
        name = assay$name,
        source = source,
        dtype = fmridataset::source_dtype(source),
        observation_digest = NULL,
        feature_digest = NULL,
        role = assay$role,
        units = assay$units,
        metadata = assay$metadata
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

#' Migrate historical HDF5 frames to FDS v1
#'
#' The migration registry is explicit and currently contains the real
#' `0.1-provisional` to FDS v1 transition. Migration reopens legacy assays as
#' lazy HDF5 sources and streams them through the atomic current writer.
#'
#' @param path Historical HDF5 frame path.
#' @param destination Output path. May equal `path` only with
#'   `overwrite = TRUE`; failed in-place migrations preserve the source file.
#' @param overwrite Replace an existing destination atomically.
#' @param ... Additional arguments passed to [write_frame_h5()].
#' @return `frame_schema_migrations()` returns the supported migration table;
#'   `migrate_frame_h5()` invisibly returns the committed destination.
#' @export
frame_schema_migrations <- function() {
  data.frame(
    from_schema = .provisional_frame_schema_version,
    to_schema = .frame_schema_id,
    to_version = .frame_schema_version,
    stringsAsFactors = FALSE
  )
}

#' @rdname frame_schema_migrations
#' @export
migrate_frame_h5 <- function(path, destination, overwrite = FALSE, ...) {
  if (missing(destination) || !is.character(destination) || length(destination) != 1L ||
    is.na(destination) || !nzchar(destination)) {
    stop("`destination` must be one non-empty path.", call. = FALSE)
  }
  source_path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  destination <- normalizePath(destination, winslash = "/", mustWork = FALSE)
  if (identical(source_path, destination) && !isTRUE(overwrite)) {
    stop("In-place migration requires `overwrite = TRUE`.", call. = FALSE)
  }
  frame <- .open_provisional_frame_h5(source_path)
  write_frame_h5(frame, destination, overwrite = overwrite, ...)
}
