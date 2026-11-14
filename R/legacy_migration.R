# Convert frozen legacy layouts (parcellated experiments, labeled volume sets)
# into FDS v1 frames, so their data outlive the legacy readers.

.legacy_volume_space <- function(mask, template = NULL) {
  support <- which(as.logical(mask@.Data))
  if (!length(support)) stop("The legacy file's mask contains no voxels.", call. = FALSE)
  fmridataset::volume_space(
    dim = as.integer(dim(mask)[1:3]),
    affine = unname(neuroim2::trans(neuroim2::space(mask))),
    support = support,
    template = template
  )
}

# Scalar atomic columns only: frame tables and feature data must be serializable
# scalars, and legacy metadata can hold anything.
.scalar_columns <- function(df) {
  if (!is.data.frame(df) || !ncol(df)) return(NULL)
  keep <- vapply(df, function(col) is.atomic(col) && is.null(dim(col)), logical(1))
  df[, keep, drop = FALSE]
}

.parcellated_frame <- function(path, kind = NULL, atlas = "fmristore-clusters", template = NULL) {
  experiment <- suppressMessages(H5ParcellatedMultiScan(path))
  on.exit(close(experiment), add = TRUE)
  runs <- experiment@runs
  if (!length(runs)) stop("The parcellated file contains no scans.", call. = FALSE)

  is_summary <- vapply(runs, methods::is, logical(1), "H5ParcellatedScanSummary")
  present <- unique(ifelse(is_summary, "summary", "full"))
  if (is.null(kind)) {
    if (length(present) > 1L) {
      stop(
        "The file mixes full and summary scans, which cannot share one feature axis. ",
        "Convert each with `kind = \"full\"` and `kind = \"summary\"`.",
        call. = FALSE
      )
    }
    kind <- present
  }
  kind <- match.arg(kind, c("summary", "full"))
  runs <- runs[is_summary == (kind == "summary")]
  if (!length(runs)) stop("The file contains no ", kind, " scans.", call. = FALSE)

  mask <- runs[[1L]]@mask
  cluster <- as.integer(clusters(experiment)@clusters)
  ids <- sort(unique(cluster))
  voxels <- .legacy_volume_space(mask, template)

  values <- lapply(runs, function(run) {
    if (kind == "summary") {
      matrix_values <- as.matrix(run)
      matrix_values[, match(ids, run@cluster_ids), drop = FALSE]
    } else {
      as.matrix(series(run, seq_along(cluster)))
    }
  })
  n_time <- vapply(values, nrow, integer(1))
  run_names <- names(runs)
  values <- do.call(rbind, values)
  dimnames(values) <- NULL

  observations <- fmridataset::axis_frame(data.frame(
    .obs_id = paste0(rep(run_names, n_time), ":", unlist(lapply(n_time, seq_len))),
    run_id = rep(run_names, n_time),
    t = unlist(lapply(n_time, seq_len)),
    stringsAsFactors = FALSE
  ))

  cluster_meta <- .scalar_columns(experiment@cluster_metadata)
  if (kind == "summary") {
    parcel_data <- data.frame(
      id = as.character(ids),
      n_voxels = as.integer(tabulate(match(cluster, ids), nbins = length(ids))),
      stringsAsFactors = FALSE
    )
    if (!is.null(cluster_meta) && "cluster_id" %in% names(cluster_meta)) {
      extra <- cluster_meta[match(ids, cluster_meta$cluster_id), setdiff(names(cluster_meta), c("cluster_id", "n_voxels")), drop = FALSE]
      parcel_data <- cbind(parcel_data, extra)
    }
    parcels <- fmridataset::parcel_space(
      voxels,
      parcel_ids = ids,
      membership = Matrix::sparseMatrix(
        i = seq_along(cluster), j = match(cluster, ids), x = 1,
        dims = c(length(cluster), length(ids))
      ),
      data = parcel_data,
      atlas = atlas
    )
    features <- fmridataset::feature_axis(
      data.frame(.feature_id = fmridataset::feature_ids(parcels)), space = parcels
    )
    assay_name <- "summary"
  } else {
    features <- fmridataset::feature_axis(
      data.frame(.feature_id = fmridataset::feature_ids(voxels), cluster_id = cluster),
      space = voxels
    )
    assay_name <- "signal"
  }

  scans <- data.frame(run_id = run_names, n_time = n_time, stringsAsFactors = FALSE)
  assays <- list(values)
  names(assays) <- assay_name
  fmridataset::fmri_frame(
    assays = assays,
    observations = observations,
    features = features,
    tables = list(scans = fmridataset::auxiliary_table(scans, key = "run_id")),
    active_assay = assay_name,
    metadata = list(
      converted_from = "H5ParcellatedMultiScan",
      source_file = basename(path),
      kind = kind
    )
  )
}

.labeled_frame <- function(path, template = NULL) {
  lvs <- read_labeled_vec(path)
  on.exit(close(lvs), add = TRUE)
  labels <- lvs@labels
  if (anyDuplicated(labels)) stop("The labeled file has duplicate labels.", call. = FALSE)
  voxels <- .legacy_volume_space(lvs@mask, template)
  support <- which(as.logical(lvs@mask@.Data))
  values <- t(vapply(labels, function(label) as.numeric(lvs[[label]]@.Data[support]), numeric(length(support))))
  dimnames(values) <- NULL

  fmridataset::fmri_frame(
    assays = list(value = values),
    observations = fmridataset::axis_frame(data.frame(.obs_id = labels, stringsAsFactors = FALSE)),
    features = fmridataset::feature_axis(
      data.frame(.feature_id = fmridataset::feature_ids(voxels)), space = voxels
    ),
    active_assay = "value",
    metadata = list(converted_from = "LabeledVolumeSet", source_file = basename(path))
  )
}

#' Convert legacy parcellated and labeled files to FDS v1 frames
#'
#' The parcellated and labeled classes are frozen and will be removed. These
#' functions read their files and write the data as FDS v1 frames with
#' [write_frame_h5()], so the files stay readable once the legacy readers are
#' gone.
#'
#' `migrate_parcellated_h5()` converts a file written by
#' [write_parcellated_experiment_h5()]. Observations are the time points of
#' every scan, in the order the scans were written, with `run_id` and `t`
#' columns; the frame's `scans` table gives each run's length.
#' - **Summary scans** become an observation x parcel assay `"summary"` over a
#'   `fmridataset::parcel_space()` whose membership records each parcel's
#'   voxels within the mask; cluster metadata become parcel columns.
#' - **Full scans** become an observation x voxel assay `"signal"` over the
#'   mask's `fmridataset::volume_space()`, with a `cluster_id` feature column.
#'
#' A frame has one feature axis, so a file that mixes full and summary scans
#' is converted one kind at a time with `kind`.
#'
#' `migrate_labeled_h5()` converts a file written by [write_labeled_vec()].
#' Observations are the labels and features are the mask voxels, in assay
#' `"value"`.
#'
#' Both read the data into memory before writing.
#'
#' @param path A legacy parcellated or labeled HDF5 file.
#' @param destination Output path for the frame.
#' @param kind For parcellated files: `"summary"` or `"full"`. Needed only when
#'   the file contains both.
#' @param atlas Atlas identity recorded on the `parcel_space` of a summary frame.
#' @param template Optional template identity recorded on the `volume_space`.
#' @param overwrite Replace an existing destination atomically.
#' @param ... Further arguments passed to [write_frame_h5()].
#' @return The committed destination, invisibly.
#' @export
migrate_parcellated_h5 <- function(path, destination, kind = NULL,
                                   atlas = "fmristore-clusters", template = NULL,
                                   overwrite = FALSE, ...) {
  .check_legacy_paths(path, destination, overwrite)
  frame <- .parcellated_frame(path, kind = kind, atlas = atlas, template = template)
  write_frame_h5(frame, destination, overwrite = overwrite, ...)
}

#' @rdname migrate_parcellated_h5
#' @export
migrate_labeled_h5 <- function(path, destination, template = NULL, overwrite = FALSE, ...) {
  .check_legacy_paths(path, destination, overwrite)
  frame <- .labeled_frame(path, template = template)
  write_frame_h5(frame, destination, overwrite = overwrite, ...)
}

.check_legacy_paths <- function(path, destination, overwrite) {
  if (!is.character(path) || length(path) != 1L || !file.exists(path)) {
    stop("`path` must be one existing file.", call. = FALSE)
  }
  if (missing(destination) || !is.character(destination) || length(destination) != 1L ||
    is.na(destination) || !nzchar(destination)) {
    stop("`destination` must be one non-empty path.", call. = FALSE)
  }
  if (file.exists(destination) && !isTRUE(overwrite)) {
    stop("Destination already exists; set `overwrite = TRUE`.", call. = FALSE)
  }
  invisible(TRUE)
}
