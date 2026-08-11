#' HDF5-backed fmridataset array source
#'
#' Creates a reconstructible descriptor for one two-dimensional HDF5 dataset.
#' The descriptor contains no live HDF5 handles; handles are opened only for a
#' read or an explicit source session.
#'
#' @param path Path to an HDF5 file.
#' @param dataset Absolute or file-relative dataset path.
#' @param physical_axes Order of the HDF5 dimensions. The logical source is
#'   always observation by feature.
#' @return A serializable `h5_array_source` implementing the
#'   [fmridataset::source_shape()] source protocol.
#' @export
#' @importFrom fmridataset source_capabilities source_chunks source_close
#'   source_dtype source_fingerprint source_open source_read source_shape
h5_array_source <- function(
  path,
  dataset,
  physical_axes = c("observation", "feature")
) {
  if (!is.character(path) || length(path) != 1L || !nzchar(path)) {
    stop("`path` must be one non-empty HDF5 file path.", call. = FALSE)
  }
  if (!file.exists(path)) {
    stop("HDF5 source file does not exist: ", path, call. = FALSE)
  }
  if (!is.character(dataset) || length(dataset) != 1L || !nzchar(dataset)) {
    stop("`dataset` must be one non-empty HDF5 dataset path.", call. = FALSE)
  }
  if (!setequal(physical_axes, c("observation", "feature")) ||
    length(physical_axes) != 2L) {
    stop(
      "`physical_axes` must contain observation and feature exactly once.",
      call. = FALSE
    )
  }

  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  info <- .h5_source_metadata(path, dataset)
  logical_order <- match(c("observation", "feature"), physical_axes)

  structure(
    list(
      uri = path,
      dataset = dataset,
      shape = as.integer(info$shape[logical_order]),
      dtype = info$dtype,
      chunks = as.integer(info$chunks[logical_order]),
      physical_axes = physical_axes,
      fingerprint = .h5_source_fingerprint(path, dataset, info),
      schema_version = 1L
    ),
    class = c("h5_array_source", "array_source")
  )
}

.h5_source_metadata <- function(path, dataset) {
  h5 <- hdf5r::H5File$new(path, mode = "r")
  on.exit(h5$close_all(), add = TRUE)
  if (!h5$exists(dataset)) {
    stop("HDF5 dataset does not exist: ", dataset, call. = FALSE)
  }
  dset <- h5[[dataset]]
  on.exit(close_h5_safely(dset), add = TRUE)
  shape <- as.integer(dset$dims)
  if (length(shape) != 2L) {
    stop("An HDF5 ArraySource dataset must be two dimensional.", call. = FALSE)
  }
  declared_dtype <- tryCatch(
    hdf5r::h5attr(dset, "fds_dtype"),
    error = function(error) NULL
  )
  list(
    shape = shape,
    dtype = declared_dtype %||% .h5_source_dtype(dset),
    chunks = .h5_source_chunks(dset, shape)
  )
}

.h5_source_dtype <- function(dset) {
  type <- dset$get_type()
  on.exit(close_h5_safely(type), add = TRUE)
  size <- as.integer(type$get_size())
  class <- type$get_class()
  if (class == hdf5r::h5const$H5T_FLOAT) {
    return(if (size <= 4L) "float32" else "float64")
  }
  if (class == hdf5r::h5const$H5T_INTEGER) {
    signed <- type$get_sign() == hdf5r::h5const$H5T_SGN_2
    return(paste0(if (signed) "int" else "uint", 8L * size))
  }
  paste0("hdf5_", size * 8L)
}

.h5_source_chunks <- function(dset, shape) {
  plist <- tryCatch(dset$get_create_plist(), error = function(e) NULL)
  if (is.null(plist)) {
    return(shape)
  }
  on.exit(close_h5_safely(plist), add = TRUE)
  chunks <- tryCatch(
    plist$get_chunk(max_ndims = length(shape)),
    error = function(e) NULL
  )
  if (is.null(chunks) || length(chunks) != 2L) shape else as.integer(chunks)
}

.h5_source_fingerprint <- function(path, dataset, metadata) {
  info <- file.info(path)
  paste(
    "h5",
    path,
    dataset,
    info$size,
    format(info$mtime, tz = "UTC", usetz = TRUE),
    paste(metadata$shape, collapse = "x"),
    metadata$dtype,
    sep = ":"
  )
}

.h5_source_index <- function(index, n, axis) {
  if (is.null(index)) {
    return(seq_len(n))
  }
  if (is.logical(index)) {
    if (length(index) != n || anyNA(index)) {
      stop("Logical ", axis, " selectors must match the axis and contain no NA.", call. = FALSE)
    }
    return(which(index))
  }
  index <- as.integer(index)
  if (anyNA(index) || any(index < 1L | index > n)) {
    stop(axis, " selector is out of bounds.", call. = FALSE)
  }
  index
}

.h5_dtype_bytes <- function(dtype) {
  sizes <- c(
    logical = 1, uint8 = 1, int8 = 1,
    uint16 = 2, int16 = 2, float16 = 2, bfloat16 = 2,
    uint32 = 4, int32 = 4, float32 = 4,
    uint64 = 8, int64 = 8, float64 = 8
  )
  value <- unname(sizes[dtype])
  if (length(value) != 1L || is.na(value)) {
    stop("Unsupported HDF5 source dtype: ", dtype, call. = FALSE)
  }
  value
}

#' Plan and enforce HDF5 read amplification
#'
#' `h5_read_plan()` estimates the physical HDF5 chunks touched by a logical
#' observation-by-feature selection. It uses descriptor metadata only and does
#' not open or read the dataset. `validate_h5_read_amplification()` turns the
#' estimate into an executable performance gate.
#'
#' @param x An `h5_array_source` descriptor.
#' @param observations,features Logical selectors accepted by
#'   [fmridataset::source_read()]. Duplicate positions do not cause additional
#'   physical reads.
#' @param plan A plan returned by `h5_read_plan()`.
#' @param max Maximum permitted ratio of physically touched values to unique
#'   requested values.
#' @return `h5_read_plan()` returns a serializable `h5_read_plan` list;
#'   `validate_h5_read_amplification()` returns `plan` invisibly or errors.
#' @export
h5_read_plan <- function(x, observations = NULL, features = NULL) {
  if (!inherits(x, "h5_array_source")) {
    stop("`x` must be an h5_array_source descriptor.", call. = FALSE)
  }
  observations <- unique(.h5_source_index(observations, x$shape[[1L]], "Observation"))
  features <- unique(.h5_source_index(features, x$shape[[2L]], "Feature"))
  bytes_per_value <- .h5_dtype_bytes(x$dtype)
  requested_values <- length(observations) * length(features)

  touched_extent <- function(index, shape, chunk) {
    if (!length(index)) return(0)
    chunk_id <- unique((index - 1L) %/% chunk)
    starts <- chunk_id * chunk + 1
    sum(pmin(chunk, shape - starts + 1))
  }
  physical_values <- touched_extent(observations, x$shape[[1L]], x$chunks[[1L]]) *
    touched_extent(features, x$shape[[2L]], x$chunks[[2L]])
  amplification <- if (requested_values == 0) 0 else physical_values / requested_values

  structure(
    list(
      shape = x$shape,
      chunks = x$chunks,
      dtype = x$dtype,
      requested_values = as.double(requested_values),
      physical_values = as.double(physical_values),
      requested_bytes = as.double(requested_values * bytes_per_value),
      physical_bytes = as.double(physical_values * bytes_per_value),
      amplification = as.double(amplification)
    ),
    class = "h5_read_plan"
  )
}

#' @rdname h5_read_plan
#' @export
validate_h5_read_amplification <- function(plan, max = 1.5) {
  if (!inherits(plan, "h5_read_plan")) {
    stop("`plan` must be an h5_read_plan.", call. = FALSE)
  }
  if (!is.numeric(max) || length(max) != 1L || is.na(max) ||
    !is.finite(max) || max < 1) {
    stop("`max` must be one finite number greater than or equal to one.", call. = FALSE)
  }
  if (plan$amplification > max) {
    stop(
      "HDF5 read amplification ", format(plan$amplification, digits = 4L),
      " exceeds the configured gate of ", format(max, digits = 4L), ".",
      call. = FALSE
    )
  }
  invisible(plan)
}

#' @export
source_shape.h5_array_source <- function(x, ...) x$shape

#' @export
source_dtype.h5_array_source <- function(x, ...) x$dtype

#' @export
source_chunks.h5_array_source <- function(x, ...) x$chunks

#' @export
source_capabilities.h5_array_source <- function(x, ...) {
  c("row_slice", "column_slice", "block_slice", "serializable")
}

#' @export
source_fingerprint.h5_array_source <- function(x, ...) x$fingerprint

#' @export
source_open.h5_array_source <- function(x, ...) {
  h5 <- hdf5r::H5File$new(x$uri, mode = "r")
  if (!h5$exists(x$dataset)) {
    h5$close_all()
    stop("HDF5 dataset does not exist: ", x$dataset, call. = FALSE)
  }
  dset <- h5[[x$dataset]]
  structure(
    list(source = x, file = h5, dataset = dset),
    class = c("h5_array_source_handle", "array_source_handle")
  )
}

#' @export
source_read.h5_array_source <- function(
  x,
  observations = NULL,
  features = NULL,
  ...
) {
  handle <- source_open(x)
  on.exit(source_close(handle), add = TRUE)
  source_read(handle, observations = observations, features = features, ...)
}

#' @export
source_close.h5_array_source <- function(x, ...) invisible(TRUE)

#' @export
source_shape.h5_array_source_handle <- function(x, ...) x$source$shape

#' @export
source_dtype.h5_array_source_handle <- function(x, ...) x$source$dtype

#' @export
source_chunks.h5_array_source_handle <- function(x, ...) x$source$chunks

#' @export
source_capabilities.h5_array_source_handle <- function(x, ...) {
  source_capabilities(x$source)
}

#' @export
source_fingerprint.h5_array_source_handle <- function(x, ...) {
  source_fingerprint(x$source)
}

#' @export
source_read.h5_array_source_handle <- function(
  x,
  observations = NULL,
  features = NULL,
  ...
) {
  if (!isTRUE(x$file$is_valid) || !isTRUE(x$dataset$is_valid)) {
    stop("HDF5 ArraySource handle is closed.", call. = FALSE)
  }
  shape <- x$source$shape
  observations <- .h5_source_index(observations, shape[1L], "Observation")
  features <- .h5_source_index(features, shape[2L], "Feature")

  logical_index <- list(observations, features)
  physical_index <- logical_index[match(
    x$source$physical_axes,
    c("observation", "feature")
  )]
  unique_index <- lapply(physical_index, function(i) sort(unique(i)))
  value <- x$dataset[
    unique_index[[1L]],
    unique_index[[2L]],
    drop = FALSE
  ]
  value <- value[
    match(physical_index[[1L]], unique_index[[1L]]),
    match(physical_index[[2L]], unique_index[[2L]]),
    drop = FALSE
  ]
  if (identical(x$source$physical_axes, c("feature", "observation"))) {
    value <- t(value)
  }
  value
}

#' @export
source_close.h5_array_source_handle <- function(x, ...) {
  close_h5_safely(x$dataset)
  close_h5_safely(x$file)
  invisible(TRUE)
}
