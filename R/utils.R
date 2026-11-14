# Small helpers shared by the frame store.

#' Null default
#'
#' Provided for compatibility with R < 4.4.0, which lacks the base version.
#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Safely close an HDF5 object
#'
#' Checks if the object is valid and not NULL before attempting to close.
#' Only the given handle is closed. For \code{H5File} objects this
#' deliberately avoids \code{$close_all()}: HDF5 shares one file structure
#' between every handle opened on the same path, so \code{close_all()} also
#' closes handles owned by unrelated objects, such as an open frame source.
#' Unreferenced child handles leaked by
#' the caller are finalized by a garbage collection before the file handle is
#' closed, so they do not keep the file open.
#'
#' @param h5obj An HDF5 object (H5File, H5Group, H5D, etc.) or NULL.
#' @keywords internal
close_h5_safely <- function(h5obj) {
  if (is.null(h5obj) || !isTRUE(tryCatch(h5obj$is_valid, error = function(e) FALSE))) {
    return(invisible(NULL))
  }
  if (inherits(h5obj, "H5File")) {
    try(h5obj$flush(), silent = TRUE)
    gc(verbose = FALSE)
  }
  try(h5obj$close(), silent = TRUE)
  invisible(NULL)
}
