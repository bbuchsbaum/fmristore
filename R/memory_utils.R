#' Configure Memory Usage Warnings
#'
#' @description
#' Control memory usage warnings for large dataset operations in fmristore.
#' These warnings help users understand when operations may consume significant memory.
#'
#' @param enable Logical. Whether to enable memory usage warnings (default: TRUE).
#' @param threshold_mb Numeric. Memory threshold in MB above which warnings are issued (default: 100).
#' @param verbose Logical. Whether to print confirmation of setting changes (default: FALSE).
#'
#' @details
#' The memory warning system estimates the size of data that will be loaded into
#' memory and issues warnings when this exceeds the specified threshold. This helps
#' users make informed decisions about subsetting large datasets.
#'
#' Settings are stored as R options and persist for the current R session:
#' - `fmristore.warn_memory`: Enable/disable warnings
#' - `fmristore.memory_threshold_mb`: Warning threshold in MB
#'
#' @examples
#' # Enable warnings for datasets larger than 50MB
#' configure_memory_warnings(enable = TRUE, threshold_mb = 50)
#'
#' # Disable memory warnings
#' configure_memory_warnings(enable = FALSE)
#'
#' # Check current settings
#' list(
#'   enabled = getOption("fmristore.warn_memory", TRUE),
#'   threshold_mb = getOption("fmristore.memory_threshold_mb", 100)
#' )
#'
#' @export
configure_memory_warnings <- function(enable = TRUE, threshold_mb = 100, verbose = FALSE) {
  if (!is.logical(enable) || length(enable) != 1) {
    stop("'enable' must be a single logical value.")
  }
  if (!is.numeric(threshold_mb) || length(threshold_mb) != 1 || threshold_mb <= 0) {
    stop("'threshold_mb' must be a single positive number.")
  }
  if (!is.logical(verbose) || length(verbose) != 1) {
    stop("'verbose' must be a single logical value.")
  }
  
  options(
    fmristore.warn_memory = enable,
    fmristore.memory_threshold_mb = threshold_mb
  )
  
  if (verbose) {
    message(sprintf(
      "[fmristore] Memory warnings %s. Threshold: %.1f MB",
      if (enable) "enabled" else "disabled",
      threshold_mb
    ))
  }
  
  invisible(list(
    enabled = enable,
    threshold_mb = threshold_mb
  ))
}

#' Get Current Memory Warning Settings
#'
#' @description
#' Retrieve the current memory warning configuration.
#'
#' @return A list with elements `enabled` and `threshold_mb`.
#' @export
#'
#' @examples
#' # Get current settings
#' get_memory_settings()
get_memory_settings <- function() {
  list(
    enabled = getOption("fmristore.warn_memory", TRUE),
    threshold_mb = getOption("fmristore.memory_threshold_mb", 100)
  )
}