#' Assert that an argument is a non-empty numeric vector
#' 
#' @param x The value to check.
#' @param arg The name of the argument being checked (as a string).
#' @param fn The name of the calling function (as a string) for error messages.
#' @return Invisible `NULL` if check passes, otherwise stops with an error.
#' @keywords internal
#' @noRd
assert_non_empty_numeric <- function(x, arg, fn) {
  if (is.null(x) || !is.numeric(x) || length(x) == 0L) {
      stop(sprintf("[%s] Argument '%s' must be a non-empty numeric vector.", 
                   fn %||% "unknown function", 
                   arg %||% "unknown argument"))
  }
  invisible(NULL)
}

#' Assert that two objects have identical dimensions (or a subset thereof)
#' 
#' @param a First object (matrix, array, or object with a dim() method).
#' @param b Second object.
#' @param dims_to_compare Numeric vector specifying which dimensions to compare. 
#'   If `NULL` (default), compares all dimensions.
#' @param msg Optional custom error message prefix.
#' @return Invisible `NULL` if dimensions match, otherwise stops with an error.
#' @keywords internal
#' @noRd
check_same_dims <- function(a, b, dims_to_compare = NULL, msg = NULL) {
  dim_a <- dim(a)
  dim_b <- dim(b)
  
  # Select dimensions to compare
  dims_a_sub <- if (is.null(dims_to_compare)) dim_a else dim_a[dims_to_compare]
  dims_b_sub <- if (is.null(dims_to_compare)) dim_b else dim_b[dims_to_compare]
  
  # Handle potential errors if dims_to_compare is out of bounds
  if (anyNA(dims_a_sub) || anyNA(dims_b_sub)) {
     stop("check_same_dims: 'dims_to_compare' contains indices out of bounds for one or both objects.")
  }

  if (!identical(dims_a_sub, dims_b_sub)) {
     # Format dimensions for error message (show full dims for context)
     dim_a_str <- if (is.null(dim_a)) paste0("[", length(a), "]") else paste(dim_a, collapse = "x")
     dim_b_str <- if (is.null(dim_b)) paste0("[", length(b), "]") else paste(dim_b, collapse = "x")
     compare_str <- if (is.null(dims_to_compare)) "all dims" else paste0("dims ", paste(dims_to_compare, collapse=","))
     
     stop(paste0(msg %||% "Dimension mismatch: ", 
                 "Object A dims [", dim_a_str, "] vs Object B dims [", dim_b_str, "]",
                 " (comparing ", compare_str, ")"))
  }
  invisible(NULL)
}

#' Validate that two objects have identical dimensions (or a subset thereof)
#' 
#' Wraps `check_same_dims` to return NULL on success or the error message on failure.
#' Useful within validation functions that collect multiple errors.
#' 
#' @param a First object (matrix, array, or object with a dim() method).
#' @param b Second object.
#' @param dims_to_compare Numeric vector specifying which dimensions to compare. 
#'   If `NULL` (default), compares all dimensions.
#' @param msg Optional custom error message prefix passed to `check_same_dims`.
#' @return `NULL` if dimensions match, otherwise the error message string from `check_same_dims`.
#' @keywords internal
#' @noRd
validate_same_dims <- function(a, b, dims_to_compare = NULL, msg = NULL) {
  result <- tryCatch({
    check_same_dims(a = a, b = b, dims_to_compare = dims_to_compare, msg = msg)
    NULL # Return NULL explicitly on success
  }, error = function(e) {
    conditionMessage(e) # Return the error message string on failure
  })
  return(result)
}

#' Helper for default value if NULL (copied from cluster_experiment.R temporarily, should be centralized)
#' 
#' @param a First object.
#' @param b Second object.
#' @return `b` if `a` is `NULL`, otherwise `a`.
#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a 
