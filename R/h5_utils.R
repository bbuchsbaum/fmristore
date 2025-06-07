#' @importFrom hdf5r H5File
NULL

#' Open an HDF5 file or use an existing handle
#'
#' Handles the boilerplate of checking if the input is a file path or an
#' existing H5File object, opening the file if necessary, and returning
#' a list containing the handle and a flag indicating ownership (whether
#' the function opened the file and is responsible for closing it).
#'
#' @param src A character string path to an HDF5 file or an existing `H5File` object.
#' @param mode The mode to open the file in (e.g., "r", "r+", "w"). Defaults to "r".
#' @return A list with two elements:
#'   - `h5`: The `H5File` object handle.
#'   - `owns`: A logical flag, `TRUE` if the function opened the file (src was a path),
#'             `FALSE` otherwise (src was already an `H5File` object).
#' @keywords internal
open_h5 <- function(src, mode = "r") {
  if (inherits(src, "H5File")) {
    # If it's already an H5File object, check if it's open
    if (!src$is_valid) {
        stop("open_h5: provided H5File object is not valid (likely closed).")
    }
    # User owns this handle, we didn't open it
    return(list(h5 = src, owns = FALSE))
  }

  if (!is.character(src) || length(src) != 1L || !nzchar(src)){
     stop("open_h5: 'src' must be a single, non-empty character string (file path) or an H5File object.")
  }

  # Only check for existence if opening in read-only mode.
  # Modes 'w', 'w+', 'a' are expected to create the file if it doesn't exist.
  if (mode == "r" && !file.exists(src)) {
      stop("open_h5: file path does not exist: ", src)
  }

  # Attempt to open the file
  h5 <- tryCatch({
      hdf5r::H5File$new(src, mode = mode)
  }, error = function(e) {
      stop("open_h5: failed to open HDF5 file '", src, "' in mode '", mode, "'. Original error: ", conditionMessage(e))
  })

  # We opened it, so we 'own' it (caller is responsible for closing)
  return(list(h5 = h5, owns = TRUE))
} 

#' Helper function to ensure a valid LogicalNeuroVol mask is available.
#' Reads from HDF5 if mask is NULL, otherwise validates the provided mask.
#'
#' @param mask A `LogicalNeuroVol` object, a file path string, or `NULL`.
#' @param h5 An open `H5File` handle (needed if mask is `NULL`).
#' @param space A `NeuroSpace` object representing the target space for validation.
#' @param path The HDF5 path where the mask should be read from (default: "/mask").
#' @return A validated `LogicalNeuroVol` object.
#' @keywords internal
ensure_mask <- function(mask, h5, space, path = "/mask") {
  assert_that(inherits(space, "NeuroSpace"), msg = "ensure_mask: 'space' must be a NeuroSpace object.")

  if (is.null(mask)) {
    if (is.null(h5) || !inherits(h5, "H5File") || !h5$is_valid) {
      stop("ensure_mask: 'h5' handle must be provided and valid when 'mask' is NULL.")
    }
    # Assume read_h5_mask_to_LogicalNeuroVol handles H5File validation and path existence
    # Use the internal helper function and pass the reference space
    m <- read_h5_mask_to_LogicalNeuroVol(h5, path, space)
  } else {
    # Validate provided mask type
    if (!is(mask, "LogicalNeuroVol")) {
      stop("Provided 'mask' must be a LogicalNeuroVol object.")
    }
    m <- mask
  }

  # Check spatial dimensions against the provided space using helper
  check_same_dims(m, space, dims_to_compare = 1:3,
                  msg = "Mask dimensions do not match space dimensions")

  # Ensure the mask's space matches the provided space object exactly
  if (!identical(space(m), space)) {
      warning("Mask's internal NeuroSpace object does not match the provided reference space. Consider updating the mask's space for consistency.")
      # Optionally, force the space to be identical?
      # space(m) <- space
      # For now, just warn, as dim check passed.
  }

  return(m)
}

#' Assert that an HDF5 path exists
#'
#' Checks for the existence of a dataset or group within an HDF5 file and
#' stops with a clear error message if the path is missing.
#'
#' @param h5 A valid open `H5File` object.
#' @param path The path to check within the file.
#' @param desc Optional description of what the path represents. This is used
#'   in the error message.
#' @return Invisible `TRUE` if the path exists.
#' @keywords internal
assert_h5_path <- function(h5, path, desc = "Path") {
  if (!inherits(h5, "H5File") || !h5$is_valid) {
    stop("assert_h5_path: 'h5' must be a valid and open H5File object.")
  }
  if (!is.character(path) || length(path) != 1 || !nzchar(path)) {
    stop("assert_h5_path: 'path' must be a single, non-empty character string.")
  }

  exists <- tryCatch(h5$exists(path), error = function(e) {
    stop(sprintf("assert_h5_path: Error checking existence of '%s': %s",
                 path, conditionMessage(e)))
  })

  if (!isTRUE(exists)) {
    fname <- tryCatch(h5$get_filename(), error = function(e) "<unknown>")
    stop(sprintf("%s not found at path '%s' in HDF5 file '%s'",
                 desc, path, fname))
  }

  invisible(TRUE)
}

#' Read data from an HDF5 dataset
#'
#' Safely reads data from a specified dataset path within an HDF5 file,
#' handling existence checks and optional return of NULL if missing.
#'
#' @param h5 An open `H5File` object handle.
#' @param path The full path to the dataset within the HDF5 file (e.g., "/data/my_dataset").
#' @param missing_ok Logical, if `TRUE`, returns `NULL` if the dataset doesn't exist. 
#'   If `FALSE` (default), stops with an error if the dataset is missing.
#' @param read_args A list of additional arguments passed to the HDF5 dataset's `$read()` method
#'   (e.g., for subsetting: `list(args = list(i = 1:10))`)
#'
#' @return The data read from the dataset, or `NULL` if `missing_ok=TRUE` and the dataset is not found.
#' @keywords internal
h5_read <- function(h5, path, missing_ok = FALSE, read_args = NULL) {
  if (!inherits(h5, "H5File") || !h5$is_valid) {
      stop("h5_read: 'h5' must be a valid and open H5File object.")
  }
  if (!is.character(path) || length(path) != 1 || !nzchar(path)) {
      stop("h5_read: 'path' must be a single, non-empty character string.")
  }
  
  path_exists <- tryCatch({
      h5$exists(path)
  }, error = function(e) {
      # Handle cases where checking existence itself fails (e.g., permission issue, intermediate group missing)
      warning(sprintf("h5_read: Suppressed HDF5 error during existence check for path '%s': %s", path, conditionMessage(e)))
      FALSE
  })
  
  if (!path_exists) {
    if (missing_ok) {
      return(NULL)
    } else {
      stop(sprintf("h5_read: Dataset or group not found at path '%s' in HDF5 file '%s'", 
                   path, h5$get_filename()))
    }
  }

  # Object exists, now try to open and read
  dset <- NULL
  data <- NULL
  tryCatch({
    # Check if it's a dataset before trying to read
    obj_info <- h5$link_info(path)
    if (obj_info$type != "H5L_TYPE_HARD") {
         # Could be a group or something else; h5[[path]] might work for groups but read() will fail.
         # Treat non-datasets as an error unless specifically handled.
         stop(sprintf("Object at path '%s' is not a dataset (type: %s).",
                      path, obj_info$type))
    }
    
    dset <- h5[[path]] # Open the dataset
    if (is.null(read_args)) {
        data <- dset$read()
    } else {
        # Pass arguments using do.call
        data <- do.call(dset$read, read_args)
    }
  }, error = function(e) {
    stop(sprintf("h5_read: Failed to read data from path '%s'. Original error: %s", 
                 path, conditionMessage(e)))
  }, finally = {
    # Ensure dataset handle is closed
    if (!is.null(dset) && inherits(dset, "H5D") && dset$is_valid) {
      close_h5_safely(dset)
    }
  })

  return(data)
}

#' Convenience wrapper to read a subset from an HDF5 dataset
#'
#' Opens the requested dataset, reads a subset using \code{$read()} and then
#' ensures the dataset handle is closed.  This avoids repetitive open/read/close
#' boilerplate when subsetting.
#'
#' @param h5 An open \code{H5File} handle.
#' @param path Character string path to the dataset within \code{h5}.
#' @param index Optional list of indices passed to \code{$read()} via the
#'   \code{index} argument.
#' @return The subset of data read from the dataset.
#' @keywords internal
h5_read_subset <- function(h5, path, index = NULL) {
  if (!inherits(h5, "H5File") || !h5$is_valid) {
    stop("h5_read_subset: 'h5' must be a valid open H5File object.")
  }
  if (!is.character(path) || length(path) != 1L || !nzchar(path)) {
    stop("h5_read_subset: 'path' must be a single, non-empty character string.")
  }

  if (!h5$exists(path)) {
    stop(sprintf("h5_read_subset: dataset '%s' not found in file '%s'",
                 path, h5$get_filename()))
  }

  dset <- NULL
  out <- NULL
  tryCatch({
    dset <- h5[[path]]
    read_args <- list()
    if (!is.null(index)) {
      if (!is.list(index)) {
        stop("h5_read_subset: 'index' must be a list when provided.")
      }
      read_args$index <- index
    }
    if (length(read_args) == 0L) {
      out <- dset$read()
    } else {
      # For subsetting, pass indices directly to read()
      if (!is.null(index)) {
        # Check dataset dimensions to handle 1D vs 2D indexing
        dims <- dset$dims
        if (length(dims) == 2 && length(index) == 1) {
          # For 2D datasets, if only row indices provided, select all columns
          # Use drop=FALSE to ensure matrix output
          out <- dset[index[[1]], , drop=FALSE]
        } else if (length(dims) == 2 && length(index) == 2) {
          # For 2D datasets with row and column indices
          out <- dset[index[[1]], index[[2]], drop=FALSE]
        } else {
          out <- do.call(dset$read, index)
        }
      } else {
        out <- do.call(dset$read, read_args)
      }
    }
  }, error = function(e) {
    stop(sprintf("h5_read_subset: failed reading subset from '%s': %s",
                 path, conditionMessage(e)))
  }, finally = {
    close_h5_safely(dset)
  })

  out
}

# -------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------

# Recursively create missing groups, preserving the leading "/"
ensure_h5_groups <- function(h5, group_path) {
  if (group_path == "/" || h5$exists(group_path)) return(invisible(TRUE))

  parts <- strsplit(sub("^/", "", group_path), "/", fixed = TRUE)[[1]]
  cur   <- ""                          # we build "/a", "/a/b", ...
  for (p in parts) {
    cur <- paste0(cur, "/", p)
    if (!h5$exists(cur)) h5$create_group(cur)
  }
  invisible(TRUE)
}

obj_dims <- function(x) {
  d <- dim(x)
  if (is.null(d)) if (length(x) == 1L) integer(0) else length(x) else d
}

# -------------------------------------------------------------------------
# Main writer
# -------------------------------------------------------------------------

h5_write <- function(h5, path, data,
                     dtype         = NULL,
                     chunk_dims    = NULL,
                     compression   = 0L,
                     overwrite     = FALSE,
                     create_parent = TRUE,
                     write_args    = list(),
                     verbose       = getOption("h5.write.verbose", FALSE)) {

  # -- 1. Validate --------------------------------------------------------
  if (!inherits(h5, "H5File")) {
    stop(sprintf("h5_write: 'h5' must be an H5File object, got class: %s", paste(class(h5), collapse=", ")))
  }
  if (!h5$is_valid) {
    stop("h5_write: H5File object is not valid (file may be closed)")
  }
  if (!is.character(path) || length(path) != 1L || !nzchar(path) || substr(path, 1, 1) != "/")
    stop("`path` must be an absolute HDF5 path (beginning with '/').")
  if (!is.numeric(compression) || compression < 0 || compression > 9)
    stop("`compression` must be an integer 0-9.")
  compression <- as.integer(compression)

  # -- 2. Ensure parent groups -------------------------------------------
  parent_path <- dirname(path)
  if (create_parent) ensure_h5_groups(h5, parent_path)
  else if (!h5$exists(parent_path))
    stop("Parent group '", parent_path, "' does not exist (create_parent = FALSE).")

  # -- 3. Overwrite logic -------------------------------------------------
  if (h5$exists(path)) {
    if (!overwrite) stop("Dataset already exists at '", path, "'.")
    h5$link_delete(path)
  }

  # -- 4. Type / dims / chunk heuristics ---------------------------------
  if (is.null(dtype)) dtype <- hdf5r::guess_dtype(data)
  if (!inherits(dtype, "H5T")) stop("`dtype` must be an 'H5T' object.")

  dims   <- obj_dims(data)
  scalar <- length(dims) == 0L
  if (scalar) { chunk_dims <- NULL; compression <- 0L }
  else if (compression > 0L && is.null(chunk_dims)) {
    ## simple heuristic: one chunk ≈ 4 MB
    dtype_size <- tryCatch(hdf5r::h5const$type_sizes[[dtype@name]], error = function(e) 8)
    nk         <- max(1L, ceiling(prod(dims) * dtype_size / (4 * 1024^2)))
    chunk_dims <- pmax(1L, floor(dims / nk))
  }

  # -- 5. Create + write --------------------------------------------------
  dset <- h5$create_dataset(name       = path,
                            dims       = dims,
                            dtype      = dtype,
                            chunk_dims = chunk_dims,
                            gzip_level = compression)

  if (length(write_args) == 0L) {
    dset$write(NULL, value = data)                # full write
  } else {
    if (is.null(write_args$robj)) write_args$robj <- data
    do.call(dset$write, write_args)
  }

  if (verbose) message("h5_write: wrote ", path)
  invisible(dset)                           # leave open; caller may close
}

# Internal helper to guess the HDF5 type from an R vector
# (already present in read_write_hdf5.R - should be centralized)
#' @keywords internal
guess_h5_type <- function(x) {
  # For vectors/arrays, check the base type
  base_type <- typeof(x)
  
  if (base_type == "character") {
    return(hdf5r::H5T_STRING$new(size = Inf))
  } else if (base_type == "integer") {
    return(hdf5r::h5types$H5T_NATIVE_INT32)
  } else if (base_type == "double") {
    return(hdf5r::h5types$H5T_NATIVE_DOUBLE)
  } else if (base_type == "logical") {
    # Using INT8 for logical as HBOOL can be problematic sometimes
    # return(hdf5r::h5types$H5T_NATIVE_HBOOL) 
    return(hdf5r::h5types$H5T_NATIVE_INT8) 
  }
  stop(sprintf("Unsupported R data type for guess_h5_type(): %s", base_type))
} 

#' Safely close an HDF5 file handle
#'
#' Checks if the object is a valid H5File handle before attempting
#' to close it, silencing potential errors from double-closing etc.
#' @param h5 An object, expected to be an H5File handle.
#' @keywords internal
safe_h5_close <- function(h5) {
  if (inherits(h5, "H5File") && h5$is_valid) {
    try(h5$close_all(), silent = TRUE)
  }
} 

# Robust dtype-to-NIfTI mapper
#' @keywords internal
#' @noRd
map_dtype <- function(h5t) {
  # Use tryCatch for safety as H5T methods might fail on invalid types
  sz <- tryCatch(h5t$get_size(), error = function(e) NA_integer_)
  cls_obj <- tryCatch(h5t$get_class(), error = function(e) NULL)

  if (is.null(cls_obj) || is.na(sz)) {
     warning("Could not reliably determine HDF5 type properties (class or size). Returning DT_UNKNOWN (0).")
     return(c(0L, 0L)) # DT_UNKNOWN, bitpix 0
  }

  # Use the H5T_CLASS object's equality comparison
  is_float <- cls_obj == hdf5r::h5const$H5T_FLOAT
  is_integer <- cls_obj == hdf5r::h5const$H5T_INTEGER
  # Close the class object handle once done
  try(cls_obj$close(), silent=TRUE)

  is_signed <- NA # Default
  if (is_integer) {
      # For integers, try check the sign property
      sign_const <- tryCatch(h5t$get_sign(), error = function(e) NA_integer_) # Use tryCatch
      if (is.na(sign_const)) {
          # This might happen if get_sign isn't applicable or fails
          warning("Could not determine sign for HDF5 integer type using get_sign(). Returning DT_UNKNOWN (0).")
          return(c(0L, 0L))
      }
      # H5T_SGN_NONE (0) means unsigned, H5T_SGN_2 (2) means signed (2's complement)
      is_signed <- (sign_const == hdf5r::h5const$H5T_SGN_2)
  } else if (is_float) {
      # Assume standard floats are signed
      is_signed <- TRUE
  } else {
      # Handle other classes (String, Compound, etc.) - currently map to unknown
      is_signed <- NA
  }

  if (is.na(is_signed)) {
     warning("Unhandled HDF5 type class or failed sign detection. Returning DT_UNKNOWN (0).")
     return(c(0L, 0L))
  }

  # Construct key based on size, float status, and signedness
  key <- paste(sz, is_float, is_signed)

  switch(key,
    # Ints (is_float = FALSE)
    "1 FALSE FALSE" = c(2L, 8L),    # DT_UINT8
    "1 FALSE TRUE"  = c(256L, 8L),  # DT_INT8
    "2 FALSE FALSE" = c(512L, 16L), # DT_UINT16
    "2 FALSE TRUE"  = c(4L, 16L),   # DT_INT16
    "4 FALSE FALSE" = c(768L, 32L), # DT_UINT32
    "4 FALSE TRUE"  = c(8L, 32L),   # DT_INT32
    "8 FALSE FALSE" = c(1280L, 64L),# DT_UINT64
    "8 FALSE TRUE"  = c(1024L, 64L),# DT_INT64
    # Floats (is_float = TRUE, assumed signed = TRUE)
    "4 TRUE TRUE"   = c(16L, 32L),  # DT_FLOAT32
    "8 TRUE TRUE"   = c(64L, 64L),  # DT_FLOAT64
    # Default case
    {
      warning(paste("Unhandled HDF5 type combination:", key, "(Size, IsFloat, IsSigned). Returning DT_UNKNOWN (0)."))
      c(0L, 0L) # DT_UNKNOWN
    }
  )
}

#' Evaluate code with an open HDF5 dataset, closing it afterwards
#'
#' Opens a dataset, passes it to a function, and ensures the dataset handle is
#' closed when finished.
#'
#' @param h5 An open `H5File` handle.
#' @param path Path to the dataset within the file.
#' @param FUN Function executed with the open dataset as its single argument.
#' @return The result of `FUN`.
#' @keywords internal
with_h5_dataset <- function(h5, path, FUN) {
  ds <- NULL
  on.exit(if (!is.null(ds)) close_h5_safely(ds), add = TRUE)
  ds <- h5[[path]]
  FUN(ds)
}

#' Obtain the dimensions of an HDF5 dataset
#'
#' @inheritParams with_h5_dataset
#' @return Integer vector of dataset dimensions.
#' @keywords internal
h5_dataset_dims <- function(h5, path) {
  with_h5_dataset(h5, path, function(ds) ds$dims)
}
