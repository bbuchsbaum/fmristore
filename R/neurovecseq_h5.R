#' Convert NeuroVecSeq to HDF5
#'
#' @description
#' Writes a NeuroVecSeq object (sequence of multiple 4D neuroimaging scans) to an HDF5 file.
#' Each scan is stored separately within the file, allowing for different time dimensions
#' while maintaining the same spatial dimensions across all scans.
#'
#' @param neurovecseq A NeuroVecSeq object containing multiple NeuroVec objects
#' @param file Path to the output HDF5 file. If NULL, uses a temporary file.
#' @param scan_names Optional character vector of names for each scan.
#'   If NULL, uses "scan_1", "scan_2", etc.
#' @param data_type Character string specifying the data type: "FLOAT" (default), "DOUBLE", or "INT"
#' @param chunk_dim Numeric vector specifying chunk dimensions for HDF5 storage.
#'   If NULL, uses time-optimized chunking (small spatial chunks, full time dimension).
#' @param compression Integer compression level 0-9. Default is 6.
#' @param scan_metadata Optional named list where each element is a list of metadata
#'   for the corresponding scan. Names should match scan_names.
#'
#' @return The file path of the created HDF5 file
#'
#' @details
#' The HDF5 file structure created is:
#' \itemize{
#'   \item \code{/} - Root group with attributes "rtype" = "NeuroVecSeq" and "n_scans"
#'   \item \code{/space/} - Group containing shared spatial information
#'   \item \code{/space/dim} - Spatial dimensions (3D)
#'   \item \code{/space/origin} - Spatial origin
#'   \item \code{/space/trans} - Spatial transformation matrix
#'   \item \code{/scans/} - Group containing all scans
#'   \item \code{/scans/<scan_name>/} - Group for each scan with "n_time" attribute
#'   \item \code{/scans/<scan_name>/data} - 4D data array for the scan
#'   \item \code{/scans/<scan_name>/metadata/} - Optional metadata group
#' }
#'
#' @examples
#' \dontrun{
#' # Create example NeuroVec objects with different time dimensions
#' vec1 <- NeuroVec(array(rnorm(10 * 10 * 5 * 20), dim = c(10, 10, 5, 20)),
#'   NeuroSpace(c(10, 10, 5, 20)))
#' vec2 <- NeuroVec(array(rnorm(10 * 10 * 5 * 30), dim = c(10, 10, 5, 30)),
#'   NeuroSpace(c(10, 10, 5, 30)))
#'
#' # Create NeuroVecSeq
#' nvs <- NeuroVecSeq(vec1, vec2)
#'
#' # Convert to HDF5
#' h5_file <- neurovecseq_to_h5(nvs,
#'   scan_names = c("rest_run1", "task_run1"),
#'   scan_metadata = list(
#'     rest_run1 = list(TR = 2.0, task = "rest"),
#'     task_run1 = list(TR = 2.0, task = "motor")
#' ))
#'
#' # File can later be read back (functionality to be implemented)
#' unlink(h5_file)
#' }
#'
#' @export
#' @importFrom hdf5r H5File h5attr
#' @importFrom neuroim2 space origin trans NeuroSpace
neurovecseq_to_h5 <- function(neurovecseq,
                              file = NULL,
                              scan_names = NULL,
                              data_type = "FLOAT",
                              chunk_dim = NULL,
                              compression = 6,
                              scan_metadata = NULL) {
  # Validate input
  if (!inherits(neurovecseq, "NeuroVecSeq")) {
    stop("Input must be a NeuroVecSeq object")
  }

  # Get list of NeuroVec objects
  vec_list <- neurovecseq@vecs
  n_scans <- length(vec_list)

  if (n_scans == 0) {
    stop("NeuroVecSeq contains no NeuroVec objects")
  }

  # Generate scan names if not provided
  if (is.null(scan_names)) {
    scan_names <- paste0("scan_", seq_len(n_scans))
  } else if (length(scan_names) != n_scans) {
    stop("Length of scan_names must match number of NeuroVec objects")
  }

  # Validate that all NeuroVec objects have same spatial dimensions
  first_space <- space(vec_list[[1]])
  for (i in seq_along(vec_list)) {
    if (!identical(dim(space(vec_list[[i]]))[1:3], dim(first_space)[1:3])) {
      stop("All NeuroVec objects must have identical spatial dimensions")
    }
  }

  # Set default file if not provided
  if (is.null(file)) {
    file <- tempfile(fileext = ".h5")
  }

  # Create/open HDF5 file
  h5obj <- hdf5r::H5File$new(file, mode = "w")
  on.exit(try(h5obj$close_all(), silent = TRUE), add = TRUE)

  # Write global attributes
  hdf5r::h5attr(h5obj, "rtype") <- "NeuroVecSeq"
  hdf5r::h5attr(h5obj, "n_scans") <- n_scans

  # Write shared spatial information from first NeuroVec
  space_grp <- h5obj$create_group("space")
  h5_write(h5obj, "/space/dim", dim(first_space)[1:3],
    dtype = hdf5r::h5types$H5T_NATIVE_INT32)
  h5_write(h5obj, "/space/origin", origin(first_space),
    dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE)
  h5_write(h5obj, "/space/trans", trans(first_space),
    dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE)

  # Create scans group
  scans_grp <- h5obj$create_group("scans")

  # Write each NeuroVec as a separate scan
  for (i in seq_along(vec_list)) {
    vec <- vec_list[[i]]
    scan_name <- scan_names[i]

    # Create scan group
    scan_grp <- scans_grp$create_group(scan_name)

    # Write time dimension as attribute
    time_dim <- dim(vec)[4]
    hdf5r::h5attr(scan_grp, "n_time") <- time_dim

    # Determine chunk dimensions if not provided
    scan_chunk_dim <- if (!is.null(chunk_dim)) {
      # Validate chunk dimensions
      vec_dims <- dim(vec)
      if (length(chunk_dim) != 4) {
        stop("chunk_dim must have exactly 4 elements")
      }
      if (any(chunk_dim > vec_dims)) {
        stop("chunk_dims must be <= dims of the data")
      }
      chunk_dim
    } else {
      # Time-optimized chunking by default
      c(10, 10, 10, time_dim)
    }

    # Write the 4D data
    dtype <- switch(data_type,
      "FLOAT" = hdf5r::h5types$H5T_IEEE_F32LE,
      "DOUBLE" = hdf5r::h5types$H5T_IEEE_F64LE,
      "INT" = hdf5r::h5types$H5T_NATIVE_INT32,
      stop("Unsupported data_type: ", data_type))

    h5_write(h5obj,
      path = file.path("/scans", scan_name, "data"),
      data = as.array(vec),
      dtype = dtype,
      chunk_dims = scan_chunk_dim,
      compression = compression)

    # Write scan-specific metadata if provided
    if (!is.null(scan_metadata) && scan_name %in% names(scan_metadata)) {
      meta <- scan_metadata[[scan_name]]
      if (is.list(meta) && length(meta) > 0) {
        meta_grp <- scan_grp$create_group("metadata")
        for (key in names(meta)) {
          h5_write(h5obj,
            path = file.path("/scans", scan_name, "metadata", key),
            data = meta[[key]],
            dtype = guess_h5_type(meta[[key]]))
        }
      }
    }
  }

  # Close the file to ensure all data is written
  h5obj$close_all()

  message("Successfully wrote NeuroVecSeq to: ", file)
  return(file)
}
