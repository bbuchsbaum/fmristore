#' @rdname as_h5-methods
#' @export
setMethod(
  f = "as_h5",
  signature = signature(object = "NeuroVec"),
  definition = function(object, file = NULL, data_type = "FLOAT",
                         chunk_dim = c(4, 4, 4, dim(object)[4]),
                         compression = 6) {
    to_nih5_vec(object, file_name = file, data_type = data_type,
                chunk_dim = chunk_dim, compression = compression)
  }
)

#' @rdname as_h5-methods
#' @export
setMethod(
  f = "as_h5",
  signature = signature(object = "LatentNeuroVec"),
  definition = function(object, file = NULL, data_type = "FLOAT",
                         compression = 6) {
    to_h5_latentvec(object, file_name = file, data_type = data_type,
                    compression = compression)
  }
)

#' @rdname as_h5-methods
#' @export
setMethod(
  f = "as_h5",
  signature = signature(object = "LabeledVolumeSet"),
  definition = function(object, file, mask, labels, compression = 4,
                         dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE,
                         chunk_size = 1024, header_values = list()) {
    write_labeled_vec(vec = object, mask = mask, labels = labels, file = file,
                      compression = compression, dtype = dtype,
                      chunk_size = chunk_size, header_values = header_values)
  }
)

#' @rdname as_h5-methods
#' @export
setMethod(
  f = "as_h5",
  signature = signature(object = "list"),
  definition = function(object, file, scan_names, mask, clusters,
                         scan_metadata, cluster_metadata = NULL,
                         summary_only = FALSE, compression = 4,
                         chunk_size = 1024) {
    write_clustered_dataset(file = file, vecs = object, scan_names = scan_names,
                           mask = mask, clusters = clusters,
                           scan_metadata = scan_metadata,
                           cluster_metadata = cluster_metadata,
                           summary_only = summary_only,
                           compression = compression,
                           chunk_size = chunk_size)
  }
)

#' @rdname as_h5-methods
#' @export
setMethod(
  f = "as_h5",
  signature = signature(object = "NeuroVecSeq"),
  definition = function(object, file = NULL, ...) {
    # Extract additional arguments
    scan_names <- list(...)$scan_names
    data_type <- list(...)$data_type %||% "FLOAT"
    chunk_dim <- list(...)$chunk_dim
    compression <- list(...)$compression %||% 6
    scan_metadata <- list(...)$scan_metadata
    
    # Validate inputs
    if (is.null(file)) {
      file <- tempfile(fileext = ".h5")
    }
    
    # Get list of NeuroVec objects from the NeuroVecSeq
    vec_list <- object@vecs
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
    
    # Return an H5NeuroVecSeq object (to be implemented) or file path
    # For now, return the file path
    message("Successfully wrote NeuroVecSeq to: ", file)
    return(file)
  }
) 