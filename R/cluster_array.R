#' @include all_class.R
#' @import hdf5r
#' @import neuroim2
#' @importFrom withr defer
#' @importFrom methods new is
#' @importFrom lifecycle deprecate_warn
#' @importFrom neuroim2 series
NULL




#' @rdname mask-methods
#' @export
setMethod("mask", "H5ClusteredArray", function(x) x@mask)

#' @rdname clusters-methods
#' @export
setMethod("clusters", "H5ClusteredArray", function(x) x@clusters)


#' @rdname h5file-methods
#' @export
setMethod("h5file", "H5ClusteredArray", function(x) x@obj)


#' Get Cluster Time Series by Mask Index (Internal Helper)
#'
#' @description
#' Internal helper function to retrieve time series data for specific voxels
#' (specified by mask-relative indices) from their corresponding cluster datasets
#' in an HDF5 file.
#'
#' This function is designed to work with objects inheriting from `H5ClusteredArray`.
#' It uses the `.dataset_path` generic method (which must be implemented by concrete subclasses)
#' to locate the correct HDF5 dataset for each cluster.
#'
#' @param x An object inheriting from `H5ClusteredArray`.
#' @param mask_indices An integer vector of voxel indices relative to the mask (1 to `x@n_voxels`).
#' @param time_indices An integer vector specifying which time points to retrieve. If `NULL`, all time points are retrieved.
#' @param n_time The total number of time points for the specific run/scan being accessed.
#'
#' @return A numeric matrix of shape `[length(mask_indices), length(time_indices)]`
#'         containing the requested time series data.
#' @keywords internal
#' @noRd
.get_cluster_timeseries_by_mask_index <- function(x, mask_indices, time_indices = NULL, n_time) {
    # 1. Validate spatial indices 'mask_indices' (relative to mask)
    nVoxMask <- x@n_voxels # Use stored voxel count from H5ClusteredArray
    if (any(mask_indices < 1 | mask_indices > nVoxMask)) {
      stop(sprintf("[.get_cluster_timeseries] Indices in 'mask_indices' are out of the valid mask range [1..%d]", nVoxMask))
    }
    mask_indices <- as.integer(mask_indices)
    n_request <- length(mask_indices)
    full_time_length <- as.integer(n_time) # Use provided n_time

    # 2. Validate optional time_indices
    if (!is.null(time_indices)) {
      if (!is.numeric(time_indices) || any(time_indices < 1 | time_indices > full_time_length)) {
        stop(sprintf("[.get_cluster_timeseries] 'time_indices' are out of the valid time range [1..%d]", full_time_length))
      }
      time_indices <- as.integer(time_indices)
      n_time_request <- length(time_indices)
    } else {
      time_indices <- seq_len(full_time_length) # Default to all time points
      n_time_request <- full_time_length
    }

    if (n_request == 0 || n_time_request == 0) {
        return(matrix(numeric(0), nrow = n_request, ncol = n_time_request))
    }

    # 3. Prepare result matrix
    result_mat <- matrix(NA_real_, nrow = n_request, ncol = n_time_request)

    # 4. Get cluster assignments for *all* requested indices ONCE
    # Accessing clusters slot from H5ClusteredArray
    clus_ids_req <- tryCatch(x@clusters@clusters[mask_indices], error = function(e) {
        stop(sprintf("[.get_cluster_timeseries] Failed to get cluster assignments for provided indices. Error: %s", e$message))
    })

    # 5. Group requested indices by their cluster ID
    index_groups <- split(seq_len(n_request), clus_ids_req)
    needed_cluster_ids <- names(index_groups)
    # Removed base_path calculation

    # 6. Loop through needed clusters and read data
    for (cid_str in needed_cluster_ids) {
        cid <- as.integer(cid_str) # Cluster ID should be integer
        row_indices_in_result <- index_groups[[cid_str]]
        mask_indices_this_cluster_req <- mask_indices[row_indices_in_result]

        ds <- NULL
        # === Use .dataset_path generic ===
        dset_path <- tryCatch(.dataset_path(x, cid),
                              error = function(e) {
                                  stop(sprintf("Failed to get dataset path for cluster %d. Subclass implementation error? Original error: %s", cid, e$message))
                              })
        # ==================================

        tryCatch({
            all_mask_indices_in_this_cluster <- which(x@clusters@clusters == cid)
            row_offsets_in_dataset <- match(mask_indices_this_cluster_req, all_mask_indices_in_this_cluster)

            if(any(is.na(row_offsets_in_dataset))) {
                stop(sprintf("Internal inconsistency: some requested voxels mapped to cluster %d but not found in its index list.", cid))
            }
            row_offsets_in_dataset <- as.integer(row_offsets_in_dataset)

            # Accessing obj slot from H5ClusteredArray
            if (!x@obj$exists(dset_path)){
                stop(sprintf("Dataset for cluster %d not found at path: %s", cid, dset_path))
            }
            ds <- x@obj[[dset_path]]
            # Simplified on.exit handling: rely on hdf5r GC or explicit closure later
            # on.exit(if (!is.null(ds) && inherits(ds, "H5D") && ds$is_valid) ds$close(), add = TRUE, after = FALSE)

            cluster_data_subset <- ds[row_offsets_in_dataset, time_indices, drop = FALSE]

            # Close dataset handle immediately after reading
            if (!is.null(ds) && inherits(ds, "H5D") && ds$is_valid) {
                 try(ds$close(), silent = TRUE)
            }

            expected_rows <- length(row_indices_in_result)
            if (!is.matrix(cluster_data_subset) ||
                 nrow(cluster_data_subset) != expected_rows ||
                 ncol(cluster_data_subset) != n_time_request) {
                 stop(sprintf("Read data dimensions [%s] mismatch expected [%d,%d] for cluster %d",
                              paste(dim(cluster_data_subset), collapse=","), expected_rows, n_time_request, cid))
            }

            result_mat[row_indices_in_result, ] <- cluster_data_subset

        }, error = function(e) {
            if (!is.null(ds) && inherits(ds, "H5D") && ds$is_valid) {
                try(ds$close(), silent = TRUE)
            }
            stop(sprintf("[.get_cluster_timeseries] Failed processing cluster %d at path '%s'. Original error: %s", cid, dset_path, e$message))
        })
    } # End loop over clusters

    # 7. Final check for NAs
    if (any(is.na(result_mat))) {
        warning("[.get_cluster_timeseries] Result matrix contains NAs. Data reading may have failed for some clusters or voxels within clusters.")
    }

    return(result_mat)
}




#' @keywords internal
#' @noRd
setMethod(".dataset_path", "H5ClusteredArray",
          function(x, cid, ...) {
            stop(sprintf("Internal Error: .dataset_path method not implemented for class '%s'. Subclass must provide a method.", class(x)[1]))
          })

#' Get HDF5 Dataset Path for H5ClusteredRunFull
#'
#' @description
#' Implementation of the `.dataset_path` generic for `H5ClusteredRunFull` objects.
#' Constructs the path to a specific cluster dataset within the run's group in the HDF5 file,
#' typically `/scans/<scan_name>/clusters/cluster_<cid>`.
#'
#' @param x An `H5ClusteredRunFull` object.
#' @param cid The cluster ID (integer).
#' @param ... Additional arguments (not used).
#'
#' @return A character string representing the HDF5 dataset path.
#' @keywords internal
#' @noRd
setMethod(".dataset_path", "H5ClusteredRunFull",
          function(x, cid) {
            # Basic validation
            if (!is.numeric(cid) || length(cid) != 1 || floor(cid) != cid || cid <= 0) {
                stop("['.dataset_path', H5ClusteredRunFull] Cluster ID 'cid' must be a single positive integer.")
            }
            if (!is.character(x@scan_name) || length(x@scan_name) != 1 || nchar(x@scan_name) == 0) {
                 stop("['.dataset_path', H5ClusteredRunFull] Invalid 'scan_name' slot.")
            }
            
            # Construct the standard path
            sprintf("/scans/%s/clusters/cluster_%d", x@scan_name, as.integer(cid))
          })



#' @keywords internal
#' @noRd
.subset_h5crunfull <- function(x, i, j, k, l, drop = TRUE) {

  # Inherits obj, mask, clusters, n_voxels from H5ClusteredArray
  # Uses scan_name, n_time from H5ClusteredRunFull
  dims   <- dim(x@mask)
  nvox   <- x@n_voxels
  nt     <- x@n_time # Use n_time from H5ClusteredRunFull

  # ---------- 1.  default / normalise inputs -----------------------
  if (missing(j) || is.null(j)) j <- NULL
  if (missing(k) || is.null(k)) k <- NULL
  if (missing(l) || is.null(l)) l <- seq_len(nt)

  # ---------- 2.  fast path : mask-based indexing ------------------
  if (is.null(j) && is.null(k)) {
    if (is.numeric(i) && all(i >= 1) && all(i <= nvox)) {
      # treat i as mask indices
      # Note: pass n_time explicitly to the helper
      dat <- .get_cluster_timeseries_by_mask_index(x, mask_indices = as.integer(i),
                                                 time_indices = l, n_time = nt)
      if (drop) return(drop(dat))
      return(dat)
    }
    # Fall through: interpret i as x-coordinate
  }

  # ---------- 3.  coordinate path ---------------------------------
  if (missing(i) || is.null(i)) i <- seq_len(dims[1])
  if (is.null(j)) j <- seq_len(dims[2])
  if (is.null(k)) k <- seq_len(dims[3])

  if (!is.numeric(i) || !is.numeric(j) || !is.numeric(k) || !is.numeric(l)) {
      stop("[',H5ClusteredRunFull'] Subscripts must be numeric or NULL.")
  }

  if (any(i < 1L | i > dims[1L]) ||
      any(j < 1L | j > dims[2L]) ||
      any(k < 1L | k > dims[3L]))
      stop("[',H5ClusteredRunFull'] Spatial subscript out of bounds")

  if (any(l < 1L | l > nt))
      stop("[',H5ClusteredRunFull'] Time subscript out of bounds")

  out <- array(0, dim = c(length(i), length(j), length(k), length(l)))

  grid <- as.matrix(expand.grid(i = i, j = j, k = k, KEEP.OUT.ATTRS = FALSE))
  lin  <- grid[,1] + (grid[,2]-1L)*dims[1] + (grid[,3]-1L)*dims[1]*dims[2]

  mask_array <- as.logical(as.array(x@mask))
  inside <- mask_array[lin]
  if (!any(inside)) return(if (drop) drop(out) else out)

  mask_lin <- which(mask_array)
  v_idx    <- match(lin[inside], mask_lin)

  # Note: pass n_time explicitly to the helper
  dat  <- .get_cluster_timeseries_by_mask_index(x, mask_indices = v_idx,
                                                time_indices = l, n_time = nt)

  idx <- grid[inside,, drop = FALSE]
  out_i <- match(idx[,1], i)
  out_j <- match(idx[,2], j)
  out_k <- match(idx[,3], k)

  for (t in seq_along(l))
    out[cbind(as.integer(out_i), as.integer(out_j), as.integer(out_k), t)] <- dat[, t]

  if (drop) drop(out) else out
}

#' Subset an H5ClusteredRunFull Object
#'
#' @description
#' Extracts data from an \code{H5ClusteredRunFull} object using array-like indexing.
#' Handles both coordinate-based and mask-index-based subsetting.
#'
#' @param x An \code{H5ClusteredRunFull} object.
#' @param i Row index (x-coordinate or mask index).
#' @param j Column index (y-coordinate).
#' @param k Slice index (z-coordinate).
#' @param l Time index.
#' @param ... Not used.
#' @param drop Logical. If \code{TRUE}, the result is coerced to the lowest possible dimension.
#'
#' @return An array or vector containing the subset of data.
#' @export
#' @family H5Clustered
setMethod("[",
  signature(x = "H5ClusteredRunFull", i = "ANY", j = "ANY", drop = "ANY"),
  function(x, i, j, k, l, ..., drop = TRUE)
      .subset_h5crunfull(x, i, j, k, l, drop))


#' @export
#' @family H5Clustered
setMethod("[",
  signature(x = "H5ClusteredRunFull", i = "ANY", j = "missing", drop = "ANY"),
  definition = function(x, i, ..., drop = TRUE) {
       # Handle cases: x[mask_indices], x[coords_matrix], x[i,j,k] (falls through if l missing)
       # Let the main helper sort it out.
      .subset_h5crunfull(x, i = i, j = NULL, k = NULL, l = NULL, drop = drop)
  })


#' @export
#' @family H5Clustered
#' @rdname dim-methods
setMethod("dim", "H5ClusteredRunFull",
  function(x) {
    n_time <- x@n_time
    if (is.null(n_time) || is.na(n_time) || length(n_time) != 1 || n_time <= 0 || floor(n_time) != n_time) {
      stop(sprintf("[dim,H5ClusteredRunFull] Slot 'n_time' (%s) is invalid.", as.character(n_time)))
    }
    c(dim(x@mask), n_time)
  }
)


#' @rdname series-methods
#' @export
#' @family H5Clustered
setMethod(
  f = "series",
  signature = signature(x = "H5ClusteredRunFull", i = "numeric"),
  definition = function(x, i, j, k, ...) {

    dims_mask <- dim(x@mask)
    n_vox_mask <- x@n_voxels
    nt <- x@n_time
    mask_indices_req <- NULL

    if (missing(j) && missing(k)) {
        # Case 1: i is numeric mask indices
        mask_indices_req <- as.integer(i)
        if (any(mask_indices_req < 1 | mask_indices_req > n_vox_mask)) {
            stop(sprintf("[series,H5ClusteredRunFull] Mask indices out of range [1..%d]", n_vox_mask))
        }
    } else if (!missing(j) && !missing(k)) {
        # Case 2: i, j, k are single coordinates
        if (length(i) != 1 || length(j) != 1 || length(k) != 1) {
             stop("[series,H5ClusteredRunFull] If providing i, j, k, they must be single values.")
        }
        coords <- matrix(c(i, j, k), nrow = 1, ncol = 3)
        # --- Convert coordinate to mask index ---
        mask_array <- as.logical(as.array(x@mask))
        global_mask_indices <- which(mask_array)
        if (any(coords[,1] < 1 | coords[,1] > dims_mask[1]) ||
            any(coords[,2] < 1 | coords[,2] > dims_mask[2]) ||
            any(coords[,3] < 1 | coords[,3] > dims_mask[3])) {
            stop("[series,H5ClusteredRunFull] Coordinates are out of bounds.")
        }
        linear_idx_full <- coords[,1] + (coords[,2]-1)*dims_mask[1] + (coords[,3]-1)*dims_mask[1]*dims_mask[2]
        if (!mask_array[linear_idx_full]) {
             warning("[series,H5ClusteredRunFull] Requested coordinate falls outside the mask.")
             return(matrix(numeric(0), nrow = nt, ncol = 0))
        }
        mask_indices_req <- match(linear_idx_full, global_mask_indices)
        if (is.na(mask_indices_req)) {
             stop("[series,H5ClusteredRunFull] Internal error: Failed to map valid coordinate to mask index.")
        }
        # --- End conversion ---
    } else {
        stop("[series,H5ClusteredRunFull] Invalid arguments. Provide numeric mask indices (i), a 3-col matrix (i), or single coordinates (i, j, k).")
    }

    if (is.null(mask_indices_req) || length(mask_indices_req) == 0) {
       return(matrix(numeric(0), nrow = nt, ncol = 0))
    }

    # Call helper with n_time
    voxel_data <- .get_cluster_timeseries_by_mask_index(x, mask_indices = mask_indices_req,
                                                        time_indices = NULL, n_time = nt)
    return(t(voxel_data))
  }
)

#' @rdname series-methods
#' @export
#' @family H5Clustered
setMethod(
  f = "series",
  signature = signature(x = "H5ClusteredRunFull", i = "matrix"),
  definition = function(x, i, ...) {
    coords <- i
    if (ncol(coords) != 3) {
        stop("[series,H5ClusteredRunFull] Coordinate matrix 'i' must have 3 columns.")
    }

    dims_mask <- dim(x@mask)
    n_vox_mask <- x@n_voxels
    nt <- x@n_time
    mask_array <- as.logical(as.array(x@mask))
    global_mask_indices <- which(mask_array)

    if (any(coords[,1] < 1 | coords[,1] > dims_mask[1]) ||
        any(coords[,2] < 1 | coords[,2] > dims_mask[2]) ||
        any(coords[,3] < 1 | coords[,3] > dims_mask[3])) {
        stop("[series,H5ClusteredRunFull] Coordinates are out of bounds.")
    }

    linear_idx_full <- coords[,1] + (coords[,2]-1)*dims_mask[1] + (coords[,3]-1)*dims_mask[1]*dims_mask[2]
    in_mask_subset <- mask_array[linear_idx_full]

    if (!all(in_mask_subset)) {
        warning("[series,H5ClusteredRunFull] Some coordinates fall outside mask and will be ignored.")
        coords <- coords[in_mask_subset, , drop = FALSE]
        linear_idx_full <- linear_idx_full[in_mask_subset]
        if (length(linear_idx_full) == 0) {
            return(matrix(numeric(0), nrow = nt, ncol = 0))
        }
    }

    mask_indices_req <- match(linear_idx_full, global_mask_indices)
    if (any(is.na(mask_indices_req))) {
         stop("[series,H5ClusteredRunFull] Internal error: Failed to map valid coordinates to mask indices.")
    }
    if (length(mask_indices_req) == 0) {
       return(matrix(numeric(0), nrow = nt, ncol = 0))
    }

    # Call helper with n_time
    voxel_data <- .get_cluster_timeseries_by_mask_index(x, mask_indices = mask_indices_req,
                                                        time_indices = NULL, n_time = nt)
    voxel_data_transposed <- t(voxel_data)

    if (nrow(coords) == ncol(voxel_data_transposed)) { # Check dimensions match
       colnames(voxel_data_transposed) <- apply(coords, 1, paste, collapse=",")
    }
    return(voxel_data_transposed)
  }
)


#' @importFrom neuroim2 linear_access NeuroSpace space 
#' @export
#' @family H5Clustered
#' @rdname linear_access-methods
#' @description Provides 4D linear access to data in an \code{H5ClusteredRunFull} object.
#' It reconstructs voxel values on the fly from the HDF5 file based on their cluster assignments.
#' Values for voxels outside the mask are returned as 0.
#'
#' @param x An \code{H5ClusteredRunFull} object.
#' @param i A numeric vector of 4D linear indices.
#' @param ... Additional arguments (not used for this method).
#'
#' @return A numeric vector of values corresponding to the provided linear indices.
#'
#' @examples
#' if (requireNamespace("neuroim2", quietly = TRUE) &&
#'     requireNamespace("hdf5r", quietly = TRUE) &&
#'     exists("H5ClusteredExperiment", where = "package:fmristore") &&
#'     exists("linear_access", where = "package:neuroim2") &&
#'     !is.null(fmristore:::create_minimal_h5_for_H5ClusteredExperiment)) {
#'
#'   temp_exp_file <- NULL
#'   exp_obj <- NULL
#'   run_full <- NULL
#' 
#'   tryCatch({
#'     # Create a minimal H5ClusteredExperiment
#'     temp_exp_file <- fmristore:::create_minimal_h5_for_H5ClusteredExperiment(
#'       master_mask_dims = c(3L, 3L, 2L), # Small dimensions
#'       num_master_clusters = 2L,
#'       n_time_run1 = 4L, # For Run1_Full
#'       n_time_run2 = 0   # No need for Run2_Summary here
#'     )
#'     exp_obj <- fmristore::H5ClusteredExperiment(file_path = temp_exp_file)
#'     
#'     # Access the H5ClusteredRunFull object (helper creates "Run1_Full")
#'     # The runs() method should give access to the list of runs
#'     available_runs <- runs(exp_obj)
#'     run_full <- available_runs[["Run1_Full"]] # Assuming helper creates this scan name
#'     
#'     if (!is.null(run_full)) {
#'       # Get dimensions: X, Y, Z, T
#'       run_dims <- dim(run_full) # Should be c(3,3,2,4)
#'       total_elements <- prod(run_dims)
#'       
#'       # Example: Access first 5 linear indices and last 5
#'       indices_to_access <- c(1:5, (total_elements-4):total_elements)
#'       # Ensure indices are within bounds if total_elements is small
#'       indices_to_access <- indices_to_access[indices_to_access <= total_elements & indices_to_access > 0]
#'       indices_to_access <- unique(indices_to_access)
#'       
#'       if (length(indices_to_access) > 0) {
#'          accessed_values <- neuroim2::linear_access(run_full, indices_to_access)
#'          cat("Accessed values for H5ClusteredRunFull:\n")
#'          print(accessed_values)
#'          cat("Number of values accessed:", length(accessed_values), "\n")
#'       } else {
#'          message("No valid indices to access for linear_access example.")
#'       }
#'     } else {
#'       message("Could not retrieve Run1_Full from the experiment for linear_access example.")
#'     }
#'     
#'   }, error = function(e) {
#'     message("linear_access example for H5ClusteredRunFull failed: ", e$message)
#'     if (!is.null(temp_exp_file)) message("Temporary file was: ", temp_exp_file)
#'   }, finally = {
#'     if (!is.null(exp_obj)) try(close(exp_obj), silent = TRUE)
#'     # run_full is part of exp_obj, its resources are managed by exp_obj$close()
#'     if (!is.null(temp_exp_file) && file.exists(temp_exp_file)) {
#'       unlink(temp_exp_file)
#'     }
#'   })
#' } else {
#'   message("Skipping linear_access H5ClusteredRunFull example: dependencies/helpers not available.")
#' }
setMethod(
  f = "linear_access",
  signature = signature(x="H5ClusteredRunFull", i="numeric"),
  definition = function(x, i, ...) {

    full_dims <- dim(x) # Uses dim method for H5ClusteredRunFull
    nt <- full_dims[4]
    n_request <- length(i)
    if (n_request == 0) return(numeric(0))

    max_elements <- prod(full_dims)
    if (any(i < 1 | i > max_elements)) {
        stop(sprintf("[linear_access,H5ClusteredRunFull] 4D indices out of range [1..%d]", max_elements))
    }

    coords4d <- arrayInd(i, .dim = full_dims)
    mask_array <- as.logical(as.array(x@mask))
    dims_mask <- dim(mask_array)
    global_mask_indices <- which(mask_array)

    coords3d <- coords4d[, 1:3, drop = FALSE]
    lin3d <- coords3d[,1] + (coords3d[,2]-1L)*dims_mask[1] + (coords3d[,3]-1L)*dims_mask[1]*dims_mask[2]
    in_mask_logical <- mask_array[lin3d]
    result_vec <- rep(NA_real_, n_request)

    if (any(in_mask_logical)) {
        valid_indices_in_i <- which(in_mask_logical)
        valid_coords4d <- coords4d[valid_indices_in_i, , drop = FALSE]
        valid_lin3d <- lin3d[valid_indices_in_i]
        mask_indices_req <- match(valid_lin3d, global_mask_indices)
        time_coords <- valid_coords4d[, 4]
        requests_by_time <- split(data.frame(orig_idx = valid_indices_in_i, mask_idx = mask_indices_req), time_coords)

        for (t_val_str in names(requests_by_time)) {
            t_val <- as.integer(t_val_str)
            current_requests <- requests_by_time[[t_val_str]]
            mask_indices_for_t <- current_requests$mask_idx
            original_indices_for_t <- current_requests$orig_idx

            # Call helper with n_time
            data_subset <- .get_cluster_timeseries_by_mask_index(x, mask_indices = mask_indices_for_t,
                                                                 time_indices = t_val, n_time = nt)
            result_vec[original_indices_for_t] <- data_subset[, 1]
        }
    }

    result_vec[is.na(result_vec)] <- 0
    return(result_vec)
  }
)


#' @rdname show-methods
#' @export
#' @family H5Clustered
setMethod(
  f = "show",
  signature = "H5ClusteredRunFull",
  definition = function(object) {
    cat("\n", crayon::bold(crayon::blue("H5ClusteredRunFull")), "\n", sep = "")
    cat(crayon::silver("────────────────────────────────────────\n"))
    cat(crayon::bold(crayon::yellow("Run Info")), "\n")
    cat(crayon::silver(" • "), crayon::green("Scan Name:"), object@scan_name, "\n")
    cat(crayon::silver(" • "), crayon::green("Time points:"), object@n_time, "\n")
    cat(crayon::bold("\nShared Info (from H5ClusteredArray)"), "\n")
    cat(crayon::silver(" • "), crayon::green("Active voxels in mask:"), object@n_voxels, "\n")
    if (!is.null(object@clusters) && length(object@clusters@clusters) > 0) {
        cluster_ids <- unique(object@clusters@clusters)
        n_clusters  <- length(cluster_ids)
        cat(crayon::silver(" • "), crayon::green("Number of clusters:"), n_clusters, "\n")
    } else {
        cat(crayon::silver(" • "), crayon::green("Number of clusters:"), "(NA)", "\n")
    }
    cat(crayon::bold("\nStorage: HDF5 File"), "\n")
    if (!is.null(object@obj) && inherits(object@obj, "H5File") && object@obj$is_valid) {
      cat(crayon::silver(" • "), "Path: ",
          crayon::magenta(object@obj$get_filename()), "\n", sep="")
      # Check if the specific run path exists
      scan_clusters_path <- tryCatch(.dataset_path(object, 1), # Use helper for base path logic
                                      error = function(e) NULL)
      if (!is.null(scan_clusters_path)) {
          scan_clusters_path <- dirname(scan_clusters_path) # Get parent dir /scans/../clusters
           cat(crayon::silver(" • "), "Run clusters path: ",
               crayon::magenta(scan_clusters_path), " (",
               ifelse(object@obj$exists(scan_clusters_path), crayon::green("exists"), crayon::red("missing")), ")\n", sep="")
      } else {
           cat(crayon::silver(" • "), "Run clusters path: unable to determine\n")
      }

    } else {
      cat(crayon::silver(" • "), "HDF5 file is ",
          crayon::red("INVALID or CLOSED"), "\n", sep="")
    }
    cat("\n")
  }
)


#' Constructor for H5ClusteredRunFull Objects
#'
#' @description
#' Creates a new `H5ClusteredRunFull` object, representing a single run of full
#' voxel-level clustered data from an HDF5 file.
#'
#' It performs necessary validation and can read metadata like `n_time` directly
#' from the HDF5 file attributes or metadata dataset.
#'
#' @param file_source Either an open `H5File` object or a character string path to the HDF5 file.
#' @param scan_name The character string name of the scan (e.g., "run1").
#' @param mask A `LogicalNeuroVol` object for the brain mask.
#' @param clusters A `ClusteredNeuroVol` object for cluster assignments.
#' @param n_time (Optional) The number of time points. If `NULL` (default), the function
#'   will attempt to read it from the HDF5 file attributes (e.g., `/scans/<scan_name>@n_time`
#'   or `/scans/<scan_name>/metadata/n_time`).
#' @param compress (Optional) Logical indicating compression status (metadata).
#'
#' @return A new `H5ClusteredRunFull` object.
#' @importFrom hdf5r H5File h5attr h5attr_names H5A 
#' @importFrom methods new is
#' @export
#' @family H5Clustered
make_run_full <- function(file_source, scan_name,
                            mask, clusters,
                            n_time = NULL, compress = FALSE) {

  lifecycle::deprecate_warn(
    when = "0.2.0", 
    what = "make_run_full()",
    with = "H5ClusteredRunFull()",
    details = "The make_run_* functions are deprecated in favor of direct H5* constructors."
  )
  
  # --- 1. Handle file source --- 
  fh <- open_h5(file_source, mode = "r")
  # Ensure file is closed if this function opened it.
  defer(if (fh$owns) fh$h5$close_all(), envir = parent.frame())
  h5obj <- fh$h5
  
  # --- 2. Validate inputs --- 
  if (!is.character(scan_name) || length(scan_name) != 1 || !nzchar(scan_name)) {
    stop("[make_run_full] 'scan_name' must be a non-empty character string.")
  }
  if (!is(mask, "LogicalNeuroVol")) {
    stop("[make_run_full] 'mask' must be a LogicalNeuroVol object.")
  }
  if (!is(clusters, "ClusteredNeuroVol")) {
     stop("[make_run_full] 'clusters' must be a ClusteredNeuroVol object.")
  }
  # Use helper to check dimensions 1-3
  check_same_dims(mask, clusters, dims_to_compare = 1:3, 
                  msg = "[make_run_full] Dimensions of 'mask' and 'clusters' must match.")
  
  n_vox <- sum(mask)
  if (!identical(length(clusters@clusters), as.integer(n_vox))) {
    stop(sprintf(
      "[make_run_full] Mismatch: clusters@clusters length (%d) != sum(mask) (%d).",
      length(clusters@clusters), n_vox
    ))
  }

  # --- 3. Determine n_time --- 
  if (is.null(n_time)) {
    # Try to read n_time from HDF5 attributes or metadata dataset
    scan_group_path <- paste0("/scans/", scan_name)
    scan_group <- NULL
    ds <- NULL # Ensure ds is NULL initially for cleanup
    
    # Check if scan group exists first, catching errors
    scan_group_exists <- tryCatch({
        h5obj$exists(scan_group_path)
    }, error = function(e) {
        warning(sprintf("[make_run_full] Suppressed HDF5 error during existence check for scan group '%s': %s", scan_group_path, conditionMessage(e)))
        FALSE
    })
    
    if (scan_group_exists) {
        tryCatch({
            scan_group <- h5obj[[scan_group_path]]
            on.exit(if(!is.null(scan_group) && inherits(scan_group, "H5Group") && scan_group$is_valid) try(scan_group$close()), add = TRUE)

            if ("n_time" %in% h5attr_names(scan_group)) {
                n_time <- h5attr(scan_group, "n_time")
            } else if (tryCatch(scan_group$exists("metadata/n_time"), error = function(e){ warning(sprintf("Suppressed HDF5 error checking metadata/n_time existence: %s", conditionMessage(e))); FALSE})) {
                 md_time <- NULL
                 tryCatch({
                     md_time <- scan_group[["metadata/n_time"]]
                     n_time <- md_time$read()
                 }, finally = {
                     if(!is.null(md_time) && inherits(md_time, "H5D") && md_time$is_valid) try(md_time$close())
                 })
            } else if (tryCatch(scan_group$exists("clusters"), error = function(e){ warning(sprintf("Suppressed HDF5 error checking clusters existence: %s", conditionMessage(e))); FALSE})) {
                # Try to infer from first cluster dataset
                first_cid <- clusters@clusters[1]
                if (!is.null(first_cid)) {
                    # Construct path directly, avoiding dummy object creation
                    dset_path_cid1 <- sprintf("/scans/%s/clusters/cluster_%d", scan_name, as.integer(first_cid))
                    
                    # Use .dataset_path to get the correct path for the first cluster
                    # dset_path_cid1 <- tryCatch(.dataset_path(new("H5ClusteredRunFull", scan_name=scan_name), first_cid), # Create dummy obj for path gen
                    #                        error=function(e) NULL)
                    if (!is.null(dset_path_cid1) && tryCatch(h5obj$exists(dset_path_cid1), error = function(e) {warning(sprintf("Suppressed HDF5 error checking %s: %s", dset_path_cid1, conditionMessage(e))); FALSE})) {
                        tryCatch({
                            ds <- h5obj[[dset_path_cid1]] # Open the dataset
                            dims <- ds$dims
                            if (length(dims) == 2) {
                                n_time <- dims[2]
                                message(sprintf("[make_run_full] Inferred n_time = %d from dataset '%s'.", n_time, dset_path_cid1))
                            } else {
                                stop(sprintf("[make_run_full] Dataset '%s' is malformed per specification: expected 2 dimensions, got %d.",
                                             dset_path_cid1, length(dims)))
                            }
                        }, error = function(e) {
                            warning(sprintf("[make_run_full] Error reading dimensions from '%s': %s", dset_path_cid1, e$message))
                        }, finally = {
                             if (!is.null(ds) && inherits(ds, "H5D") && ds$is_valid) try(ds$close(), silent = TRUE) # Close dataset handle
                        })
                    } else {
                       warning(sprintf("[make_run_full] Could not find or access dataset for first cluster (%d) at path '%s' to infer n_time.", first_cid, dset_path_cid1 %||% "(path generation failed)"))
                    }
                } else {
                     warning("[make_run_full] No clusters found in provided 'clusters' object to infer n_time.")
                }
            }
        }, error = function(e) {
            # Ensure scan_group is closed if error occurred after opening it
            if (!is.null(scan_group) && inherits(scan_group, "H5Group") && scan_group$is_valid) try(scan_group$close(), silent=TRUE)
            stop(sprintf("[make_run_full] Error accessing or validating scan group '%s'. Original error: %s", scan_group_path, e$message))
        })
    } else {
        warning(sprintf("[make_run_full] Scan group path '%s' not found. Cannot automatically determine n_time.", scan_group_path))
    }

    if (is.null(n_time)) {
      stop(sprintf("[make_run_full] Could not determine 'n_time' for scan '%s'. Provide it explicitly or ensure it exists in HDF5 attributes/metadata.", scan_name))
    }
  }

  # Final validation of n_time
  if (!is.numeric(n_time) || length(n_time) != 1 || n_time <= 0 || floor(n_time) != n_time) {
    stop(sprintf("[make_run_full] 'n_time' (%s) must be a single positive integer.", as.character(n_time)))
  }
  n_time <- as.integer(n_time)

  # --- 4. Create the object --- 
  run_obj <- new("H5ClusteredRunFull",
                 obj       = h5obj,
                 scan_name = scan_name,
                 mask      = mask,
                 clusters  = clusters,
                 n_voxels  = as.integer(n_vox),
                 n_time    = n_time,
                 compress  = as.logical(compress)
                )

  # --- 5. Important: Detach the object handle if we opened the file --- 
  # The returned object now owns the handle if fh$owns is FALSE.
  # If fh$owns is TRUE, the defer() handler will close it when the function exits.
  # We need to prevent the created object from trying to close a handle
  # that the calling function might still own or that defer will close.
  # If the user passed an open handle (fh$owns == FALSE), run_obj keeps it.
  # If we opened it (fh$owns == TRUE), run_obj's handle will be invalid after exit, 
  # which is the desired behavior unless the user manages the file handle lifecycle outside.
  if (fh$owns) {
    # If we opened the file, the created object should not manage the handle.
    # The validity checks within methods should handle this.
    # The 'defer' takes care of closing. We don't need to nullify run_obj@obj.
    # The key is that the caller understands the handle's lifecycle.
    
    # Alternative: If the design requires the returned object *always* has a valid 
    # handle *if* the source was a path, we'd need a different approach (e.g., 
    # clone the handle, manage externally). Sticking with simple defer for now.
  }

  # Return the created object
  return(run_obj)
}



#' @rdname as.matrix-methods
#' @export
#' @family H5Clustered
setMethod(
  f = "as.matrix",
  signature = signature(x="H5ClusteredRunSummary"),
  definition = function(x) {
    summary_grp <- NULL
    ds <- NULL
    mat_data <- NULL

    # Validate inputs from the object
    if (is.null(x@obj) || !x@obj$is_valid) stop("[as.matrix,H5ClusteredRunSummary] HDF5 file handle is invalid or closed.")
    if (!is.character(x@scan_name) || !nzchar(x@scan_name)) stop("[as.matrix,H5ClusteredRunSummary] Invalid 'scan_name' slot.")
    if (!is.character(x@summary_dset) || !nzchar(x@summary_dset)) stop("[as.matrix,H5ClusteredRunSummary] Invalid 'summary_dset' slot.")

    # Construct path to the specific summary dataset
    dset_path <- file.path("/scans", x@scan_name, "clusters_summary", x@summary_dset)

    tryCatch({
        if (!x@obj$exists(dset_path)) {
            stop(sprintf("[as.matrix,H5ClusteredRunSummary] Summary dataset not found at path: %s", dset_path))
        }
        ds <- x@obj[[dset_path]]
        on.exit(if (!is.null(ds) && inherits(ds, "H5D") && ds$is_valid) try(ds$close(), silent = TRUE), add = TRUE)

        mat_data <- ds$read() # Use $read() instead of []

    }, error = function(e) {
        # Ensure ds handle is closed on error
        if (!is.null(ds) && inherits(ds, "H5D") && ds$is_valid) {
            try(ds$close(), silent = TRUE)
        }
        stop(sprintf("[as.matrix,H5ClusteredRunSummary] Failed to read summary data for scan '%s' from '%s'. Original error: %s",
                     x@scan_name, dset_path, e$message))
    })

    if (is.null(mat_data)) {
        stop(sprintf("[as.matrix,H5ClusteredRunSummary] Failed to retrieve summary data matrix for scan '%s'. Result is NULL.", x@scan_name))
    }
    
    # Optional: Validate dimensions against n_time if available?
    # if (!is.na(x@n_time) && nrow(mat_data) != x@n_time) { ... warning ... }
    
    # Set column names if available and dimensions match
    if (length(x@cluster_names) > 0) {
        if (length(x@cluster_names) == ncol(mat_data)) {
            colnames(mat_data) <- x@cluster_names
        } else {
            warning(sprintf("[as.matrix,H5ClusteredRunSummary] Length of cluster_names (%d) != number of columns (%d) for scan '%s'. Names not set.",
                          length(x@cluster_names), ncol(mat_data), x@scan_name))
        }
    }

    mat_data
  }
)

#' @rdname as.data.frame-methods
#' @export
#' @family H5Clustered
setMethod(
  f = "as.data.frame",
  signature = signature(x="H5ClusteredRunSummary"),
  definition = function(x, row.names=NULL, optional=FALSE, ...) {
    # Calls the as.matrix method defined above for H5ClusteredRunSummary
    mat_data <- as.matrix(x)
    df <- as.data.frame(mat_data, row.names=row.names, optional=optional, ...)
    df
  }
)

# Override voxel-level accessors to prevent misuse


#' @export
#' @family H5Clustered
setMethod("[",
  signature(x = "H5ClusteredRunSummary", i = "ANY", j = "ANY", drop = "ANY"),
  function(x, i, j, k, l, ..., drop = TRUE) {
    stop("Voxel-level subsetting ([i,j,k,l]) is not available for H5ClusteredRunSummary objects. Use as.matrix() or as.data.frame() to access summary data.")
  })


#' @export
#' @family H5Clustered
setMethod("[",
  signature(x = "H5ClusteredRunSummary", i = "ANY", j = "missing", drop = "ANY"),
  definition = function(x, i, ..., drop = TRUE) {
    stop("Voxel-level subsetting (e.g., x[indices]) is not available for H5ClusteredRunSummary objects. Use as.matrix() or as.data.frame() to access summary data.")
  })

#' @rdname series-methods
#' @export
#' @family H5Clustered
setMethod(
  f = "series",
  signature = signature(x = "H5ClusteredRunSummary", i = "ANY"), # Catch numeric or matrix
  definition = function(x, i, ...) {
    stop("Voxel-level time series extraction (series()) is not available for H5ClusteredRunSummary objects. Use as.matrix() or as.data.frame() to access summary data.")
  }
)


#' @rdname linear_access-methods
#' @export
#' @family H5Clustered
#' @importFrom neuroim2 linear_access
setMethod(
  f = "linear_access",
  signature = signature(x="H5ClusteredRunSummary", i="numeric"),
  definition = function(x, i, ...) {
    stop("Voxel-level linear access (linear_access()) is not available for H5ClusteredRunSummary objects. Use as.matrix() or as.data.frame() to access summary data.")
  }
)

#' Constructor for H5ClusteredRunSummary Objects
#'
#' @description
#' Creates a new `H5ClusteredRunSummary` object, representing a single run of
#' summary cluster time-series data from an HDF5 file.
#'
#' It performs validation, checks for the existence of the summary dataset,
#' and determines `n_time` from the dataset dimensions.
#'
#' @param file_source Either an open `H5File` object or a character string path to the HDF5 file.
#' @param scan_name The character string name of the scan (e.g., "run1").
#' @param mask A `LogicalNeuroVol` object for the brain mask (used for reference/consistency).
#' @param clusters A `ClusteredNeuroVol` object for cluster assignments. This is required
#'   for consistency and reference, even though summary data is accessed.
#' @param cluster_names (Optional) Character vector of names for the clusters.
#' @param cluster_ids (Optional) Integer vector of IDs for the clusters.
#' @param summary_dset (Optional) The name of the dataset within the run's summary group
#'   (default: "summary_data").
#'
#' @return A new `H5ClusteredRunSummary` object.
#' @importFrom hdf5r H5File H5D
#' @importFrom methods new is
#' @importFrom lifecycle deprecate_warn
#' @export
#' @family H5Clustered
make_run_summary <- function(file_source, scan_name,
                               mask, clusters,
                               cluster_names = character(), cluster_ids = integer(),
                               summary_dset = "summary_data") {

  lifecycle::deprecate_warn(
    when = "0.2.0", 
    what = "make_run_summary()",
    with = "H5ClusteredRunSummary()",
    details = "The make_run_* functions are deprecated in favor of direct H5* constructors."
  )

  # --- 1. Handle file source ---
  fh <- open_h5(file_source, mode = "r")
  # Ensure file is closed if this function opened it.
  defer(if (fh$owns) fh$h5$close_all(), envir = parent.frame())
  h5obj <- fh$h5

  # --- 2. Validate common inputs ---
  if (!is.character(scan_name) || length(scan_name) != 1 || !nzchar(scan_name)) {
    stop("[make_run_summary] \'scan_name\' must be a non-empty character string.")
  }
  if (!is(mask, "LogicalNeuroVol")) {
    stop("[make_run_summary] \'mask\' must be a LogicalNeuroVol object.")
  }
  n_vox <- sum(mask)
  # Validate clusters if provided, but allow NULL
  if (!is.null(clusters)) {
      if (!is(clusters, "ClusteredNeuroVol")) {
        stop("[make_run_summary] \'clusters\' must be a ClusteredNeuroVol object.")
      }
      # Use helper to check dimensions 1-3
      check_same_dims(mask, clusters, dims_to_compare = 1:3,
                      msg = "[make_run_summary] Dimensions of 'mask' and provided 'clusters' must match.")
      
      if (!identical(length(clusters@clusters), as.integer(n_vox))) {
        stop(sprintf(
          "[make_run_summary] Mismatch: provided clusters@clusters length (%d) != sum(mask) (%d).",
          length(clusters@clusters), n_vox
        ))
      }
  } else {
       # If clusters is NULL, we still need n_voxels, which we got from sum(mask)
  }
  
  if (!is.character(summary_dset) || length(summary_dset) != 1 || !nzchar(summary_dset)) {
      stop("[make_run_summary] \'summary_dset\' must be a non-empty character string.")
  }

  # --- 3. Validate summary dataset and get n_time ---
  dset_path <- file.path("/scans", scan_name, "clusters_summary", summary_dset)
  n_time <- NA_integer_
  n_clusters_in_dset <- NA_integer_
  ds <- NULL

  # Check existence, catching potential HDF5 errors if path exists up to parent but not final object
  path_exists <- tryCatch({
      h5obj$exists(dset_path)
  }, error = function(e) {
      # If H5Lexists fails because the link itself is missing within an existing group,
      # treat it as FALSE for our purpose (dataset not found)
      # More specific error checking could be added here if needed
      warning(sprintf("[make_run_summary] Suppressed HDF5 error during existence check for '%s': %s", dset_path, conditionMessage(e)))
      FALSE
  })

  if (!path_exists) {
    stop(sprintf("[make_run_summary] Summary dataset not found at expected path: %s", dset_path))
  }

  tryCatch({
      ds <- h5obj[[dset_path]]
      dims <- ds$dims
      if (length(dims) != 2) {
          stop(sprintf("[make_run_summary] Summary dataset at \'%s\' must have 2 dimensions [nTime, nClusters], but found %d dimensions.",
                       dset_path, length(dims)))
      }
      n_time <- as.integer(dims[1])
      n_clusters_in_dset <- as.integer(dims[2])

      # Validate n_time
      if (is.na(n_time) || n_time <= 0) {
          stop(sprintf("[make_run_summary] Invalid time dimension (%d) read from summary dataset \'%s\'.", n_time, dset_path))
      }

      # Validate cluster info against dataset dimensions
      # Determine the source of truth for cluster IDs
      has_provided_ids <- length(cluster_ids) > 0
      has_clusters_object <- !is.null(clusters) && inherits(clusters, "ClusteredNeuroVol") && length(clusters@clusters) > 0
      
      # 1. Determine actual cluster IDs
      if (has_provided_ids) {
          actual_cluster_ids <- cluster_ids
      } else if (has_clusters_object) {
          actual_cluster_ids <- unique(clusters@clusters)
      } else {
          # No IDs provided and no clusters object: cannot determine IDs reliably
          actual_cluster_ids <- integer() 
          warning(sprintf("[make_run_summary] Cannot determine cluster IDs for scan '%s' without 'cluster_ids' or a valid 'clusters' object.", scan_name))
      }
      
      # 2. Determine cluster names: Prioritize provided, then generate if needed
      has_provided_names <- length(cluster_names) > 0
      final_cluster_names <- character()
      generated_names <- FALSE
      
      if (has_provided_names) {
          final_cluster_names <- cluster_names
      } else if (length(actual_cluster_ids) > 0) {
          # Generate names from IDs if names are missing but IDs are known
          final_cluster_names <- paste0("Clus_", actual_cluster_ids)
          # If we couldn't determine IDs above, we can't generate ID-based names
      } else if (n_clusters_in_dset > 0) { 
          # Fallback: Generate names like "Col_X" if dataset has columns but no names/IDs known
          final_cluster_names <- paste0("Col_", seq_len(n_clusters_in_dset))
          warning(sprintf("[make_run_summary] No cluster names or IDs provided/derivable for scan '%s'; generated default column names (Col_X).", scan_name))
          generated_names <- TRUE
      } else {
          # No names provided, no IDs known, dataset has 0 columns -> empty names is correct
          final_cluster_names <- character()
      }

      # 3. Check consistency of final names/IDs against dataset dimensions
      if (length(final_cluster_names) != n_clusters_in_dset) {
           warning(sprintf("[make_run_summary] Final number of cluster names (%d) does not match dataset columns (%d) for scan '%s'. Check inputs/clusters object/dataset.",
                        length(final_cluster_names), n_clusters_in_dset, scan_name))
           # If names mismatched and we *generated* Col_X names, this indicates an issue
           # If names were provided or derived from IDs, this warning stands.
      }
      
      # Check IDs consistency - only if we didn't just generate Col_X names (as IDs might be unknown)
      if (!generated_names && length(actual_cluster_ids) != n_clusters_in_dset) {
           warning(sprintf("[make_run_summary] Final number of cluster IDs (%d) does not match dataset columns (%d) for scan '%s'. Check inputs/clusters object/dataset.",
                         length(actual_cluster_ids), n_clusters_in_dset, scan_name))
            # If IDs mismatch, should we force sequential IDs? Risky. Let's warn.
           # For safety, if using generated Col_X names, ensure IDs are sequential too.
           # Let's assign default sequential IDs if names were generated OR if ID length mismatches
           if (generated_names || length(actual_cluster_ids) != n_clusters_in_dset) {
                actual_cluster_ids <- seq_len(n_clusters_in_dset)
                warning(sprintf("[make_run_summary] Assigning default sequential IDs (1..%d) for scan '%s' due to mismatch or generated names.", n_clusters_in_dset, scan_name))
           }
      }

      # If cluster_ids were initially empty and we *didn't* generate Col_X names, but derived IDs from cluster object
      # ensure we now have IDs matching the columns, generating if necessary
      if (!has_provided_ids && !generated_names && length(actual_cluster_ids) != n_clusters_in_dset && n_clusters_in_dset > 0) {
           actual_cluster_ids <- seq_len(n_clusters_in_dset)
           warning(sprintf("[make_run_summary] Generating default sequential IDs (1..%d) for scan '%s' as derived IDs didn't match dataset columns.", n_clusters_in_dset, scan_name))
      }
      
      # ... after generating default names ...
      if (generated_names && length(actual_cluster_ids) == 0 && n_clusters_in_dset > 0) {
           actual_cluster_ids <- seq_len(n_clusters_in_dset)
           warning(sprintf("[make_run_summary] Generated default IDs 1..%d for scan '%s' to match dataset columns.", n_clusters_in_dset, scan_name))
       }

  }, error = function(e) {
      stop(sprintf("[make_run_summary] Error accessing or validating summary dataset \'%s\'. Original error: %s", dset_path, e$message))
  }, finally = {
      if (!is.null(ds) && inherits(ds, "H5D") && ds$is_valid) {
          try(ds$close(), silent = TRUE)
      }
  })

  # --- 4. Create the object ---
  run_obj <- new("H5ClusteredRunSummary",
                 obj           = h5obj,
                 scan_name     = scan_name,
                 mask          = mask,
                 clusters      = clusters,
                 n_voxels      = as.integer(n_vox),
                 n_time        = n_time,
                 cluster_names = as.character(final_cluster_names),
                 cluster_ids   = as.integer(actual_cluster_ids),
                 summary_dset  = summary_dset
                )

  # Similar handle lifecycle considerations as make_run_full
  # If fh$owns is TRUE, defer() will close the handle. The returned object's
  # handle will become invalid upon function exit.

  return(run_obj)
}

#' Show Method for H5ClusteredRunSummary
#' @param object An H5ClusteredRunSummary object
#' @importFrom cli cli_h1 cli_h2 cli_li col_blue col_green col_yellow col_silver col_red col_magenta symbol cli_text cli_alert_info
#' @importFrom utils head
#' @importFrom methods show
#' @export
#' @family H5Clustered
#' @rdname show-methods
setMethod(
  f = "show",
  signature = "H5ClusteredRunSummary",
  definition = function(object) {

    cli::cli_h1(cli::col_blue("H5ClusteredRunSummary"))

    cli::cli_h2(cli::col_yellow("Run Summary Info"))
    cli::cli_li(items = c(
        "{cli::col_green(\'Scan Name\')} : {object@scan_name}",
        "{cli::col_green(\'Time points\')} : {object@n_time}",
        "{cli::col_green(\'Summary Dset\')} : {object@summary_dset}"
    ))

    # Display cluster names/IDs (show first few)
    n_clus_names <- length(object@cluster_names)
    n_clus_ids   <- length(object@cluster_ids)
    cli::cli_text("{cli::col_green(\'Cluster Names\')} ({n_clus_names}):")
    if (n_clus_names > 0) {
        name_preview <- utils::head(object@cluster_names, 5)
        cli::cli_alert_info("  {paste(name_preview, collapse=\', \')}{ifelse(n_clus_names > 5, \' ...\', \'\')}")
    } else {
        cli::cli_alert_info("  (No cluster names provided)")
    }
    cli::cli_text("{cli::col_green(\'Cluster IDs\')} ({n_clus_ids}):")
     if (n_clus_ids > 0) {
        id_preview <- utils::head(object@cluster_ids, 5)
        cli::cli_alert_info("  {paste(id_preview, collapse=\', \')}{ifelse(n_clus_ids > 5, \' ...\', \'\')}")
    } else {
        cli::cli_alert_info("  (No cluster IDs provided)")
    }


    cli::cli_h2(cli::col_yellow("Shared Info (from H5ClusteredArray)"))
    # Check if clusters slot is valid and has data
    has_clusters <- !is.null(object@clusters) && inherits(object@clusters, "ClusteredNeuroVol") && length(object@clusters@clusters) > 0
    n_clusters_from_slot <- if(has_clusters) length(unique(object@clusters@clusters)) else NA

    cli::cli_li(items = c(
        "{cli::col_green(\'Active voxels in mask\')} : {object@n_voxels}",
        "{cli::col_green(\'Number of clusters (map)\')} : {ifelse(has_clusters, n_clusters_from_slot, \'(NA or Map not loaded)\')}"
    ))

    cli::cli_h2(cli::col_yellow("Storage: HDF5 File"))
    h5_valid <- !is.null(object@obj) && inherits(object@obj, "H5File") && object@obj$is_valid
    if (h5_valid) {
      filename <- tryCatch(object@obj$get_filename(), error = function(e) "Error getting name")
      summary_path <- file.path("/scans", object@scan_name, "clusters_summary", object@summary_dset)
      path_exists <- tryCatch(object@obj$exists(summary_path), error = function(e) FALSE)
      cli::cli_li(items = c(
        "{cli::col_green(\'Path\')} : {cli::col_magenta(filename)}",
        "{cli::col_green(\'Status\')} : {cli::col_green(\'Valid and Open\')}",
        "{cli::col_green(\'Summary Dataset\')} : {cli::col_magenta(summary_path)} ({ifelse(path_exists, cli::col_green(cli::symbol$tick), cli::col_red(cli::symbol$cross))} exists)"
      ))
    } else {
      cli::cli_li(items = c(
        "{cli::col_green(\'Status\')} : {cli::col_red(\'INVALID or CLOSED\')}"
      ))
    }
    cli::cli_text("") # Add a newline at the end
  }
)
