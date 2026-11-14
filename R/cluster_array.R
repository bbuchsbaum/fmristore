#' @include all_class.R
#' @importFrom hdf5r H5File is.h5file
#' @importFrom neuroim2 NeuroSpace space
#' @importFrom withr defer
#' @importFrom methods new is
#' @importFrom lifecycle deprecate_warn
#' @importFrom neuroim2 series
NULL




#' @rdname mask-methods
#' @export
setMethod("mask", "H5ParcellatedArray", function(x) x@mask)

#' @rdname clusters-methods
#' @export
setMethod("clusters", "H5ParcellatedArray", function(x) x@clusters)


#' @rdname h5file-methods
#' @export
setMethod("h5file", "H5ParcellatedArray", function(x) {
  if (is.null(x@obj)) {
    stop("H5File object is NULL")
  }
  x@obj$filename
})


#' Get Cluster Time Series by Mask Index (Internal Helper)
#'
#' @description
#' Internal helper function to retrieve time series data for specific voxels
#' (specified by mask-relative indices) from their corresponding cluster datasets
#' in an HDF5 file.
#'
#' This function is designed to work with objects inheriting from `H5ParcellatedArray`.
#' It uses the `.dataset_path` generic method (which must be implemented by concrete subclasses)
#' to locate the correct HDF5 dataset for each cluster.
#'
#' @param x An object inheriting from `H5ParcellatedArray`.
#' @param mask_indices An integer vector of voxel indices relative to the mask (1 to `x@n_voxels`).
#' @param time_indices An integer vector specifying which time points to retrieve.
#'   If `NULL`, all time points are retrieved.
#' @param n_time The total number of time points for the specific run/scan being accessed.
#'
#' @return A numeric matrix of shape `[length(mask_indices), length(time_indices)]`
#'         containing the requested time series data.
#' @keywords internal
#' @noRd
.get_cluster_data_cached <- function(obj, cluster_ids) {
  # Implementation would go here but this function seems incomplete
  stop("Function not implemented")
}

#' Pre-compute cluster index mappings for efficient access
#'
#' @description
#' Creates a mapping from cluster IDs to their corresponding voxel indices.
#' This avoids repeated calls to `which()` for each cluster during data retrieval.
#'
#' @param clusters_vector Integer vector of cluster assignments
#' @return Named list where names are cluster IDs and values are integer vectors of indices
#' @keywords internal
#' @noRd
.precompute_cluster_indices <- function(clusters_vector) {
  split(seq_along(clusters_vector), clusters_vector)
}

.get_cluster_timeseries_by_mask_index <- function(x, mask_indices, time_indices = NULL, n_time) {
  # 1. Validate spatial indices 'mask_indices' (relative to mask)
  nVoxMask <- x@n_voxels # Use stored voxel count from H5ParcellatedArray
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
  # Accessing clusters slot from H5ParcellatedArray
  clus_ids_req <- tryCatch(x@clusters@clusters[mask_indices], error = function(e) {
    stop(sprintf("[.get_cluster_timeseries] Failed to get cluster assignments for provided indices. Error: %s", e$message))
  })

  # 5. Group requested indices by their cluster ID
  index_groups <- split(seq_len(n_request), clus_ids_req)
  needed_cluster_ids <- names(index_groups)

  # 6. Pre-compute cluster index map once for all clusters
  cluster_index_map <- .precompute_cluster_indices(x@clusters@clusters)

  # 7. Loop through needed clusters and read data
  for (cid_str in needed_cluster_ids) {
    cid <- as.integer(cid_str) # Cluster ID should be integer
    row_indices_in_result <- index_groups[[cid_str]]
    mask_indices_this_cluster_req <- mask_indices[row_indices_in_result]

    # === Use .dataset_path generic ===
    dset_path <- tryCatch(.dataset_path(x, cid),
      error = function(e) {
        stop(sprintf("Failed to get dataset path for cluster %d. Subclass implementation error? Original error: %s", cid, e$message))
      }
    )
    # ==================================

    tryCatch(
      {
        # Use pre-computed cluster indices instead of which()
        all_mask_indices_in_this_cluster <- cluster_index_map[[cid_str]]
        if (is.null(all_mask_indices_in_this_cluster)) {
          stop(sprintf("Cluster %d not found in pre-computed index map", cid))
        }
        row_offsets_in_dataset <- match(mask_indices_this_cluster_req, all_mask_indices_in_this_cluster)

        if (any(is.na(row_offsets_in_dataset))) {
          stop(sprintf("Internal inconsistency: some requested voxels mapped to cluster %d but not found in its index list.", cid))
        }
        row_offsets_in_dataset <- as.integer(row_offsets_in_dataset)

        # If time_indices spans all time points, we can read without column subsetting
        if (identical(time_indices, seq_len(full_time_length))) {
          # Read only rows, all columns
          cluster_data_subset <- h5_read_subset(
            x@obj,
            dset_path,
            list(row_offsets_in_dataset)
          )
        } else {
          # Read specific rows and columns
          cluster_data_subset <- h5_read_subset(
            x@obj,
            dset_path,
            list(row_offsets_in_dataset, time_indices)
          )
        }

        expected_rows <- length(row_indices_in_result)
        if (!is.matrix(cluster_data_subset) ||
          nrow(cluster_data_subset) != expected_rows ||
          ncol(cluster_data_subset) != n_time_request) {
          stop(sprintf(
            "Read data dimensions [%s] mismatch expected [%d,%d] for cluster %d",
            paste(dim(cluster_data_subset), collapse = ","), expected_rows, n_time_request, cid
          ))
        }

        result_mat[row_indices_in_result, ] <- cluster_data_subset
      },
      error = function(e) {
        stop(sprintf("[.get_cluster_timeseries] Failed processing cluster %d at path '%s'. Original error: %s", cid, dset_path, e$message))
      }
    )
  } # End loop over clusters

  # 7. Final check for NAs
  if (any(is.na(result_mat))) {
    warning("[.get_cluster_timeseries] Result matrix contains NAs. Data reading may have failed for some clusters or voxels within clusters.")
  }

  return(result_mat)
}




#' @keywords internal
#' @noRd
setMethod(
  ".dataset_path", "H5ParcellatedArray",
  function(x, cid, ...) {
    stop(sprintf("Internal Error: .dataset_path method not implemented for class '%s'. Subclass must provide a method.", class(x)[1]))
  }
)

#' Get HDF5 Dataset Path for H5ParcellatedScan
#'
#' @description
#' Implementation of the `.dataset_path` generic for `H5ParcellatedScan` objects.
#' Constructs the path to a specific cluster dataset within the run's group in the HDF5 file,
#' typically `/scans/<scan_name>/clusters/cluster_<cid>`.
#'
#' @param x An `H5ParcellatedScan` object.
#' @param cid The cluster ID (integer).
#' @param ... Additional arguments (not used).
#'
#' @return A character string representing the HDF5 dataset path.
#' @keywords internal
#' @noRd
setMethod(
  ".dataset_path", "H5ParcellatedScan",
  function(x, cid) {
    # Basic validation
    if (!is.numeric(cid) || length(cid) != 1 || floor(cid) != cid || cid <= 0) {
      stop("['.dataset_path', H5ParcellatedScan] Cluster ID 'cid' must be a single positive integer.")
    }
    if (!is.character(x@scan_name) || length(x@scan_name) != 1 || nchar(x@scan_name) == 0) {
      stop("['.dataset_path', H5ParcellatedScan] Invalid 'scan_name' slot.")
    }

    # Construct the standard path using constants
    sprintf(H5_PATHS$CLUSTER_DSET_TPL, x@scan_name, as.integer(cid))
  }
)



#' @keywords internal
#' @noRd
.subset_h5crunfull <- function(x, i, j, k, l, drop = TRUE) {
  # Inherits obj, mask, clusters, n_voxels from H5ParcellatedArray
  # Uses scan_name, n_time from H5ParcellatedScan
  dims <- dim(x@mask)
  nvox <- x@n_voxels
  nt <- x@n_time # Use n_time from H5ParcellatedScan

  # ---------- 1.  default / normalise inputs -----------------------
  if (missing(j) || is.null(j)) j <- NULL
  if (missing(k) || is.null(k)) k <- NULL
  if (missing(l) || is.null(l)) l <- seq_len(nt)

  # ---------- 2.  fast path : mask-based indexing ------------------
  if (is.null(j) && is.null(k)) {
    if (is.numeric(i) && all(i >= 1) && all(i <= nvox)) {
      # treat i as mask indices
      # Note: pass n_time explicitly to the helper
      dat <- .get_cluster_timeseries_by_mask_index(x,
        mask_indices = as.integer(i),
        time_indices = l, n_time = nt
      )
      if (drop) {
        return(drop(dat))
      }
      return(dat)
    }
    # Fall through: interpret i as x-coordinate
  }

  # ---------- 3.  coordinate path ---------------------------------
  if (missing(i) || is.null(i)) i <- seq_len(dims[1])
  if (is.null(j)) j <- seq_len(dims[2])
  if (is.null(k)) k <- seq_len(dims[3])

  if (!is.numeric(i) || !is.numeric(j) || !is.numeric(k) || !is.numeric(l)) {
    stop("[',H5ParcellatedScan'] Subscripts must be numeric or NULL.")
  }

  if (any(i < 1L | i > dims[1L]) ||
    any(j < 1L | j > dims[2L]) ||
    any(k < 1L | k > dims[3L])) {
    stop("[',H5ParcellatedScan'] Spatial subscript out of bounds")
  }

  if (any(l < 1L | l > nt)) {
    stop("[',H5ParcellatedScan'] Time subscript out of bounds")
  }

  out <- array(0, dim = c(length(i), length(j), length(k), length(l)))

  grid <- as.matrix(expand.grid(i = i, j = j, k = k, KEEP.OUT.ATTRS = FALSE))
  lin <- grid[, 1] + (grid[, 2] - 1L) * dims[1] + (grid[, 3] - 1L) * dims[1] * dims[2]

  mask_array <- as.logical(as.array(x@mask))
  inside <- mask_array[lin]
  if (!any(inside)) {
    return(if (drop) drop(out) else out)
  }

  mask_lin <- which(mask_array)
  v_idx <- match(lin[inside], mask_lin)

  # Note: pass n_time explicitly to the helper
  dat <- .get_cluster_timeseries_by_mask_index(x,
    mask_indices = v_idx,
    time_indices = l, n_time = nt
  )

  idx <- grid[inside, , drop = FALSE]
  out_i <- match(idx[, 1], i)
  out_j <- match(idx[, 2], j)
  out_k <- match(idx[, 3], k)

  for (t in seq_along(l)) {
    out[cbind(as.integer(out_i), as.integer(out_j), as.integer(out_k), t)] <- dat[, t]
  }

  if (drop) drop(out) else out
}

#' Subset an H5ParcellatedScan Object
#'
#' @description
#' Extracts data from an \code{H5ParcellatedScan} object using array-like indexing.
#' Handles both coordinate-based and mask-index-based subsetting.
#'
#' @param x An \code{H5ParcellatedScan} object.
#' @param i Row index (x-coordinate or mask index).
#' @param j Column index (y-coordinate).
#' @param k Slice index (z-coordinate).
#' @param l Time index.
#' @param ... Not used.
#' @param drop Logical. If \code{TRUE}, the result is coerced to the lowest possible dimension.
#'
#' @return An array or vector containing the subset of data.
#' @export
#' @family H5Cluster
setMethod(
  "[",
  signature(x = "H5ParcellatedScan", i = "ANY", j = "ANY", drop = "ANY"),
  function(x, i, j, k, l, ..., drop = TRUE) {
    .subset_h5crunfull(x, i, j, k, l, drop)
  }
)


#' Extract data from H5ParcellatedScan with missing j parameter
#'
#' @description
#' Extracts data from an \code{H5ParcellatedScan} object when only the first index is provided.
#' Can handle mask indices, coordinate matrices, or pass through to coordinate-based indexing.
#'
#' @param x An \code{H5ParcellatedScan} object
#' @param i Row index (x-coordinate), mask indices, or coordinate matrix
#' @param j Missing (not provided)
#' @param k Slice index (z-coordinate) - optional, passed via ...
#' @param l Time index - optional, passed via ...
#' @param ... Additional arguments
#' @param drop Logical. If \code{TRUE}, the result is coerced to the lowest possible dimension
#'
#' @return An array or vector containing the subset of data
#'
#' @export
#' @family H5Cluster
#' @rdname extract-methods
setMethod("[",
  signature(x = "H5ParcellatedScan", i = "ANY", j = "missing", drop = "ANY"),
  definition = function(x, i, j, ..., drop = TRUE) {
    # x[i] (optionally with drop =): a single subscript, i.e. mask indices or
    # coordinates. Let the main helper sort it out. nargs() counts x and every
    # subscript, blank ones included, so x[i] gives 2 and x[i, , , l] gives 5.
    n_args <- nargs() - as.integer(!missing(drop))
    if (n_args <= 2) {
      return(.subset_h5crunfull(x, i = i, j = NULL, k = NULL, l = NULL, drop = drop))
    }
    # x[i, , k, l] and similar: j is blank, so it spans the whole y axis.
    # k and l arrive through `...` and may be blank too.
    dims <- dim(x@mask)
    kl <- .blank_subscripts_kl(...)
    .subset_h5crunfull(x,
      i = if (missing(i)) seq_len(dims[1]) else i,
      j = seq_len(dims[2]),
      k = kl$k, l = kl$l, drop = drop
    )
  }
)

# Collect the k and l subscripts forwarded through `...`; blank ones become NULL.
.blank_subscripts_kl <- function(k, l, ...) {
  list(
    k = if (missing(k)) NULL else k,
    l = if (missing(l)) NULL else l
  )
}


#' Dimensions of an H5ParcellatedScan
#'
#' Returns the 4D dimensions of the scan: the three spatial dimensions of the
#' mask followed by the number of time points.
#'
#' @param x An \code{H5ParcellatedScan}.
#' @return An integer vector \code{c(x, y, z, time)}.
#' @export
#' @family H5Cluster
#' @rdname dim-methods
setMethod(
  "dim", "H5ParcellatedScan",
  function(x) {
    n_time <- x@n_time
    if (is.null(n_time) || is.na(n_time) || length(n_time) != 1 || n_time <= 0 || floor(n_time) != n_time) {
      stop(sprintf("[dim,H5ParcellatedScan] Slot 'n_time' (%s) is invalid.", as.character(n_time)))
    }
    c(dim(x@mask), n_time)
  }
)


#' @rdname series-methods
#' @export
#' @family H5Cluster
setMethod(
  f = "series",
  signature = signature(x = "H5ParcellatedScan", i = "numeric"),
  definition = function(x, i, j, k, ...) {
    dims_mask <- dim(x@mask)
    n_vox_mask <- x@n_voxels
    nt <- x@n_time
    mask_indices_req <- NULL

    if (missing(j) && missing(k)) {
      # Case 1: i is numeric mask indices
      mask_indices_req <- as.integer(i)
      if (any(mask_indices_req < 1 | mask_indices_req > n_vox_mask)) {
        stop(sprintf("[series,H5ParcellatedScan] Mask indices out of range [1..%d]", n_vox_mask))
      }
    } else if (!missing(j) && !missing(k)) {
      # Case 2: i, j, k are single coordinates
      if (length(i) != 1 || length(j) != 1 || length(k) != 1) {
        stop("[series,H5ParcellatedScan] If providing i, j, k, they must be single values.")
      }
      coords <- matrix(c(i, j, k), nrow = 1, ncol = 3)
      # --- Convert coordinate to mask index ---
      mask_array <- as.logical(as.array(x@mask))
      global_mask_indices <- which(mask_array)
      if (any(coords[, 1] < 1 | coords[, 1] > dims_mask[1]) ||
        any(coords[, 2] < 1 | coords[, 2] > dims_mask[2]) ||
        any(coords[, 3] < 1 | coords[, 3] > dims_mask[3])) {
        stop("[series,H5ParcellatedScan] Coordinates are out of bounds.")
      }
      linear_idx_full <- coords[, 1] + (coords[, 2] - 1) * dims_mask[1] + (coords[, 3] - 1) * dims_mask[1] * dims_mask[2]
      if (!mask_array[linear_idx_full]) {
        warning("[series,H5ParcellatedScan] Requested coordinate falls outside the mask.")
        return(matrix(numeric(0), nrow = nt, ncol = 0))
      }
      mask_indices_req <- match(linear_idx_full, global_mask_indices)
      if (is.na(mask_indices_req)) {
        stop("[series,H5ParcellatedScan] Internal error: Failed to map valid coordinate to mask index.")
      }
      # --- End conversion ---
    } else {
      stop("[series,H5ParcellatedScan] Invalid arguments. Provide numeric mask indices (i), a 3-col matrix (i), or single coordinates (i, j, k).")
    }

    if (is.null(mask_indices_req) || length(mask_indices_req) == 0) {
      return(matrix(numeric(0), nrow = nt, ncol = 0))
    }

    # Call helper with n_time
    voxel_data <- .get_cluster_timeseries_by_mask_index(x,
      mask_indices = mask_indices_req,
      time_indices = NULL, n_time = nt
    )
    return(t(voxel_data))
  }
)

#' @rdname series-methods
#' @export
#' @family H5Cluster
setMethod(
  f = "series",
  signature = signature(x = "H5ParcellatedScan", i = "matrix"),
  definition = function(x, i, ...) {
    coords <- i
    if (ncol(coords) != 3) {
      stop("[series,H5ParcellatedScan] Coordinate matrix 'i' must have 3 columns.")
    }

    dims_mask <- dim(x@mask)
    n_vox_mask <- x@n_voxels
    nt <- x@n_time
    mask_array <- as.logical(as.array(x@mask))
    global_mask_indices <- which(mask_array)

    if (any(coords[, 1] < 1 | coords[, 1] > dims_mask[1]) ||
      any(coords[, 2] < 1 | coords[, 2] > dims_mask[2]) ||
      any(coords[, 3] < 1 | coords[, 3] > dims_mask[3])) {
      stop("[series,H5ParcellatedScan] Coordinates are out of bounds.")
    }

    linear_idx_full <- coords[, 1] + (coords[, 2] - 1) * dims_mask[1] + (coords[, 3] - 1) * dims_mask[1] * dims_mask[2]
    in_mask_subset <- mask_array[linear_idx_full]

    if (!all(in_mask_subset)) {
      warning("[series,H5ParcellatedScan] Some coordinates fall outside mask and will be ignored.")
      coords <- coords[in_mask_subset, , drop = FALSE]
      linear_idx_full <- linear_idx_full[in_mask_subset]
      if (length(linear_idx_full) == 0) {
        return(matrix(numeric(0), nrow = nt, ncol = 0))
      }
    }

    mask_indices_req <- match(linear_idx_full, global_mask_indices)
    if (any(is.na(mask_indices_req))) {
      stop("[series,H5ParcellatedScan] Internal error: Failed to map valid coordinates to mask indices.")
    }
    if (length(mask_indices_req) == 0) {
      return(matrix(numeric(0), nrow = nt, ncol = 0))
    }

    # Call helper with n_time
    voxel_data <- .get_cluster_timeseries_by_mask_index(x,
      mask_indices = mask_indices_req,
      time_indices = NULL, n_time = nt
    )
    voxel_data_transposed <- t(voxel_data)

    if (nrow(coords) == ncol(voxel_data_transposed)) { # Check dimensions match
      colnames(voxel_data_transposed) <- apply(coords, 1, paste, collapse = ",")
    }
    return(voxel_data_transposed)
  }
)


#' @importFrom neuroim2 linear_access NeuroSpace space
#' @export
#' @family H5Cluster
#' @rdname linear_access-methods
#' @description Provides 4D linear access to data in an \code{H5ParcellatedScan} object.
#' It reconstructs voxel values on the fly from the HDF5 file based on their cluster assignments.
#' Values for voxels outside the mask are returned as 0.
#'
#' @param x An \code{H5ParcellatedScan} object.
#' @param i A numeric vector of 4D linear indices.
#' @param ... Additional arguments (not used for this method).
#'
#' @return A numeric vector of values corresponding to the provided linear indices.
#'
#' @examples
#' \dontrun{
#' if (requireNamespace("neuroim2", quietly = TRUE) &&
#'   requireNamespace("hdf5r", quietly = TRUE) &&
#'   exists("H5ParcellatedMultiScan", where = "package:fmristore") &&
#'   exists("linear_access", where = "package:neuroim2") &&
#'   !is.null(fmristore:::create_minimal_h5_for_H5ParcellatedMultiScan)) {
#'   temp_exp_file <- NULL
#'   exp_obj <- NULL
#'   run_full <- NULL
#'
#'   tryCatch({
#'     # Create a minimal H5ParcellatedMultiScan
#'     temp_exp_file <- fmristore:::create_minimal_h5_for_H5ParcellatedMultiScan(
#'       master_mask_dims = c(3L, 3L, 2L), # Small dimensions
#'       num_master_clusters = 2L,
#'       n_time_run1 = 4L, # For Run1_Full
#'       n_time_run2 = 0 # No need for Run2_Summary here
#'     )
#'     exp_obj <- fmristore::H5ParcellatedMultiScan(file_path = temp_exp_file)
#'
#'     # Access the H5ParcellatedScan object (helper creates "Run1_Full")
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
#'       indices_to_access <- c(1:5, (total_elements - 4):total_elements)
#'       # Ensure indices are within bounds if total_elements is small
#'       indices_to_access <- indices_to_access[indices_to_access <= total_elements &
#'         indices_to_access > 0]
#'       indices_to_access <- unique(indices_to_access)
#'
#'       if (length(indices_to_access) > 0) {
#'         accessed_values <- neuroim2::linear_access(run_full, indices_to_access)
#'         cat("Accessed values for H5ParcellatedScan:\n")
#'         print(accessed_values)
#'         cat("Number of values accessed:", length(accessed_values), "\n")
#'       } else {
#'         message("No valid indices to access for linear_access example.")
#'       }
#'     } else {
#'       message("Could not retrieve Run1_Full from the experiment for linear_access example.")
#'     }
#'   }, error = function(e) {
#'     message("linear_access example for H5ParcellatedScan failed: ", e$message)
#'     if (!is.null(temp_exp_file)) message("Temporary file was: ", temp_exp_file)
#'   }, finally = {
#'     if (!is.null(exp_obj)) try(close(exp_obj), silent = TRUE)
#'     # run_full is part of exp_obj, its resources are managed by exp_obj$close()
#'     if (!is.null(temp_exp_file) && file.exists(temp_exp_file)) {
#'       unlink(temp_exp_file)
#'     }
#'   })
#' } else {
#'   message("Skipping linear_access H5ParcellatedScan example: dependencies/helpers not available.")
#' }
#' }
setMethod(
  f = "linear_access",
  signature = signature(x = "H5ParcellatedScan", i = "numeric"),
  definition = function(x, i, ...) {
    full_dims <- dim(x) # Uses dim method for H5ParcellatedScan
    nt <- full_dims[4]
    n_request <- length(i)
    if (n_request == 0) {
      return(numeric(0))
    }

    max_elements <- prod(full_dims)
    if (any(i < 1 | i > max_elements)) {
      stop(sprintf("[linear_access,H5ParcellatedScan] 4D indices out of range [1..%d]", max_elements))
    }

    coords4d <- arrayInd(i, .dim = full_dims)
    mask_array <- as.logical(as.array(x@mask))
    dims_mask <- dim(mask_array)
    global_mask_indices <- which(mask_array)

    coords3d <- coords4d[, 1:3, drop = FALSE]
    lin3d <- coords3d[, 1] + (coords3d[, 2] - 1L) * dims_mask[1] + (coords3d[, 3] - 1L) * dims_mask[1] * dims_mask[2]
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
        data_subset <- .get_cluster_timeseries_by_mask_index(x,
          mask_indices = mask_indices_for_t,
          time_indices = t_val, n_time = nt
        )
        result_vec[original_indices_for_t] <- data_subset[, 1]
      }
    }

    result_vec[is.na(result_vec)] <- 0
    return(result_vec)
  }
)


#' @rdname show-methods
#' @export
#' @family H5Cluster
setMethod(
  f = "show",
  signature = "H5ParcellatedScan",
  definition = function(object) {
    cat("\n", cli::style_bold(cli::col_blue("H5ParcellatedScan")), "\n", sep = "")
    cat(cli::col_silver("----------------------------------------\n"))
    cat(cli::style_bold(cli::col_yellow("Run Info")), "\n")
    cat(cli::col_silver(" * "), cli::col_green("Scan Name:"), object@scan_name, "\n")
    cat(cli::col_silver(" * "), cli::col_green("Time points:"), object@n_time, "\n")
    cat(cli::style_bold("\nShared Info (from H5ParcellatedArray)"), "\n")
    cat(cli::col_silver(" * "), cli::col_green("Active voxels in mask:"), object@n_voxels, "\n")
    if (!is.null(object@clusters) && length(object@clusters@clusters) > 0) {
      cluster_ids <- unique(object@clusters@clusters)
      n_clusters <- length(cluster_ids)
      cat(cli::col_silver(" * "), cli::col_green("Number of clusters:"), n_clusters, "\n")
    } else {
      cat(cli::col_silver(" * "), cli::col_green("Number of clusters:"), "(NA)", "\n")
    }
    cat(cli::style_bold("\nStorage: HDF5 File"), "\n")
    if (!is.null(object@obj) && inherits(object@obj, "H5File") && object@obj$is_valid) {
      cat(cli::col_silver(" * "), "Path: ",
        cli::col_magenta(object@obj$get_filename()), "\n",
        sep = ""
      )
      scan_clusters_path <- tryCatch(.dataset_path(object, 1),
        error = function(e) NULL
      )
      if (!is.null(scan_clusters_path)) {
        scan_clusters_path <- dirname(scan_clusters_path)
        cat(cli::col_silver(" * "), "Run clusters path: ",
          cli::col_magenta(scan_clusters_path), " (",
          ifelse(object@obj$exists(scan_clusters_path), cli::col_green("exists"), cli::col_red("missing")), ")\n",
          sep = ""
        )
      } else {
        cat(cli::col_silver(" * "), "Run clusters path: unable to determine\n")
      }
    } else {
      cat(cli::col_silver(" * "), "HDF5 file is ",
        cli::col_red("INVALID or CLOSED"), "\n",
        sep = ""
      )
    }
    cat("\n")
  }
)


#' @rdname as.matrix-methods
#' @export
#' @family H5Cluster
setMethod(
  f = "as.matrix",
  signature = signature(x = "H5ParcellatedScanSummary"),
  definition = function(x) {
    summary_grp <- NULL
    ds <- NULL
    mat_data <- NULL

    # Validate inputs from the object
    if (is.null(x@obj) || !x@obj$is_valid) stop("[as.matrix,H5ParcellatedScanSummary] HDF5 file handle is invalid or closed.")
    if (!is.character(x@scan_name) || !nzchar(x@scan_name)) stop("[as.matrix,H5ParcellatedScanSummary] Invalid 'scan_name' slot.")
    if (!is.character(x@summary_dset) || !nzchar(x@summary_dset)) stop("[as.matrix,H5ParcellatedScanSummary] Invalid 'summary_dset' slot.")

    # Construct path to the specific summary dataset
    dset_path <- sprintf(H5_PATHS$SUMMARY_DSET_TPL, x@scan_name, x@summary_dset)

    tryCatch(
      {
        assert_h5_path(
          x@obj, dset_path,
          "[as.matrix,H5ParcellatedScanSummary] Summary dataset"
        )
        ds <- x@obj[[dset_path]]
        on.exit(if (!is.null(ds) && inherits(ds, "H5D") && ds$is_valid) close_h5_safely(ds), add = TRUE)

        mat_data <- ds$read() # Use $read() instead of []
      },
      error = function(e) {
        # Ensure ds handle is closed on error
        if (!is.null(ds) && inherits(ds, "H5D") && ds$is_valid) {
          close_h5_safely(ds)
        }
        stop(sprintf(
          "[as.matrix,H5ParcellatedScanSummary] Failed to read summary data for scan '%s' from '%s'. Original error: %s",
          x@scan_name, dset_path, e$message
        ))
      }
    )

    if (is.null(mat_data)) {
      stop(sprintf("[as.matrix,H5ParcellatedScanSummary] Failed to retrieve summary data matrix for scan '%s'. Result is NULL.", x@scan_name))
    }

    # Optional: Validate dimensions against n_time if available?
    # if (!is.na(x@n_time) && nrow(mat_data) != x@n_time) { ... warning ... }

    # Set column names if available and dimensions match
    if (length(x@cluster_names) > 0) {
      if (length(x@cluster_names) == ncol(mat_data)) {
        colnames(mat_data) <- x@cluster_names
      } else {
        warning(sprintf(
          "[as.matrix,H5ParcellatedScanSummary] Length of cluster_names (%d) != number of columns (%d) for scan '%s'. Names not set.",
          length(x@cluster_names), ncol(mat_data), x@scan_name
        ))
      }
    }

    mat_data
  }
)

#' @rdname as.data.frame-methods
#' @export
#' @family H5Cluster
setMethod(
  f = "as.data.frame",
  signature = signature(x = "H5ParcellatedScanSummary"),
  definition = function(x, row.names = NULL, optional = FALSE, ...) {
    # Calls the as.matrix method defined above for H5ParcellatedScanSummary
    mat_data <- as.matrix(x)
    df <- as.data.frame(mat_data, row.names = row.names, optional = optional, ...)
    df
  }
)

# Override voxel-level accessors to prevent misuse


#' @export
#' @family H5Cluster
#' @rdname extract-methods
setMethod(
  "[",
  signature(x = "H5ParcellatedScanSummary", i = "ANY", j = "ANY", drop = "ANY"),
  function(x, i, j, k, l, ..., drop = TRUE) {
    stop("Voxel-level subsetting ([i,j,k,l]) is not available for H5ParcellatedScanSummary objects. Use as.matrix() or as.data.frame() to access summary data.")
  }
)


#' @export
#' @family H5Cluster
#' @rdname extract-methods
setMethod("[",
  signature(x = "H5ParcellatedScanSummary", i = "ANY", j = "missing", drop = "ANY"),
  definition = function(x, i, ..., drop = TRUE) {
    stop("Voxel-level subsetting (e.g., x[indices]) is not available for H5ParcellatedScanSummary objects. Use as.matrix() or as.data.frame() to access summary data.")
  }
)

#' @rdname series-methods
#' @export
#' @family H5Cluster
setMethod(
  f = "series",
  signature = signature(x = "H5ParcellatedScanSummary", i = "ANY"), # Catch numeric or matrix
  definition = function(x, i, ...) {
    stop("Voxel-level time series extraction (series()) is not available for H5ParcellatedScanSummary objects. Use as.matrix() or as.data.frame() to access summary data.")
  }
)


#' @rdname linear_access-methods
#' @export
#' @family H5Cluster
#' @importFrom neuroim2 linear_access
setMethod(
  f = "linear_access",
  signature = signature(x = "H5ParcellatedScanSummary", i = "numeric"),
  definition = function(x, i, ...) {
    stop("Voxel-level linear access (linear_access()) is not available for H5ParcellatedScanSummary objects. Use as.matrix() or as.data.frame() to access summary data.")
  }
)

#' Show Method for H5ParcellatedScanSummary
#' @param object An H5ParcellatedScanSummary object
#' @importFrom cli cli_h1 cli_h2 cli_li col_blue col_green col_yellow col_silver
#' @importFrom cli col_red col_magenta symbol cli_text cli_alert_info
#' @importFrom utils head
#' @importFrom methods show
#' @export
#' @family H5Cluster
#' @rdname show-methods
setMethod(
  f = "show",
  signature = "H5ParcellatedScanSummary",
  definition = function(object) {
    cli::cli_h1(cli::col_blue("H5ParcellatedScanSummary"))

    cli::cli_h2(cli::col_yellow("Run Summary Info"))
    cli::cli_li(items = c(
      "{cli::col_green(\'Scan Name\')} : {object@scan_name}",
      "{cli::col_green(\'Time points\')} : {object@n_time}",
      "{cli::col_green(\'Summary Dset\')} : {object@summary_dset}"
    ))

    # Display cluster names/IDs (show first few)
    n_clus_names <- length(object@cluster_names)
    n_clus_ids <- length(object@cluster_ids)
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


    cli::cli_h2(cli::col_yellow("Shared Info (from H5ParcellatedArray)"))
    # Check if clusters slot is valid and has data
    has_clusters <- !is.null(object@clusters) && inherits(object@clusters, "ClusteredNeuroVol") && length(object@clusters@clusters) > 0
    n_clusters_from_slot <- if (has_clusters) length(unique(object@clusters@clusters)) else NA

    cli::cli_li(items = c(
      "{cli::col_green(\'Active voxels in mask\')} : {object@n_voxels}",
      "{cli::col_green(\'Number of clusters (map)\')} : {ifelse(has_clusters, n_clusters_from_slot, \'(NA or Map not loaded)\')}"
    ))

    cli::cli_h2(cli::col_yellow("Storage: HDF5 File"))
    h5_valid <- !is.null(object@obj) && inherits(object@obj, "H5File") && object@obj$is_valid
    if (h5_valid) {
      filename <- tryCatch(object@obj$get_filename(), error = function(e) "Error getting name")
      summary_path <- sprintf(H5_PATHS$SUMMARY_DSET_TPL, object@scan_name, object@summary_dset)
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
