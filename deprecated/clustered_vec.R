#' @include all_class.R
#' @importMethodsFrom neuroim2 series
NULL

# =========================================================================
# DEPRECATED FILE
# -------------------------------------------------------------------------
# The contents of this file are part of a deprecated framework
# (H5ClusteredVec, H5ReducedClusteredVec, H5ClusteredVecSeq).
# These classes and their methods are being replaced by the
# H5ClusteredExperiment framework and its associated run objects
# (H5ClusteredRunFull, H5ClusteredRunSummary).
#
# This file and its contents will be removed in a future version.
# Please update your code to use H5ClusteredExperiment.
# =========================================================================

#' Create an H5ClusteredVec object (DEPRECATED)
#'
#' @description
#' This constructor is deprecated. Use `make_run_full()` instead to create
#' an `H5ClusteredRunFull` object.
#' Constructs a single-scan \code{H5ClusteredVec} object referencing a clustered
#' time-series dataset stored in an HDF5 file.
#'
#' @param obj An \code{\link[hdf5r]{H5File}} object referencing the open HDF5 file.
#' @param scan_name \code{character}, the scan identifier (subgroup under \code{/scans/}).
#' @param mask A \code{\link[neuroim2]{LogicalNeuroVol}} representing the 3D brain mask.
#' @param clusters A \code{\link[neuroim2]{ClusteredNeuroVol}} describing cluster assignments.
#' @param n_time The number of time points (integer) for this scan.
#'
#' @return A new \code{H5ClusteredRunFull} instance (via `make_run_full`).
#' @importFrom methods new
#' @importFrom lifecycle deprecate_warn
#' @export
H5ClusteredVec <- function(obj, scan_name, mask, clusters, n_time) {

  # Issue deprecation warning
  lifecycle::deprecate_warn(
    when = "0.1.0", # Replace with the version number where deprecation happens
    what = "H5ClusteredVec()",
    with = "make_run_full()",
    details = "The H5ClusteredVec class is being replaced by H5ClusteredRunFull and H5ClusteredRunSummary within the new H5ClusteredExperiment framework."
  )
  
  # Call the new constructor
  # Need to handle potential errors from make_run_full if called directly
  # make_run_full might close the file handle if it opened it, which is different from H5ClusteredVec's behavior
  # Let's call it but note the potential difference in handle management
  make_run_full(file_source = obj,
                  scan_name = scan_name,
                  mask = mask,
                  clusters = clusters,
                  n_time = n_time)
  
  # Original code is removed/commented out as we now delegate
  # # Validate input types and values
  # ... (rest of original validation and new("H5ClusteredVec", ...) call)
}







#' --- Helper: Get cluster timeseries by mask index ---
#' This function replaces the non-compliant linear_access method for H5ClusteredVec.
#' It adheres to the specific data access pattern needed for the HDF5 cluster structure.
.get_cluster_timeseries_by_mask_index <- function(x, mask_indices, time_indices=NULL) {
    # 1. Validate spatial indices \'mask_indices\' (relative to mask)
    nVoxMask <- x@n_voxels # Use stored voxel count
    if (any(mask_indices < 1 | mask_indices > nVoxMask)) {
      stop(sprintf("[.get_cluster_timeseries] Some indices in 'mask_indices' are out of the valid mask range [1..%d]", nVoxMask))
    }
    mask_indices <- as.integer(mask_indices)
    n_request <- length(mask_indices)
    full_time_length <- x@n_time # Use stored time length

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
    clus_ids_req <- tryCatch(x@clusters@clusters[mask_indices], error = function(e) {
        stop(sprintf("[.get_cluster_timeseries] Failed to get cluster assignments for provided indices. Error: %s", e$message))
    })

    # 5. Group requested indices by their cluster ID
    # Map from original index in `mask_indices` to row in result_mat is simply 1:n_request
    index_groups <- split(seq_len(n_request), clus_ids_req)
    needed_cluster_ids <- names(index_groups)
    base_path <- paste0("/scans/", x@scan_name, "/clusters/")

    # 6. Loop through needed clusters and read data
    for (cid_str in needed_cluster_ids) {
        cid <- as.integer(cid_str) # Cluster ID should be integer
        # Rows in the *result matrix* corresponding to this cluster
        row_indices_in_result <- index_groups[[cid_str]]
        # Original mask indices requested that belong to this cluster
        mask_indices_this_cluster_req <- mask_indices[row_indices_in_result]

        ds <- NULL # Initialize ds to NULL
        dset_path <- paste0(base_path, "cluster_", cid) # sprintf might be slightly slower here

        tryCatch({
            # --- Optimized Cluster Processing ---
            # Get the full list of mask indices for this cluster ONCE
            all_mask_indices_in_this_cluster <- which(x@clusters@clusters == cid)

            # Map requested mask indices to row offsets within this cluster\'s dataset using match
            row_offsets_in_dataset <- match(mask_indices_this_cluster_req, all_mask_indices_in_this_cluster)
            # --- End Optimization ---

            if(any(is.na(row_offsets_in_dataset))) {
                stop(sprintf("Internal inconsistency: some requested voxels mapped to cluster %d but not found in its index list.", cid))
            }
            row_offsets_in_dataset <- as.integer(row_offsets_in_dataset) # Ensure integer

            # Check dataset existence
            if (!x@obj$exists(dset_path)){
                stop(sprintf("Dataset for cluster %d not found at path: %s", cid, dset_path))
            }
            ds <- x@obj[[dset_path]] # Get dataset handle
            on.exit(if (!is.null(ds) && inherits(ds, "H5D") && ds$is_valid) ds$close(), add = TRUE, after = FALSE)

            # Read the required rows and SPECIFIED time columns (or all if time_indices was NULL)
            # Use drop = FALSE to ensure result is always a matrix
            cluster_data_subset <- ds[row_offsets_in_dataset, time_indices, drop = FALSE]

            # Dimension check before assignment
            expected_rows <- length(row_indices_in_result)
            if (!is.matrix(cluster_data_subset) ||
                 nrow(cluster_data_subset) != expected_rows ||
                 ncol(cluster_data_subset) != n_time_request) {
                 stop(sprintf("Read data dimensions [%s] mismatch expected [%d,%d] for cluster %d",
                              paste(dim(cluster_data_subset), collapse=","), expected_rows, n_time_request, cid))
            }

            # Assign to result matrix
            result_mat[row_indices_in_result, ] <- cluster_data_subset

        }, error = function(e) {
            # Ensure ds handle is attempted to be closed even on error
            if (!is.null(ds) && inherits(ds, "H5D") && ds$is_valid) {
                try(ds$close(), silent = TRUE)
            }
            stop(sprintf("[.get_cluster_timeseries] Failed processing cluster %d. Original error: %s", cid, e$message))
        })
    } # End loop over clusters

    # 7. Final check for NAs (should not happen if all clusters read successfully)
    if (any(is.na(result_mat))) {
        warning("[.get_cluster_timeseries] Result matrix contains NAs. Data reading may have failed for some clusters despite no error being thrown.")
    }

    return(result_mat)
  }

### --- helper shared by every signature ---------------------------------------
.subset_h5cvec <- function(x, i, j, k, l, drop = TRUE) {

  dims   <- dim(x@mask)
  nvox   <- x@n_voxels
  nt     <- x@n_time

  # ---------- 1.  default / normalise inputs -----------------------
  # treat "missing *or* NULL" uniformly
  if (missing(j) || is.null(j)) j <- NULL
  if (missing(k) || is.null(k)) k <- NULL
  if (missing(l) || is.null(l)) l <- seq_len(nt)

  # ---------- 2.  fast path : mask-based indexing ------------------
  if (is.null(j) && is.null(k)) {                # only one spatial arg
    if (is.numeric(i) && all(i %in% seq_len(nvox)))
    {
      # treat i as mask indices
      # message(sprintf("\\n>>> Debug .subset_h5cvec: FAST PATH calling .get_cluster_timeseries with mask_indices = [%s]\\n", paste(i, collapse=", "))) # DEBUG removed
      dat <- .get_cluster_timeseries_by_mask_index(x, mask_indices = as.integer(i), time_indices = l)  # [nVox, nT]
      if (drop) return(drop(dat))
      return(dat)
    }
    # … otherwise fall through and interpret i as x-coordinate
    # message(sprintf("\\n>>> Debug .subset_h5cvec: FAST PATH FALLTHROUGH for i = [%s]\\n", paste(i, collapse=", "))) # DEBUG removed
  }

  # ---------- 3.  coordinate path ---------------------------------
  if (missing(i) || is.null(i)) i <- seq_len(dims[1]) # Default i if missing
  if (is.null(j)) j <- seq_len(dims[2])
  if (is.null(k)) k <- seq_len(dims[3])


  if (!is.numeric(i) || !is.numeric(j) || !is.numeric(k) || !is.numeric(l)) {
      stop("Subscripts must be numeric or NULL.")
  }

  if (any(i < 1L | i > dims[1L]) ||
      any(j < 1L | j > dims[2L]) ||
      any(k < 1L | k > dims[3L]))
      stop("spatial subscript out of bounds")

  if (any(l < 1L | l > nt))
      stop("time subscript out of bounds")

  # build 4-D result filled with 0
  out <- array(0,
               dim = c(length(i), length(j), length(k), length(l)))

  # map requested (i,j,k) to linear positions in whole volume
  grid <- as.matrix(expand.grid(i = i, j = j, k = k, KEEP.OUT.ATTRS = FALSE))
  lin  <- grid[,1] + (grid[,2]-1L)*dims[1] + (grid[,3]-1L)*dims[1]*dims[2]

  mask <- as.logical(as.array(x@mask))
  inside <- mask[lin]
  if (!any(inside)) return(if (drop) drop(out) else out)

  # translate to 1..n_voxels inside mask
  mask_lin <- which(mask)
  v_idx    <- match(lin[inside], mask_lin)

  # message(sprintf("\\n>>> Debug .subset_h5cvec: COORD PATH calling .get_cluster_timeseries with v_idx = [%s]\\n", paste(v_idx, collapse=", "))) # DEBUG removed
  dat  <- .get_cluster_timeseries_by_mask_index(x, mask_indices = v_idx, time_indices = l)   # [nInside, nT]

  # scatter back
  idx <- grid[inside,, drop = FALSE]
  out_i <- match(idx[,1], i)
  out_j <- match(idx[,2], j)
  out_k <- match(idx[,3], k)

  for (t in seq_along(l))
    out[cbind(as.integer(out_i), as.integer(out_j), as.integer(out_k), t)] <- dat[, t]

  if (drop) drop(out) else out
}

### --- universal [ method -----------------------------------------------------
setMethod("[",
  signature(x = "H5ClusteredVec", i = "ANY", j = "ANY", drop = "ANY"),
  function(x, i, j, k, l, ..., drop = TRUE)
      .subset_h5cvec(x, i, j, k, l, drop))

### --- single-argument [ method ------------------------------------------------
setMethod("[",
  signature(x = "H5ClusteredVec", i = "ANY", j = "missing"),
  definition = function(x, i, ..., drop = TRUE) {
      .subset_h5cvec(x, i = i, j = NULL, k = NULL, l = NULL, drop = drop)
  })

#' Get dimensions of an object
#' @param x the object
#' @return A numeric vector of dims
#' @export
#' @rdname dim-methods
setMethod("dim", "H5ClusteredVec",
  function(x) {
    n_time <- x@n_time # Should be integer and positive from constructor
    # Add validation check here
    if (is.null(n_time) || is.na(n_time) || length(n_time) != 1 || n_time <= 0 || floor(n_time) != n_time) {
      stop(sprintf("[dim,H5ClusteredVec] Time dimension (%s) is invalid (must be a single positive integer).", as.character(n_time)))
    }
    c(dim(x@mask), n_time)
  }
)

#' Extract Time Series from H5ClusteredVec
#'
#' @description
#' Extracts time series data for specified voxels from an \code{H5ClusteredVec} object.
#' Voxels can be specified using mask-based linear indices, a 3-column coordinate matrix, or individual x,y,z coordinates.
#'
#' @param x An \code{H5ClusteredVec} object.
#' @param i Either a numeric vector of mask-based indices (1 to sum(mask)),
#'   a numeric matrix with 3 columns (x, y, z coordinates),
#'   or the x-coordinate if j and k are also provided.
#' @param j Optional y-coordinate (if i is x-coordinate).
#' @param k Optional z-coordinate (if i is x-coordinate).
#' @param ... Not used.
#'
#' @return A matrix where rows are time points and columns are the requested voxels \code{[nTime, nVoxels]}.
#' @export
#' @rdname series-methods
setMethod(
  f = "series",
  signature = signature(x = "H5ClusteredVec"),
  definition = function(x, i, j, k, ...) {

    dims_mask <- dim(x@mask)
    n_vox_mask <- x@n_voxels
    mask_indices_req <- NULL

    # Determine how voxels are specified
    if (!missing(j) && !missing(k)) {
      # Case 1: i, j, k are individual coordinates
      coords <- matrix(c(i, j, k), nrow = length(i), ncol = 3) # Allow vector i,j,k
    } else if (is.matrix(i) && ncol(i) == 3) {
      # Case 2: i is a coordinate matrix
      coords <- i
    } else if (is.numeric(i) && (missing(j) || is.null(j)) && (missing(k) || is.null(k))) {
       # Case 3: i is assumed to be mask-based linear indices (and j, k are missing/NULL)
       coords <- NULL
       mask_indices_req <- as.integer(i)
       # Validate mask indices
       if (any(mask_indices_req < 1 | mask_indices_req > n_vox_mask)) {
           stop(sprintf("[series,H5ClusteredVec] Mask-based indices are out of range [1..%d]", n_vox_mask))
       }
    } else {
      stop("[series,H5ClusteredVec] Invalid specification. Provide mask indices (i), a 3-col coordinate matrix (i), or coordinate vectors (i, j, k).")
    }


    # If coordinates were provided, convert them to mask indices
    if (!is.null(coords)) {
        # Validate coordinate bounds
        if (any(coords[,1] < 1 | coords[,1] > dims_mask[1]) ||
            any(coords[,2] < 1 | coords[,2] > dims_mask[2]) ||
            any(coords[,3] < 1 | coords[,3] > dims_mask[3])) {
            stop("[series,H5ClusteredVec] Coordinates are out of bounds.")
        }

        # Convert coords to linear indices in the full 3D grid
        linear_idx_full <- coords[,1] + (coords[,2]-1)*dims_mask[1] + (coords[,3]-1)*dims_mask[1]*dims_mask[2]

        # Find which of these linear indices are actually in the mask
        mask_array <- as.logical(as.array(x@mask))
        in_mask_subset <- mask_array[linear_idx_full]

        if (!all(in_mask_subset)) {
            warning("[series,H5ClusteredVec] Some requested coordinates fall outside the mask and will be ignored.")
            linear_idx_full <- linear_idx_full[in_mask_subset]
            if (length(linear_idx_full) == 0) {
                return(matrix(numeric(0), nrow = x@n_time, ncol = 0)) # Return empty matrix if no valid voxels
            }
        }

        # Get the corresponding indices *within the mask* (1..n_voxels)
        global_mask_indices <- which(mask_array)
        mask_indices_req <- match(linear_idx_full, global_mask_indices)

        if (any(is.na(mask_indices_req))) {
             stop("[series,H5ClusteredVec] Internal error: Failed to map some valid coordinates to mask indices.")
        }
    }

    # Ensure we have valid indices before proceeding
    if (is.null(mask_indices_req) || length(mask_indices_req) == 0) {
       return(matrix(numeric(0), nrow = x@n_time, ncol = 0))
    }

    # Call the RENAMED helper function to get data [nVoxels, nTime]
    voxel_data <- .get_cluster_timeseries_by_mask_index(x, mask_indices = mask_indices_req, time_indices = NULL)

    # Transpose to return [nTime, nVoxels]
    return(t(voxel_data))

  }
)

#' Pretty Printer for H5ClusteredVec (Revised)
#' @importFrom crayon bold blue silver yellow green magenta red
#' @importFrom methods show
#' @export
setMethod(
  f = "show",
  signature = "H5ClusteredVec",
  definition = function(object) {
    cat("\n", crayon::bold(crayon::blue("H5ClusteredVec")), "\n", sep = "")
    cat(crayon::silver("────────────────────────────────────────\n"))
    cat(crayon::bold(crayon::yellow("Basic Info")), "\n")
    cat(crayon::silver(" • "), crayon::green("Scan Name:"), object@scan_name, "\n")
    # Use n_voxels slot
    cat(crayon::silver(" • "), crayon::green("Active voxels in mask:"), object@n_voxels, "\n")
    cluster_ids <- unique(object@clusters@clusters)
    n_clusters  <- length(cluster_ids)
    cat(crayon::silver(" • "), crayon::green("Number of clusters:"), n_clusters, "\n")
    # Use n_time slot
    cat(crayon::silver(" • "), crayon::green("Time points:        "), object@n_time, "\n")
    cat(crayon::bold("\nStorage:"), "\n")
    if (inherits(object@obj, "H5File") && object@obj$is_valid) {
      cat(crayon::silver(" • "), "HDF5 file: ",
          crayon::magenta(object@obj$get_filename()), "\n", sep="")
    } else {
      cat(crayon::silver(" • "), "HDF5 file is ",
          crayon::red("INVALID or CLOSED"), "\n", sep="")
    }
    cat("\n")
  }
)




#' @rdname linear_access-methods
#' @importFrom neuroim2 NeuroSpace space
#' @export
setMethod(
  f = "linear_access",
  signature = signature(x="H5ClusteredVec", i="numeric"),
  definition = function(x, i, ...) {

    full_dims <- dim(x) # Get 4D dimensions
    n_request <- length(i)
    if (n_request == 0) return(numeric(0))

    # Validate 4D linear indices
    max_elements <- prod(full_dims)
    if (any(i < 1 | i > max_elements)) {
        stop(sprintf("[linear_access,H5ClusteredVec] 4D linear indices are out of range [1..%d]", max_elements))
    }

    # Convert 4D linear indices to 4D coordinates
    coords4d <- arrayInd(i, .dim = full_dims)

    # Prepare mask and related info
    mask_array <- as.logical(as.array(x@mask))
    dims_mask <- dim(mask_array) # 3D dims
    global_mask_indices <- which(mask_array) # Mask indices in 3D space

    # Calculate 3D linear indices from the first 3 columns of coords4d
    coords3d <- coords4d[, 1:3, drop = FALSE]
    lin3d <- coords3d[,1] + (coords3d[,2]-1L)*dims_mask[1] + (coords3d[,3]-1L)*dims_mask[1]*dims_mask[2]

    # Identify which requested indices fall within the mask
    in_mask_logical <- mask_array[lin3d]

    # Initialize result vector (use NA to distinguish missing from zero)
    result_vec <- rep(NA_real_, n_request)

    # Process only the indices that are inside the mask
    if (any(in_mask_logical)) {
        valid_indices_in_i <- which(in_mask_logical)
        valid_coords4d <- coords4d[valid_indices_in_i, , drop = FALSE]
        valid_lin3d <- lin3d[valid_indices_in_i]

        # Calculate mask-relative indices (1..n_voxels)
        mask_indices_req <- match(valid_lin3d, global_mask_indices)

        # Group requests by time coordinate (column 4)
        time_coords <- valid_coords4d[, 4]
        requests_by_time <- split(data.frame(orig_idx = valid_indices_in_i, mask_idx = mask_indices_req), time_coords)

        # Process each required time point
        for (t_val_str in names(requests_by_time)) {
            t_val <- as.integer(t_val_str)
            current_requests <- requests_by_time[[t_val_str]]
            mask_indices_for_t <- current_requests$mask_idx
            original_indices_for_t <- current_requests$orig_idx

            # Fetch data for these voxels at this time point using the helper
            # Helper returns [length(mask_indices_for_t), 1] matrix
            data_subset <- .get_cluster_timeseries_by_mask_index(x, mask_indices = mask_indices_for_t, time_indices = t_val)

            # Assign fetched values to the correct positions in the result vector
            result_vec[original_indices_for_t] <- data_subset[, 1]
        }
    }

    # Replace remaining NAs (outside mask) with 0, consistent with sparse interpretation
    result_vec[is.na(result_vec)] <- 0

    return(result_vec)
  }
)

# --- series methods for H5ClusteredVec ---

#' @description
#' Extracts time series data for specified voxels from an \code{H5ClusteredVec} object.
#' Voxels can be specified either using mask-based linear indices (numeric vector)
#' or a 3-column coordinate matrix.
#'
#' @param x An \code{H5ClusteredVec} object.
#' @param i Either a numeric vector of mask-based indices (1 to sum(mask)),
#'   or a numeric matrix with 3 columns (x, y, z coordinates).
#' @param ... Not used.
#'
#' @return A matrix where rows are time points and columns are the requested voxels \code{[nTime, nVoxels]}.
#' @rdname series-methods


# Method 1: i = numeric (handles mask indices OR single i,j,k coordinates)
#' @export
setMethod(
  f = "series",
  signature = signature(x = "H5ClusteredVec", i = "numeric"),
  definition = function(x, i, j, k, ...) { # Define with j, k to check missing status

    mask_indices_req <- NULL

    if (missing(j) && missing(k)) {
        # --- Case 1: i is numeric mask indices --- 
        n_vox_mask <- x@n_voxels
        mask_indices_req <- as.integer(i)

        # Validate mask indices
        if (any(mask_indices_req < 1 | mask_indices_req > n_vox_mask)) {
            stop(sprintf("[series,H5ClusteredVec] Mask-based indices are out of range [1..%d]", n_vox_mask))
        }

    } else if (!missing(j) && !missing(k)) {
        # --- Case 2: i, j, k are single coordinates ---
        if (length(i) != 1 || length(j) != 1 || length(k) != 1) {
             stop("[series,H5ClusteredVec] When providing i, j, k, they must each be single numeric values.")
        }
        coords <- matrix(c(i, j, k), nrow = 1, ncol = 3)
        
        # Convert single coordinate to mask index (reuse logic from matrix method)
        dims_mask <- dim(x@mask)
        mask_array <- as.logical(as.array(x@mask))
        global_mask_indices <- which(mask_array)

        # Validate coordinate bounds
        if (any(coords[,1] < 1 | coords[,1] > dims_mask[1]) ||
            any(coords[,2] < 1 | coords[,2] > dims_mask[2]) ||
            any(coords[,3] < 1 | coords[,3] > dims_mask[3])) {
            stop("[series,H5ClusteredVec] Coordinates are out of bounds.")
        }

        linear_idx_full <- coords[,1] + (coords[,2]-1)*dims_mask[1] + (coords[,3]-1)*dims_mask[1]*dims_mask[2]

        if (!mask_array[linear_idx_full]) {
             warning("[series,H5ClusteredVec] Requested coordinate falls outside the mask.")
             return(matrix(numeric(0), nrow = x@n_time, ncol = 0)) # Return empty if outside
        }
        
        mask_indices_req <- match(linear_idx_full, global_mask_indices)
        if (is.na(mask_indices_req)) {
             stop("[series,H5ClusteredVec] Internal error: Failed to map valid coordinate to mask index.")
        }
        
    } else {
        # Invalid combination (e.g., only i and j provided)
        stop("[series,H5ClusteredVec] Invalid arguments. Provide numeric mask indices (i), a 3-column matrix (i), or single coordinates (i, j, k).")
    }

    # --- Common section: Fetch data using mask indices --- 
    if (is.null(mask_indices_req) || length(mask_indices_req) == 0) {
       return(matrix(numeric(0), nrow = x@n_time, ncol = 0))
    }

    voxel_data <- .get_cluster_timeseries_by_mask_index(x, mask_indices = mask_indices_req, time_indices = NULL)
    return(t(voxel_data)) # Transpose to [nTime, nVoxels]
  }
)

# Method 2: i = matrix (coordinate matrix)
#' @export
setMethod(
  f = "series",
  signature = signature(x = "H5ClusteredVec", i = "matrix"),
  definition = function(x, i, ...) { 
    # Removed check for j, k here as numeric method handles i,j,k case
    
    coords <- i
    if (ncol(coords) != 3) {
        stop("[series,H5ClusteredVec] Coordinate matrix 'i' must have 3 columns.")
    }

    dims_mask <- dim(x@mask)
    mask_array <- as.logical(as.array(x@mask))
    global_mask_indices <- which(mask_array)

    # Validate coordinate bounds
    if (any(coords[,1] < 1 | coords[,1] > dims_mask[1]) ||
        any(coords[,2] < 1 | coords[,2] > dims_mask[2]) ||
        any(coords[,3] < 1 | coords[,3] > dims_mask[3])) {
        stop("[series,H5ClusteredVec] Coordinates are out of bounds.")
    }

    # Convert coords to linear indices in the full 3D grid
    linear_idx_full <- coords[,1] + (coords[,2]-1)*dims_mask[1] + (coords[,3]-1)*dims_mask[1]*dims_mask[2]

    # Find which of these linear indices are actually in the mask
    in_mask_subset <- mask_array[linear_idx_full]
    original_row_indices <- seq_len(nrow(coords))

    if (!all(in_mask_subset)) {
        warning("[series,H5ClusteredVec] Some requested coordinates fall outside the mask and will be ignored.")
        valid_original_indices <- original_row_indices[in_mask_subset]
        linear_idx_full <- linear_idx_full[in_mask_subset]
        coords <- coords[in_mask_subset, , drop = FALSE] 
        if (length(linear_idx_full) == 0) {
            return(matrix(numeric(0), nrow = x@n_time, ncol = 0))
        }
    } else {
        valid_original_indices <- original_row_indices
    }

    # Get the corresponding indices *within the mask* (1..n_voxels)
    mask_indices_req <- match(linear_idx_full, global_mask_indices)

    if (any(is.na(mask_indices_req))) {
         stop("[series,H5ClusteredVec] Internal error: Failed to map some valid coordinates to mask indices.")
    }

    if (length(mask_indices_req) == 0) {
       return(matrix(numeric(0), nrow = x@n_time, ncol = 0))
    }

    # Call the internal helper function
    voxel_data <- .get_cluster_timeseries_by_mask_index(x, mask_indices = mask_indices_req, time_indices = NULL)
    
    # Transpose to [nTime, nReturnedVoxels]
    voxel_data_transposed <- t(voxel_data)
    
    # Ensure column names match original valid coordinates if possible
    if (nrow(coords) == ncol(voxel_data)) { # Check if dimensions match after helper call
       colnames(voxel_data_transposed) <- apply(coords, 1, paste, collapse=",")
    }
    return(voxel_data_transposed)
  }
)
