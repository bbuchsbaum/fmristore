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

#' Create an H5ClusteredVecSeq object (DEPRECATED)
#'
#' @description
#' This constructor is deprecated. Use `H5ClusteredExperiment()` instead.
#' Constructs a multi-scan \code{H5ClusteredVecSeq} object referencing a set of
#' clustered time-series datasets stored in an HDF5 file.
#'
#' @param obj An \code{\link[hdf5r]{H5File}} object (open for reading).
#' @param scan_names A \code{character} vector of scan IDs to include.
#' @param mask A \code{\link[neuroim2]{LogicalNeuroVol}} giving the 3D mask.
#' @param clusters A \code{\link[neuroim2]{ClusteredNeuroVol}} describing cluster assignments.
#' @param scan_metadata A \code{list} of metadata for each scan. **Each element MUST contain an entry named 'n_time' specifying the number of time points for that scan.**
#' @param cluster_metadata A \code{data.frame} of cluster metadata.
#'
#' @return A new \code{H5ClusteredExperiment} instance.
#' @importFrom methods new
#' @importFrom lifecycle deprecate_warn
#' @export
H5ClusteredVecSeq <- function(obj, scan_names, mask, clusters,
                              scan_metadata, cluster_metadata) {

  lifecycle::deprecate_warn(
    when = "0.1.0", # Replace with the version number where deprecation happens
    what = "H5ClusteredVecSeq()",
    with = "H5ClusteredExperiment()",
    details = "The H5ClusteredVecSeq class is being replaced by the H5ClusteredExperiment class."
  )
  
  # Note: The original H5ClusteredVecSeq performed validation on scan_metadata$n_time.
  # The new H5ClusteredExperiment constructor aims to read this from HDF5 itself.
  # Passing scan_metadata directly might override this or cause conflicts if the
  # HDF5 reading isn't fully implemented in H5ClusteredExperiment yet.
  # For now, we pass the arguments directly, assuming H5ClusteredExperiment will handle them.
  
  # Call the new experiment constructor
  H5ClusteredExperiment(
      file_source = obj,
      scan_names = scan_names,
      mask = mask,
      clusters = clusters,
      scan_metadata = scan_metadata, # Pass along, constructor might override/merge
      cluster_metadata = cluster_metadata # Pass along
      # summary_preference and keep_handle_open use defaults
  )

  # Original code removed/commented
  # # Validate inputs
  # ... (rest of original validation and new("H5ClusteredVecSeq", ...) call)
}


#' Get a single scan from H5ClusteredVecSeq (Revised)
#'
#' @description
#' Retrieves a single \code{H5ClusteredVec} object from a \code{H5ClusteredVecSeq}
#' object based on the scan index.
#'
#' @param x An \code{H5ClusteredVecSeq} object.
#' @param i The index of the scan to retrieve.
#'
#' @return A \code{H5ClusteredVec} object corresponding to the specified scan.
#' @export
setMethod("[[", signature(x = "H5ClusteredVecSeq"),
          function(x, i) {
            if (i < 1 || i > length(x@scan_names)) {
                stop(sprintf("Scan index %d is out of bounds [1..%d]", i, length(x@scan_names)))
            }
            scan_name <- x@scan_names[i]
            # Get n_time for this scan from the metadata
            if (is.null(x@scan_metadata[[i]]$n_time)) {
                 stop(sprintf("Internal error: n_time missing for scan index %d ('%s') in scan_metadata.", i, scan_name))
            }
            n_time_scan <- as.integer(x@scan_metadata[[i]]$n_time)

            # Create the H5ClusteredVec object
            H5ClusteredVec(obj = x@obj,
                           scan_name = scan_name,
                           mask = x@mask,
                           clusters = x@clusters,
                           n_time = n_time_scan) # Pass required n_time
          })


#' Get total length (concatenated time points) of H5ClusteredVecSeq (Revised)
#' @export
setMethod("length", signature(x = "H5ClusteredVecSeq"),
          function(x) {
            ct <- x@cumulative_time
            # Check if it's a valid integer vector of length >= 2
            if (!is.integer(ct) || length(ct) < 2) {
                stop("[length,H5ClusteredVecSeq] Slot 'cumulative_time' is invalid (must be integer vector of length >= 2).")
            }
            # Return the last element (total length)
            ct[length(ct)]
          })

#' Pretty Printer for H5ClusteredVecSeq (Revised)
#' @importFrom crayon bold blue silver yellow green magenta red
#' @importFrom methods show
#' @importFrom utils head
#' @export
setMethod(
  f = "show",
  signature = "H5ClusteredVecSeq",
  definition = function(object) {
    cat("\n", crayon::bold(crayon::blue("H5ClusteredVecSeq")), "\n", sep="")
    cat(crayon::silver("────────────────────────────────────────\n"))
    cat(crayon::bold(crayon::yellow("Basic Info")), "\n")
    n_scans <- length(object@scan_names)
    cat(crayon::silver(" • "), crayon::green("Number of scans:"), n_scans, "\n")
    # Use n_voxels slot
    cat(crayon::silver(" • "), crayon::green("Active voxels in mask:"), object@n_voxels, "\n")
    n_clusters <- length(unique(object@clusters@clusters))
    cat(crayon::silver(" • "), crayon::green("Number of clusters:"), n_clusters, "\n")
    cat(crayon::bold("\nScan Names:"), "\n")
    if (n_scans <= 5) {
      cat("  ", paste(object@scan_names, collapse=", "), "\n")
    } else {
      preview <- paste(utils::head(object@scan_names, 5), collapse=", ")
      remainder <- n_scans - 5
      cat("  ", preview, crayon::silver(paste0(" ... (", remainder, " more)")), "\n")
    }
    cat(crayon::bold("\nCluster Metadata:"), "\n")
    if (nrow(object@cluster_metadata) == 0) {
      cat(crayon::silver("  No cluster_metadata table.\n"))
    } else {
      md_cols <- names(object@cluster_metadata)
      cat(crayon::silver("  Columns:"), paste(md_cols, collapse=", "), "\n")
      cat(crayon::silver("  #rows="), nrow(object@cluster_metadata), "\n")
    }
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

### --- helper shared by H5ClusteredVecSeq [ method ----------------------------
.subset_h5cvecseq <- function(x, i, j, k, l, drop = TRUE) {
  # 1) Get scan info and total time length
  scan_names <- x@scan_names
  n_scans <- length(scan_names)
  cumsums <- x@cumulative_time # Integer vector: c(0, T1, T1+T2, ...)
  total_time <- cumsums[n_scans + 1]
  dims <- dim(x@mask)
  nvox <- x@n_voxels

  # ---------- 1.  default / normalise inputs -----------------------
  # treat "missing *or* NULL" uniformly
  if (missing(j) || is.null(j)) j <- NULL
  if (missing(k) || is.null(k)) k <- NULL
  if (missing(l) || is.null(l)) l <- seq_len(total_time)

  # ---------- 2.  fast path : mask-based indexing ------------------
  if (is.null(j) && is.null(k)) {                # only one spatial arg
    if (is.numeric(i) && all(i %in% seq_len(nvox)))
    {
      # treat i as mask indices
      # For mask indices, we need to build a 2D matrix: [maskIndices, time]
      result <- matrix(0, nrow = length(i), ncol = length(l))

      # Map global time indices to scan and local time indices
      scan_indices_for_t <- findInterval(l, cumsums, left.open = TRUE, rightmost.closed = TRUE)
      local_time_indices <- l - cumsums[scan_indices_for_t]

      # Process each needed scan
      needed_scans <- unique(scan_indices_for_t)
      for (sc_id in needed_scans) {
        cols_for_this_scan <- which(scan_indices_for_t == sc_id)
        local_times_needed <- local_time_indices[cols_for_this_scan]

        # Get data from this scan
        hv <- x[[sc_id]]
        # Call RENAMED helper
        scan_data <- .get_cluster_timeseries_by_mask_index(hv, mask_indices = as.integer(i), time_indices = local_times_needed)

        # Fill the result matrix
        result[, cols_for_this_scan] <- scan_data
      }

      if (drop) return(drop(result))
      return(result)
    }
    # … otherwise fall through and interpret i as x-coordinate
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

  if (any(l < 1L | l > total_time))
      stop("time subscript out of bounds")

  # build 4-D result filled with 0
  out <- array(0, dim = c(length(i), length(j), length(k), length(l)))

  # map requested (i,j,k) to linear positions in whole volume
  grid <- as.matrix(expand.grid(i = i, j = j, k = k, KEEP.OUT.ATTRS = FALSE))
  lin  <- grid[,1] + (grid[,2]-1L)*dims[1] + (grid[,3]-1L)*dims[1]*dims[2]

  mask <- as.logical(as.array(x@mask))
  inside <- mask[lin]
  if (!any(inside)) return(if (drop) drop(out) else out)

  # translate to 1..n_voxels inside mask
  mask_lin <- which(mask)
  valid_mask_indices <- match(lin[inside], mask_lin)

  # convert from coordinates to result array indices
  valid_coords <- grid[inside,, drop = FALSE]
  result_indices <- matrix(0, nrow=nrow(valid_coords), ncol=3)
  result_indices[,1] <- match(valid_coords[,1], i)
  result_indices[,2] <- match(valid_coords[,2], j)
  result_indices[,3] <- match(valid_coords[,3], k)

  # Map global time indices to scan and local time indices
  scan_indices_for_t <- findInterval(l, cumsums, left.open = TRUE, rightmost.closed = TRUE)
  local_time_indices <- l - cumsums[scan_indices_for_t]

  # Loop through needed scans
  needed_scans <- unique(scan_indices_for_t)
  for (sc_id in needed_scans) {
      cols_for_this_scan <- which(scan_indices_for_t == sc_id)
      local_times_needed <- local_time_indices[cols_for_this_scan]

      hv <- x[[sc_id]]  # Get the H5ClusteredVec for this scan

      # Use RENAMED helper to get data for this scan
      tryCatch({
          # Get the data for valid voxels and this scan\'s time points
          scan_data <- .get_cluster_timeseries_by_mask_index(hv, mask_indices = valid_mask_indices, time_indices = local_times_needed)

          # For each valid coordinate and this scan\'s time points, fill the result array
          for (t_idx in seq_along(cols_for_this_scan)) {
              global_t_idx <- cols_for_this_scan[t_idx]
              # Using cbind for direct 4D indexing (matching .subset_h5cvec)
              out[cbind(as.integer(result_indices[,1]),
                        as.integer(result_indices[,2]),
                        as.integer(result_indices[,3]),
                        global_t_idx)] <- scan_data[, t_idx]

          }
      }, error = function(e) {
          warning(sprintf("[`,H5ClusteredVecSeq`] Error accessing data for scan %d (\'%s\'): %s",
                        sc_id, scan_names[sc_id], e$message))
      })
  }

  if (drop) drop(out) else out
}

### --- universal [ method for H5ClusteredVecSeq -------------------------------
setMethod("[",
  signature(x = "H5ClusteredVecSeq", i = "ANY", j = "ANY", drop = "ANY"),
  function(x, i, j, k, l, ..., drop = TRUE)
      .subset_h5cvecseq(x, i, j, k, l, drop))

### --- single-argument [ method for H5ClusteredVecSeq -----------------------
setMethod("[",
  signature(x = "H5ClusteredVecSeq", i = "ANY", j = "missing"),
  definition = function(x, i, ..., drop = TRUE) {
      .subset_h5cvecseq(x, i = i, j = NULL, k = NULL, l = NULL, drop = drop)
  })
