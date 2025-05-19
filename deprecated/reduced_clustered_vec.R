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

#' Convert an H5ReducedClusteredVec to a matrix
#'
#' @description
#' Reads the summary dataset \code{/scans/<scan_name>/clusters_summary/summary_data}
#' from HDF5, returning it as a numeric matrix \code{[nTime, nClusters]}. If
#' \code{cluster_names} are available, they become the column names.
#'
#' @param x An \code{H5ReducedClusteredVec} instance.
#' @return A numeric matrix of dimension \code{[nTime, nClusters]}.
#' @export
setMethod(
  f = "as.matrix",
  signature = signature(x="H5ReducedClusteredVec"),
  definition = function(x) {
    summary_grp <- NULL
    ds <- NULL
    mat_data <- NULL
    
    tryCatch({
        # Build path => e.g. /scans/<scan_name>/clusters_summary/summary_data
        summary_grp_path <- paste0("/scans/", x@scan_name, "/clusters_summary")

        # Check if group exists
        if (!x@obj$exists(summary_grp_path)) {
            stop(paste0("[as.matrix,H5ReducedClusteredVec] Summary group not found at path: ", summary_grp_path))
        }
        summary_grp <- x@obj[[summary_grp_path]]
        
        # Check if dataset exists within group
        dset_name <- "summary_data"
        if (!summary_grp$exists(dset_name)) {
            stop(paste0("[as.matrix,H5ReducedClusteredVec] Summary dataset '", dset_name, "' not found in group: ", summary_grp_path))
        }
        ds <- summary_grp[[dset_name]] 
        on.exit(if (exists("ds") && !is.null(ds) && inherits(ds, "H5D") && ds$is_valid) ds$close(), add = TRUE)

        # Read entire dataset => shape [nTime, nClusters]
        mat_data <- ds[]  

    }, error = function(e) {
        # Ensure ds handle is attempted to be closed even on error before this point
        if (exists("ds") && !is.null(ds) && inherits(ds, "H5D") && ds$is_valid) {
            try(ds$close(), silent = TRUE) # Attempt closure, suppress errors during cleanup
        }
        stop(paste0("[as.matrix,H5ReducedClusteredVec] Failed to read summary data for scan '", x@scan_name, "'. Original error: ", e$message))
    })

    # If cluster_names are available, set column names
    if (!is.null(mat_data) && !is.null(x@cluster_names) &&
        length(x@cluster_names) == ncol(mat_data)) {
      colnames(mat_data) <- x@cluster_names
    } else if (!is.null(mat_data) && !is.null(x@cluster_names) && length(x@cluster_names) != ncol(mat_data)){
        warning(paste0("[as.matrix,H5ReducedClusteredVec] Length of cluster_names (", 
                      length(x@cluster_names), ") does not match number of columns (", 
                      ncol(mat_data), ") in summary data for scan '", x@scan_name, "'. Column names not set."))
    }
    
    # Check if mat_data was successfully read
    if (is.null(mat_data)) {
        stop("[as.matrix,H5ReducedClusteredVec] Failed to retrieve summary data matrix for scan '", x@scan_name, "'. Result is NULL.")
    }

    mat_data
  }
)

#' Convert an H5ReducedClusteredVec to a data.frame
#'
#' @description
#' Calls \code{as.matrix(x)} and then converts to \code{data.frame}.
#'
#' @param x An \code{H5ReducedClusteredVec} instance.
#' @param row.names Optional row names (usually \code{NULL}).
#' @param optional Passed to \code{\link{as.data.frame}} internals.
#' @param ... Additional arguments for \code{as.data.frame}.
#'
#' @return A \code{data.frame} with \code{nTime} rows and \code{nClusters} columns.
#' @export
setMethod(
  f = "as.data.frame",
  signature = signature(x="H5ReducedClusteredVec"),
  definition = function(x, row.names=NULL, optional=FALSE, ...) {
    mat_data <- as.matrix(x)
    df <- as.data.frame(mat_data, row.names=row.names, optional=optional, ...)
    df
  }
)

#' Convert an H5ReducedClusteredVecSeq to matrices
#'
#' @description
#' For each scan in the sequence, constructs an \code{H5ReducedClusteredVec} on the fly,
#' reads its summary dataset, and returns the result(s).
#'
#' @param x An \code{H5ReducedClusteredVecSeq} object (multiple scans).
#' @param bind_scans \code{logical}. If \code{FALSE}, returns a \code{list} of matrices,
#'   one per scan. If \code{TRUE} and all matrices have the same number of columns
#'   (clusters), row-binds them into one large matrix.
#'
#' @return A \code{list} of matrices (if \code{bind_scans=FALSE}) or a single matrix
#'   (if \code{bind_scans=TRUE}).
#'
#' @export
setMethod(
  f = "as.matrix",
  signature = signature(x="H5ReducedClusteredVecSeq"),
  definition = function(x, bind_scans=FALSE) {
    out_list <- list()
    has_failures <- FALSE
    for (sname in x@scan_names) {
      tmp_obj <- NULL # Ensure tmp_obj is defined in this scope
      tryCatch({
          # Construct an H5ReducedClusteredVec on-the-fly for that scan
          tmp_obj <- new("H5ReducedClusteredVec",
                         obj           = x@obj,
                         scan_name     = sname,
                         mask          = x@mask,
                         clusters      = x@clusters,
                         cluster_names = x@cluster_names,
                         cluster_ids   = x@cluster_ids,
                         n_time        = NA_real_) # n_time not used by as.matrix
          
          # Call the (now robust) as.matrix for the single scan vec
          out_list[[sname]] <- as.matrix(tmp_obj)
          
      }, error = function(e) {
          warning(paste0("[as.matrix,H5ReducedClusteredVecSeq] Failed to process scan '", 
                       sname, "'. Skipping. Original error: ", e$message))
          has_failures <<- TRUE # Use <<- to modify variable in parent environment
          # Do not add to out_list if it failed
      })
    }
    
    # If no scans were processed successfully, return an empty list or handle appropriately
    if(length(out_list) == 0){
        warning("[as.matrix,H5ReducedClusteredVecSeq] No scans processed successfully.")
        if (!bind_scans) {
            return(list()) # Return empty list
        } else {
            # Return an empty matrix? Or maybe stop? Let's return empty matrix for now.
            return(matrix(numeric(), nrow = 0, ncol = 0)) 
        }
    }

    if (!bind_scans) {
      if (has_failures) {
           warning("[as.matrix,H5ReducedClusteredVecSeq] Returning list of matrices, but some scans failed to process.")
      }
      return(out_list)
    } else {
      # Attempt to row-bind. Must have the same # of columns (#clusters).
      # Filter out any potential NULLs if error handling changed (though currently we just skip)
      valid_matrices <- Filter(Negate(is.null), out_list)
      if (length(valid_matrices) == 0) {
           warning("[as.matrix,H5ReducedClusteredVecSeq] No valid matrices to bind after processing scans.")
           return(matrix(numeric(), nrow = 0, ncol = 0))
      }
      
      ncols_vec <- sapply(valid_matrices, ncol)
      if (length(unique(ncols_vec)) != 1) {
        stop("[as.matrix,H5ReducedClusteredVecSeq] Cannot bind scans: differing number of cluster columns found among successfully processed scans.")
      }
      
      # rbind only the valid matrices
      full_mat <- do.call(rbind, valid_matrices)
      
      if (has_failures) {
          warning("[as.matrix,H5ReducedClusteredVecSeq] Returning bound matrix, but some scans failed to process and were excluded.")
      }
      return(full_mat)
    }
  }
)

#' Convert an H5ReducedClusteredVecSeq to data.frames
#'
#' @description
#' Similar to \code{as.matrix} for \code{H5ReducedClusteredVecSeq}, but each matrix
#' is turned into a \code{data.frame}. If \code{bind_scans=TRUE}, row-binds them into
#' one combined \code{data.frame} (adding a \code{ScanName} column for identification).
#'
#' @param x An \code{H5ReducedClusteredVecSeq} object.
#' @param bind_scans \code{logical} indicating whether to combine all scans into
#'   one big data.frame. If \code{FALSE}, returns a list of data.frames, one per scan.
#' @param ... Further arguments passed to \code{as.data.frame()}.
#'
#' @return A \code{list} of data.frames if \code{bind_scans=FALSE}, or a single
#'   combined data.frame if \code{bind_scans=TRUE}.
#'
#' @export
setMethod(
  f = "as.data.frame",
  signature = signature(x="H5ReducedClusteredVecSeq"),
  definition = function(x, bind_scans=FALSE, ...) {
    mat_list <- as.matrix(x, bind_scans=FALSE)
    df_list  <- lapply(mat_list, as.data.frame, ...)

    if (!bind_scans) {
      return(df_list)
    } else {
      # row-bind. We'll add a 'ScanName' column for each
      out_dfs <- list()
      for (nm in names(df_list)) {
        df <- df_list[[nm]]
        df$ScanName <- nm
        out_dfs[[nm]] <- df
      }
      big_df <- do.call(rbind, out_dfs)
      # Optionally reorder columns => e.g., c("ScanName", rest)
      return(big_df)
    }
  }
)

#' @rdname show-methods
#' @importFrom crayon bold blue silver yellow green magenta
#' @importFrom methods show
#' @export
setMethod(
  f = "show",
  signature = "H5ReducedClusteredVecSeq",
  definition = function(object) {
    # Header line
    cat("\n", crayon::bold(crayon::blue("H5ReducedClusteredVecSeq")), "\n", sep="")

    # Basic info block: how many scans, how many clusters, shape, etc.
    n_scans <- length(object@scan_names)
    n_clusters <- length(object@cluster_ids)

    cat(crayon::silver("────────────────────────────────────────\n"))
    cat(crayon::bold(crayon::yellow("Basic Info")), "\n")
    cat(crayon::silver(" • "), crayon::green("Number of scans:"), n_scans, "\n")
    cat(crayon::silver(" • "), crayon::green("Number of clusters:"), n_clusters, "\n")

    # Time points (n_time slot is a numeric vector specifying time lengths per scan)
    # If all scans share the same time length, we show one number; otherwise we show a range
    if (length(unique(object@n_time)) == 1) {
      cat(crayon::silver(" • "), crayon::green("Time points per scan:"),
          unique(object@n_time), "\n")
    } else {
      cat(crayon::silver(" • "), crayon::green("Time points (per scan):"),
          paste(object@n_time, collapse=", "), "\n")
    }

    # Possibly show a short listing of scan names
    cat(crayon::bold("\nScan Names"), ":\n")
    if (n_scans <= 5) {
      cat("  ", paste(object@scan_names, collapse=", "), "\n")
    } else {
      preview <- paste(head(object@scan_names, 5), collapse=", ")
      remainder <- n_scans - 5
      cat("  ", preview, crayon::silver(paste0("... (", remainder, " more)")), "\n")
    }

    # Show a snippet about the cluster metadata
    cat(crayon::bold("\nCluster Metadata"), ":\n")
    if (nrow(object@cluster_metadata) == 0) {
      cat(crayon::silver("  No cluster_metadata table.\n"))
    } else {
      md_cols <- names(object@cluster_metadata)
      cat(crayon::silver("  Columns:"), paste(md_cols, collapse=", "), "\n")
      # Optionally show first few rows
      cat(crayon::silver("  #rows="), nrow(object@cluster_metadata), "\n")
    }

    # Show if we have cluster_names
    if (length(object@cluster_names) > 0) {
      cat(crayon::bold("\nCluster Names"), ":\n")
      if (length(object@cluster_names) <= 6) {
        cat("  ", paste(object@cluster_names, collapse=", "), "\n")
      } else {
        preview <- paste(head(object@cluster_names, 6), collapse=", ")
        cat("  ", preview, " ...\n")
      }
    }

    # Show info about the underlying file
    cat(crayon::bold("\nStorage:"), "\n")
    if (object@obj$is_valid) {
      cat(crayon::silver(" • "), "HDF5 file: ", crayon::magenta(object@obj$get_filename()), "\n", sep="")
    } else {
      cat(crayon::silver(" • "), "HDF5 file is ", crayon::red("CLOSED"), "\n", sep="")
    }

    cat("\n")
  }
)

#' Create an H5ReducedClusteredVec object (DEPRECATED)
#'
#' @description
#' This constructor is deprecated. Use `make_run_summary()` instead to create
#' an `H5ClusteredRunSummary` object.
#'
#' Originally, this class represented summary data for a single scan.
#' This function now acts as a compatibility wrapper.
#'
#' @param obj An H5File object.
#' @param scan_name The name of the scan.
#' @param mask A LogicalNeuroVol mask.
#' @param clusters A ClusteredNeuroVol object (can be NULL for summary).
#' @param cluster_names Optional character vector of cluster names.
#' @param cluster_ids Optional integer vector of cluster IDs.
#' @param n_time Optional integer number of time points (will be inferred by `make_run_summary` if NULL).
#'
#' @return A new `H5ClusteredRunSummary` instance (via `make_run_summary`).
#' @importFrom lifecycle deprecate_warn
#' @export
H5ReducedClusteredVec <- function(obj, scan_name, mask, clusters = NULL,
                                  cluster_names = character(), cluster_ids = integer(),
                                  n_time = NULL) { # n_time is not directly used by make_run_summary but kept for signature compatibility
                                  
  lifecycle::deprecate_warn(
    when = "0.1.0", # Replace with the version number
    what = "H5ReducedClusteredVec()",
    with = "make_run_summary()",
    details = "The H5ReducedClusteredVec class is being replaced by H5ClusteredRunSummary."
  )

  # Call the new constructor for the summary run
  # Note: n_time is not directly passed as make_run_summary infers it from the dataset
  make_run_summary(
    file_source = obj,
    scan_name = scan_name,
    mask = mask,
    clusters = clusters, # Pass along, might be NULL
    cluster_names = cluster_names,
    cluster_ids = cluster_ids
    # summary_dset uses default "summary_data"
  )
}

#' Create an H5ReducedClusteredVecSeq object (DEPRECATED)
#'
#' @description
#' This constructor is deprecated. Use `H5ClusteredExperiment()` instead, which
#' can represent experiments containing summary-only runs.
#'
#' Originally, this class represented a sequence of summary-only scans.
#' This function now acts as a compatibility wrapper.
#'
#' @param obj An H5File object.
#' @param scan_names A character vector of scan names.
#' @param mask A LogicalNeuroVol mask.
#' @param clusters A ClusteredNeuroVol object (can be NULL for summary).
#' @param cluster_names Optional character vector of cluster names.
#' @param cluster_ids Optional integer vector of cluster IDs.
#' @param n_time (DEPRECATED, ignored) Vector of time points per scan.
#' @param scan_metadata List of metadata per scan.
#' @param cluster_metadata Data frame of cluster metadata.
#' @param cumulative_time (DEPRECATED, ignored) Cumulative time points.
#'
#' @return A new `H5ClusteredExperiment` instance.
#' @importFrom lifecycle deprecate_warn
#' @export
H5ReducedClusteredVecSeq <- function(obj, scan_names, mask, clusters = NULL,
                                       cluster_names = character(), cluster_ids = integer(),
                                       n_time = NULL, # Deprecated arg
                                       scan_metadata = list(),
                                       cluster_metadata = data.frame(),
                                       cumulative_time = NULL) { # Deprecated arg

  lifecycle::deprecate_warn(
    when = "0.1.0", # Replace with the version number
    what = "H5ReducedClusteredVecSeq()",
    with = "H5ClusteredExperiment(summary_preference='require')", # Suggest preference
    details = "The H5ReducedClusteredVecSeq class is replaced by H5ClusteredExperiment. Use summary_preference to control run types."
  )

  # Call the H5ClusteredExperiment constructor
  # Set summary_preference to ensure summary runs are expected/loaded
  # Pass other relevant arguments
  H5ClusteredExperiment(
    file_source = obj,
    scan_names = scan_names,
    mask = mask,
    clusters = clusters, # Pass along, might be NULL
    scan_metadata = scan_metadata,
    cluster_metadata = cluster_metadata,
    summary_preference = "require" # Or "prefer", depending on desired strictness
    # H5ClusteredExperiment constructor ignores n_time and cumulative_time
  )
}

