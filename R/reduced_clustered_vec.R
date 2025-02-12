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
    # Build path => e.g. /scans/<scan_name>/clusters_summary/summary_data
    summary_grp_path <- paste0("/scans/", x@scan_name, "/clusters_summary")
    summary_grp <- x@obj$get(summary_grp_path)

    # Typically "summary_data" is the dataset name
    dset <- summary_grp$get("summary_data")
    mat_data <- dset[]  # entire dataset => shape [nTime, nClusters]

    # If cluster_names are available, set column names
    if (!is.null(x@cluster_names) &&
        length(x@cluster_names) == ncol(mat_data)) {
      colnames(mat_data) <- x@cluster_names
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
    for (sname in x@scan_names) {
      # Construct an H5ReducedClusteredVec on-the-fly for that scan
      tmp_obj <- new("H5ReducedClusteredVec",
                     obj           = x@obj,
                     scan_name     = sname,
                     mask          = x@mask,
                     clusters      = x@clusters,
                     cluster_names = x@cluster_names,
                     cluster_ids   = x@cluster_ids,
                     n_time        = NA_real_)
      out_list[[sname]] <- as.matrix(tmp_obj)
    }

    if (!bind_scans) {
      return(out_list)
    } else {
      # Attempt to row-bind. Must have the same # of columns (#clusters).
      ncols_vec <- sapply(out_list, ncol)
      if (length(unique(ncols_vec)) != 1) {
        stop("Cannot bind scans with differing number of cluster columns.")
      }
      full_mat <- do.call(rbind, out_list)
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

#' Pretty Printer for H5ReducedClusteredVecSeq
#'
#' @description
#' Displays a concise but nicely formatted summary of an \code{H5ReducedClusteredVecSeq},
#' including information about the scans, clusters, and data location.
#'
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

