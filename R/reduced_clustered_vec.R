#' Convert an H5ReducedClusteredVec to a matrix
#'
#' @param x An H5ReducedClusteredVec instance
#' @return A numeric matrix of shape [nTime, nClusters]
#'   with columns named by \code{x@cluster_names} if available.
#' @export
setMethod(
  f = "as.matrix",
  signature = signature(x="H5ReducedClusteredVec"),
  definition = function(x) {
    # Build path: e.g. /scans/<scan_name>/clusters_summary/summary_data
    summary_grp_path <- paste0("/scans/", x@scan_name, "/clusters_summary")
    summary_grp <- x@obj$get(summary_grp_path)

    # The dataset name can be "summary_data" or whatever you choose:
    dset <- summary_grp$get("summary_data")
    mat_data <- dset[]  # read entire dataset => shape [nTime, nClusters]

    # Assign column names if cluster_names are available
    if (!is.null(x@cluster_names) &&
        length(x@cluster_names) == ncol(mat_data)) {
      colnames(mat_data) <- x@cluster_names
    }
    mat_data
  }
)

#' Convert an H5ReducedClusteredVec to a data.frame
#'
#' @param x An H5ReducedClusteredVec instance
#' @param row.names Optional row names (usually NULL)
#' @param optional Passed to \code{as.data.frame} internals
#' @param ... Additional args for \code{as.data.frame}
#' @return A \code{data.frame} with nTime rows, nClusters columns
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



#' Convert an H5ReducedClusteredVecSeq to a list of matrices or a single bound matrix
#'
#' @param x An \code{H5ReducedClusteredVecSeq}
#' @param bind_scans Logical: if TRUE, attempt to row-bind the scans if they have the same number of columns.
#' @return Either a list of matrices or a single matrix if \code{bind_scans=TRUE}.
#' @export
setMethod(
  f = "as.matrix",
  signature = signature(x="H5ReducedClusteredVecSeq"),
  definition = function(x, bind_scans=FALSE) {
    out_list <- list()
    for (sname in x@scan_names) {
      # Construct an H5ReducedClusteredVec on-the-fly for each scan:
      tmp_obj <- new("H5ReducedClusteredVec",
                     obj = x@obj,
                     scan_name = sname,
                     mask = x@mask,
                     clusters = x@clusters,
                     cluster_names = x@cluster_names,
                     cluster_ids = x@cluster_ids,
                     n_time = NA_real_)
      out_list[[sname]] <- as.matrix(tmp_obj)
    }
    if (!bind_scans) {
      return(out_list)
    } else {
      # Attempt to row-bind. We assume each matrix has same ncol => #clusters
      # We add a column "ScanName" or we do direct row bind
      # e.g.:
      # check col counts
      ncols_vec <- sapply(out_list, ncol)
      if (length(unique(ncols_vec)) != 1) {
        stop("Cannot bind scans with differing number of cluster columns.")
      }
      # row-bind:
      full_mat <- do.call(rbind, out_list)
      return(full_mat)
    }
  }
)

#' Convert an H5ReducedClusteredVecSeq to a list of data.frames or a single combined data.frame
#'
#' @param x An H5ReducedClusteredVecSeq
#' @param bind_scans Logical: if TRUE, we row-bind the scans.
#' @param ... further args to as.data.frame
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
      # row-bind. We'll add a 'ScanName' column first for identification
      out_dfs <- list()
      for (nm in names(df_list)) {
        df <- df_list[[nm]]
        df$ScanName <- nm
        out_dfs[[nm]] <- df
      }
      # now do.call(rbind, out_dfs)
      big_df <- do.call(rbind, out_dfs)
      # optionally reorder columns => c("ScanName", rest)
      # ...
      return(big_df)
    }
  }
)


