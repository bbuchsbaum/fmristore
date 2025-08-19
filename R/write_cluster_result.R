#' Write Clustering Results to HDF5
#'
#' A convenience function for writing clustering results (e.g., from supervoxels())
#' along with the original neuroimaging data to an HDF5 file. This function serves
#' as a user-friendly interface to \code{\link{write_dataset}}.
#'
#' @param cluster_result Either a \code{ClusteredNeuroVol} object or a list containing:
#'   \itemize{
#'     \item \code{clusters}: Numeric vector of cluster assignments for each voxel within the mask
#'     \item \code{mask}: A \code{LogicalNeuroVol} defining which voxels were included in clustering
#'     \item \code{metadata}: (Optional) Additional cluster information
#'   }
#' @param neurovec A \code{NeuroVec} object containing the 4D fMRI time series data
#' @param filename Character string specifying the output HDF5 file path
#' @param scan_name Character string naming this scan in the HDF5 file (default: "scan_001")
#' @param type Character string, either "full" for voxel-level data or "summary" for
#'   cluster-averaged data
#' @param cluster_metadata Optional data.frame containing metadata for each cluster.
#'   Should have a column named \code{cluster_id} matching the unique cluster IDs.
#' @param overwrite Logical, whether to overwrite an existing file (default: FALSE)
#' @param compress Logical or numeric. If logical: TRUE uses compression level 4,
#'   FALSE uses no compression. If numeric: specify compression level 0-9.
#' @param verbose Logical, whether to print progress messages (default: TRUE)
#' @param ... Additional arguments passed to \code{\link{write_dataset}}
#'
#' @return Invisibly returns NULL. Called for its side effect of creating an HDF5 file.
#'
#' @details
#' This function bridges clustering results from various sources (e.g., \code{supervoxels()},
#' custom clustering algorithms) with the fmristore HDF5 storage system. It automatically
#' converts the clustering information into the required \code{ClusteredNeuroVol} format
#' and delegates to \code{write_dataset()} for efficient HDF5 storage.
#'
#' When \code{type = "full"}, the complete voxel-level time series are stored organized
#' by cluster. When \code{type = "summary"}, only the cluster-averaged time series are
#' stored, resulting in much smaller file sizes.
#'
#' @examples
#' \dontrun{
#' # Example with clustering results from supervoxels()
#' clust_result <- supervoxels(brain_data, n_clusters = 100)
#'
#' # Write full voxel-level data
#' write_cluster_result(
#'   cluster_result = clust_result,
#'   neurovec = my_4d_data,
#'   filename = "clustered_data.h5",
#'   type = "full"
#' )
#'
#' # Write only cluster averages (more compact)
#' write_cluster_result(
#'   cluster_result = clust_result,
#'   neurovec = my_4d_data,
#'   filename = "clustered_summary.h5",
#'   type = "summary",
#'   compress = 6 # Higher compression
#' )
#'
#' # With cluster metadata
#' cluster_info <- data.frame(
#'   cluster_id = 1:100,
#'   region = rep(c("frontal", "parietal", "temporal", "occipital"), 25),
#'   size = sapply(1:100, function(i) sum(clust_result$clusters == i))
#' )
#'
#' write_cluster_result(
#'   cluster_result = clust_result,
#'   neurovec = my_4d_data,
#'   filename = "clustered_with_metadata.h5",
#'   cluster_metadata = cluster_info,
#'   type = "full"
#' )
#' }
#'
#' @seealso
#' \code{\link{write_dataset}} for the underlying implementation
#' \code{\link{read_dataset}} for loading the saved data
#' \code{\link[neuroim2]{ClusteredNeuroVol}} for the cluster representation
#'
#' @export
#' @importFrom neuroim2 ClusteredNeuroVol
#' @importFrom methods is
write_cluster_result <- function(cluster_result,
                                 neurovec,
                                 filename,
                                 scan_name = "scan_001",
                                 type = c("full", "summary"),
                                 cluster_metadata = NULL,
                                 overwrite = FALSE,
                                 compress = TRUE,
                                 verbose = TRUE,
                                 ...) {
  # Validate type argument
  type <- match.arg(type)

  # Validate neurovec
  if (!is(neurovec, "NeuroVec")) {
    stop("neurovec must be a NeuroVec object")
  }

  # Handle cluster_result structure
  if (is(cluster_result, "ClusteredNeuroVol")) {
    # Already a ClusteredNeuroVol - use directly
    clusters <- cluster_result
  } else if (is.list(cluster_result)) {
    # Convert from list format (e.g., from supervoxels() or custom clustering)
    if (!all(c("clusters", "mask") %in% names(cluster_result))) {
      stop("cluster_result list must contain 'clusters' and 'mask' elements")
    }

    # Validate mask
    if (!is(cluster_result$mask, "LogicalNeuroVol")) {
      stop("cluster_result$mask must be a LogicalNeuroVol object")
    }

    # Validate clusters vector
    if (!is.numeric(cluster_result$clusters)) {
      stop("cluster_result$clusters must be a numeric vector")
    }

    # Check that clusters vector length matches mask
    n_mask_voxels <- sum(cluster_result$mask@.Data)
    if (length(cluster_result$clusters) != n_mask_voxels) {
      stop(sprintf(
        "Length of clusters vector (%d) doesn't match number of TRUE voxels in mask (%d)",
        length(cluster_result$clusters), n_mask_voxels
      ))
    }

    # Create ClusteredNeuroVol from components
    clusters <- ClusteredNeuroVol(
      mask = cluster_result$mask,
      clusters = cluster_result$clusters
    )

    # Transfer any metadata if present
    if (!is.null(cluster_result$metadata) && is.null(cluster_metadata)) {
      cluster_metadata <- cluster_result$metadata
    }
  } else {
    stop("cluster_result must be a ClusteredNeuroVol object or a list with 'clusters' and 'mask' elements")
  }

  # Handle compression parameter
  if (is.logical(compress)) {
    compression <- ifelse(compress, 4, 0)
  } else if (is.numeric(compress)) {
    compression <- as.integer(compress)
    if (compression < 0 || compression > 9) {
      stop("Compression level must be between 0 and 9")
    }
  } else {
    stop("compress must be logical (TRUE/FALSE) or numeric (0-9)")
  }

  # Check file overwrite
  if (file.exists(filename) && !overwrite) {
    stop("File already exists: ", filename, "\nSet overwrite=TRUE to replace")
  }

  # Progress messages
  if (verbose) {
    message("Writing cluster results to: ", filename)
    message("  Scan name: ", scan_name)
    message("  Type: ", type, " (", ifelse(type == "full", "voxel-level", "cluster-averaged"), " data)")
    message("  Compression: ", if (compression > 0) paste0("level ", compression) else "none")

    # Report cluster statistics
    unique_clusters <- unique(clusters@clusters[clusters@clusters > 0])
    message("  Number of clusters: ", length(unique_clusters))
    message("  Total voxels: ", sum(clusters@mask@.Data))
  }

  # Prepare additional arguments
  dots <- list(...)

  # Remove any duplicate arguments that we're explicitly setting
  dots$cluster_metadata <- NULL
  dots$overwrite <- NULL
  dots$verbose <- NULL

  # Call write_dataset with appropriate parameters
  write_dataset(
    x = neurovec,
    file = filename,
    clusters = clusters,
    scan_name = scan_name,
    summary = (type == "summary"),
    compression = compression,
    cluster_metadata = cluster_metadata,
    overwrite = overwrite,
    verbose = verbose,
    ...
  )

  if (verbose) {
    message("Successfully wrote cluster results to: ", filename)
  }

  invisible(NULL)
}
