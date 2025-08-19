#' @include all_generic.R
#' @include all_class.R
NULL

# Helper Functions ----

#' Summarize Data by Clusters
#' 
#' Apply a summary function to voxel data within each cluster
#' 
#' @param nvec A NeuroVec object containing the voxel data
#' @param clusters A ClusteredNeuroVol object defining the clusters
#' @param fun The summary function to apply (default: mean)
#' @return A matrix with clusters as rows and time points as columns
#' @keywords internal
summarize_by_clusters <- function(nvec, clusters, fun = mean) {
  # Get mask from clusters since it defines the voxel space
  nvec_mask <- clusters@mask
  
  # Extract voxel data - need to handle the masking properly
  # series() returns a matrix with time as rows, voxels as columns
  voxel_indices <- which(nvec_mask)
  voxel_data <- series(nvec, voxel_indices)
  # Transpose to get voxels as rows, time as columns
  voxel_data <- t(voxel_data)
  
  # Get unique clusters (excluding 0 which represents non-clustered voxels)
  unique_clusters <- sort(unique(clusters@clusters[clusters@clusters > 0]))
  
  # Apply summary function to each cluster
  # Result should be time x clusters
  cluster_summaries <- sapply(unique_clusters, function(cid) {
    cluster_voxels <- which(clusters@clusters == cid)
    if (length(cluster_voxels) > 0) {
      # Apply function across voxels for each time point
      apply(voxel_data[cluster_voxels, , drop = FALSE], 2, fun, na.rm = TRUE)
    } else {
      rep(NA, ncol(voxel_data))
    }
  })
  
  # Return as time x clusters matrix
  return(cluster_summaries)
}

#' Extract Full Cluster Data
#' 
#' Extract voxel-level data organized by clusters
#' 
#' @param nvec A NeuroVec object containing the voxel data
#' @param clusters A ClusteredNeuroVol object defining the clusters
#' @return A list with one matrix per cluster (voxels x time)
#' @keywords internal
extract_cluster_data <- function(nvec, clusters) {
  # Get mask from clusters since it defines the voxel space
  nvec_mask <- clusters@mask
  
  # Extract voxel data
  voxel_indices <- which(nvec_mask)
  voxel_data <- series(nvec, voxel_indices)
  # Transpose to get voxels as rows, time as columns
  voxel_data <- t(voxel_data)
  
  # Get unique clusters (excluding 0)
  unique_clusters <- sort(unique(clusters@clusters[clusters@clusters > 0]))
  
  # Extract data for each cluster
  cluster_list <- lapply(unique_clusters, function(cid) {
    # cluster_voxels are indices within the clusters@clusters vector
    cluster_mask_indices <- which(clusters@clusters == cid)
    if (length(cluster_mask_indices) > 0) {
      # We need to map from cluster indices to voxel_data row indices
      # Since voxel_data has all mask voxels in order, and clusters@clusters
      # also follows mask order, the indices should match directly
      voxel_data[cluster_mask_indices, , drop = FALSE]
    } else {
      matrix(numeric(0), nrow = 0, ncol = ncol(voxel_data))
    }
  })
  
  names(cluster_list) <- paste0("cluster_", unique_clusters)
  return(cluster_list)
}

# Write Dataset Methods ----

#' Write Dataset Method for NeuroVec
#' 
#' Write a single NeuroVec object with clustering to HDF5
#' 
#' @param x A NeuroVec object
#' @param file Output HDF5 file path
#' @param clusters A ClusteredNeuroVol object defining spatial clusters
#' @param mask Optional LogicalNeuroVol mask (extracted from x if not provided)
#' @param summary Logical, whether to summarize data by clusters (default: FALSE)
#' @param summary_fun Function to use for summarization (default: mean)
#' @param scan_name Name for the scan (default: "scan_001")
#' @param compression Compression level 0-9 (default: 4)
#' @param ... Additional arguments passed to write_parcellated_experiment_h5
#' @return Invisible NULL
#' 
#' @examples
#' \dontrun{
#' # Write full voxel-level data
#' write_dataset(nvec, file = "output.h5", clusters = cvol)
#' 
#' # Write summarized data (mean per cluster)
#' write_dataset(nvec, file = "output_summary.h5", clusters = cvol, summary = TRUE)
#' 
#' # Use median for summarization
#' write_dataset(nvec, file = "output_median.h5", clusters = cvol, 
#'               summary = TRUE, summary_fun = median)
#' }
#' 
#' @export
#' @rdname write_dataset
setMethod("write_dataset", 
  signature(x = "NeuroVec"),
  function(x, file, clusters, mask = NULL, summary = FALSE, 
           summary_fun = mean, scan_name = "scan_001",
           compression = 4, ...) {
    
    # Validate inputs
    if (missing(file)) stop("'file' argument is required")
    if (missing(clusters)) stop("'clusters' argument is required")
    if (!is(clusters, "ClusteredNeuroVol")) {
      stop("'clusters' must be a ClusteredNeuroVol object")
    }
    
    # Extract mask if not provided
    if (is.null(mask)) {
      # Use the mask from the clusters object
      mask <- clusters@mask
      
      if (!is(mask, "LogicalNeuroVol")) {
        stop("Could not extract a valid LogicalNeuroVol mask from the clusters object")
      }
    }
    
    # Convert to runs_data format
    if (summary) {
      # Apply summary function to each cluster
      data <- summarize_by_clusters(x, clusters, summary_fun)
      type <- "summary"
    } else {
      # Keep full voxel-level data
      data <- extract_cluster_data(x, clusters)
      type <- "full"
    }
    
    # Create runs_data structure
    runs_data <- list(
      list(
        scan_name = scan_name,
        type = type,
        data = data
      )
    )
    
    # Delegate to write_parcellated_experiment_h5
    # Note: write_parcellated_experiment_h5 uses 'compress' not 'compression'
    # compression levels: 0 = none, 1-9 = increasing compression
    # convert to boolean for backward compatibility
    write_parcellated_experiment_h5(
      filepath = file,
      mask = mask,
      clusters = clusters,
      runs_data = runs_data,
      compress = compression > 0,
      ...
    )
    
    invisible(NULL)
  }
)

#' Write Dataset Method for List of NeuroVecs
#' 
#' Write multiple NeuroVec objects with clustering to HDF5
#' 
#' @param x A list of NeuroVec objects
#' @param file Output HDF5 file path
#' @param clusters A ClusteredNeuroVol object defining spatial clusters
#' @param mask Optional LogicalNeuroVol mask (extracted from first element if not provided)
#' @param summary Logical, whether to summarize data by clusters (default: FALSE)
#' @param summary_fun Function to use for summarization (default: mean)
#' @param scan_names Character vector of scan names (auto-generated if not provided)
#' @param compression Compression level 0-9 (default: 4)
#' @param ... Additional arguments passed to write_parcellated_experiment_h5
#' @return Invisible NULL
#' 
#' @examples
#' \dontrun{
#' # Write multiple scans
#' scans <- list(nvec1, nvec2, nvec3)
#' write_dataset(scans, file = "multi_scan.h5", clusters = cvol)
#' 
#' # With custom names and summary
#' write_dataset(scans, file = "multi_summary.h5", clusters = cvol,
#'               summary = TRUE, scan_names = c("rest", "task1", "task2"))
#' }
#' 
#' @export
#' @rdname write_dataset
setMethod("write_dataset",
  signature(x = "list"),
  function(x, file, clusters, mask = NULL, summary = FALSE,
           summary_fun = mean, scan_names = NULL, compression = 4, ...) {
    
    # Validate inputs
    if (missing(file)) stop("'file' argument is required")
    if (missing(clusters)) stop("'clusters' argument is required")
    if (!is(clusters, "ClusteredNeuroVol")) {
      stop("'clusters' must be a ClusteredNeuroVol object")
    }
    
    # Validate all elements are NeuroVec
    if (length(x) == 0) stop("List cannot be empty")
    if (!all(sapply(x, function(v) is(v, "NeuroVec")))) {
      stop("All elements must be NeuroVec objects")
    }
    
    # Check dimensions match (first 3 dimensions should be spatial)
    dims <- lapply(x, dim)
    spatial_dims <- lapply(dims, function(d) d[1:3])
    if (!all(sapply(spatial_dims[-1], identical, spatial_dims[[1]]))) {
      stop("All NeuroVec objects must have the same spatial dimensions")
    }
    
    # Extract mask from first element if not provided
    if (is.null(mask)) {
      mask <- mask(x[[1]])
      if (!is(mask, "LogicalNeuroVol")) {
        stop("Could not extract a valid LogicalNeuroVol mask from the first NeuroVec object")
      }
    }
    
    # Generate scan names if not provided
    if (is.null(scan_names)) {
      scan_names <- sprintf("scan_%03d", seq_along(x))
    } else if (length(scan_names) != length(x)) {
      stop("Length of scan_names must match length of input list")
    }
    
    # Convert each NeuroVec to runs_data format
    runs_data <- lapply(seq_along(x), function(i) {
      if (summary) {
        data <- summarize_by_clusters(x[[i]], clusters, summary_fun)
        type <- "summary"
      } else {
        data <- extract_cluster_data(x[[i]], clusters)
        type <- "full"
      }
      
      list(
        scan_name = scan_names[i],
        type = type,
        data = data
      )
    })
    
    # Delegate to write_parcellated_experiment_h5
    # Note: write_parcellated_experiment_h5 uses 'compress' not 'compression'
    # compression levels: 0 = none, 1-9 = increasing compression
    # convert to boolean for backward compatibility
    write_parcellated_experiment_h5(
      filepath = file,
      mask = mask,
      clusters = clusters,
      runs_data = runs_data,
      compress = compression > 0,
      ...
    )
    
    invisible(NULL)
  }
)

#' Write Dataset Method for H5ParcellatedMultiScan
#' 
#' Export an existing H5ParcellatedMultiScan object to a new HDF5 file
#' 
#' @param x An H5ParcellatedMultiScan object
#' @param file Output HDF5 file path
#' @param overwrite Logical, whether to overwrite existing file (default: FALSE)
#' @param ... Additional arguments (currently unused)
#' @return Invisible file path
#' 
#' @examples
#' \dontrun{
#' # Export existing H5 object to new file
#' h5_obj <- H5ParcellatedMultiScan("existing.h5")
#' write_dataset(h5_obj, file = "copy.h5")
#' }
#' 
#' @export
#' @rdname write_dataset
setMethod("write_dataset",
  signature(x = "H5ParcellatedMultiScan"),
  function(x, file, overwrite = FALSE, ...) {
    
    if (missing(file)) stop("'file' argument is required")
    
    # Get the source file
    src_file <- h5file(x)
    
    # Check if source and destination are the same
    if (normalizePath(src_file, mustWork = FALSE) == 
        normalizePath(file, mustWork = FALSE)) {
      message("Object already backed by target file: ", file)
      return(invisible(file))
    }
    
    # Check overwrite
    if (!overwrite && file.exists(file)) {
      stop("Output file already exists: ", file, 
           "\nUse overwrite=TRUE to replace")
    }
    
    message("Copying H5ParcellatedMultiScan from ", src_file, " to ", file)
    
    # Use file.copy for now (could implement more sophisticated H5-to-H5 copy)
    success <- file.copy(src_file, file, overwrite = overwrite)
    
    if (!success) {
      stop("Failed to copy file from ", src_file, " to ", file)
    }
    
    invisible(file)
  }
)

#' Write Dataset Method for H5NeuroVec
#' 
#' Export an existing H5NeuroVec object to a new HDF5 file
#' 
#' @param x An H5NeuroVec object  
#' @param file Output HDF5 file path
#' @param overwrite Logical, whether to overwrite existing file (default: FALSE)
#' @param ... Additional arguments (currently unused)
#' @return Invisible file path
#' 
#' @export
#' @rdname write_dataset
setMethod("write_dataset",
  signature(x = "H5NeuroVec"),
  function(x, file, overwrite = FALSE, ...) {
    
    if (missing(file)) stop("'file' argument is required")
    
    # Get the source file
    src_file <- h5file(x)
    
    # Check if source and destination are the same
    if (normalizePath(src_file, mustWork = FALSE) == 
        normalizePath(file, mustWork = FALSE)) {
      message("Object already backed by target file: ", file)
      return(invisible(file))
    }
    
    # Check overwrite
    if (!overwrite && file.exists(file)) {
      stop("Output file already exists: ", file, 
           "\nUse overwrite=TRUE to replace")
    }
    
    message("Copying H5NeuroVec from ", src_file, " to ", file)
    
    # Use file.copy for now
    success <- file.copy(src_file, file, overwrite = overwrite)
    
    if (!success) {
      stop("Failed to copy file from ", src_file, " to ", file)
    }
    
    invisible(file)
  }
)

#' Write Dataset Method for H5NeuroVol
#' 
#' Export an existing H5NeuroVol object to a new HDF5 file
#' 
#' @param x An H5NeuroVol object
#' @param file Output HDF5 file path
#' @param overwrite Logical, whether to overwrite existing file (default: FALSE)
#' @param ... Additional arguments (currently unused)
#' @return Invisible file path
#' 
#' @export
#' @rdname write_dataset
setMethod("write_dataset",
  signature(x = "H5NeuroVol"),
  function(x, file, overwrite = FALSE, ...) {
    
    if (missing(file)) stop("'file' argument is required")
    
    # Get the source file
    src_file <- h5file(x)
    
    # Check if source and destination are the same
    if (normalizePath(src_file, mustWork = FALSE) == 
        normalizePath(file, mustWork = FALSE)) {
      message("Object already backed by target file: ", file)
      return(invisible(file))
    }
    
    # Check overwrite
    if (!overwrite && file.exists(file)) {
      stop("Output file already exists: ", file, 
           "\nUse overwrite=TRUE to replace")
    }
    
    message("Copying H5NeuroVol from ", src_file, " to ", file)
    
    # Use file.copy for now
    success <- file.copy(src_file, file, overwrite = overwrite)
    
    if (!success) {
      stop("Failed to copy file from ", src_file, " to ", file)
    }
    
    invisible(file)
  }
)

#' Write Dataset Method for LatentNeuroVec
#' 
#' Write a LatentNeuroVec object to HDF5
#' 
#' @param x A LatentNeuroVec object
#' @param file Output HDF5 file path
#' @param compression Compression level 0-9 (default: 6)
#' @param ... Additional arguments passed to as_h5
#' @return Invisible NULL
#' 
#' @examples
#' \dontrun{
#' # Write latent representation
#' write_dataset(latent_nvec, file = "latent.h5")
#' }
#' 
#' @export
#' @rdname write_dataset
setMethod("write_dataset",
  signature(x = "LatentNeuroVec"),
  function(x, file, compression = 6, ...) {
    
    if (missing(file)) stop("'file' argument is required")
    
    # Delegate to existing as_h5 method
    result <- as_h5(x, file = file, compression = compression, ...)
    
    # Close the file handle if returned
    if (is(result, "H5File")) {
      result$close_all()
    }
    
    invisible(NULL)
  }
)