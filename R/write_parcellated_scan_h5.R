#' Write Single Parcellated Scan to HDF5
#'
#' @description
#' Writes a single parcellated/clustered scan to an HDF5 file. Unlike
#' \code{write_parcellated_experiment_h5} which creates a multi-scan container,
#' this function creates a file optimized for a single scan.
#'
#' @param filepath Path to the output HDF5 file
#' @param mask A \code{LogicalNeuroVol} object defining the brain mask
#' @param clusters A \code{ClusteredNeuroVol} object defining the parcellation
#' @param scan_data Either a named list (for full data) or matrix (for summary data):
#'   - For full data: list with cluster names as keys, each containing a voxels x time matrix
#'   - For summary data: time x clusters matrix
#' @param scan_name Name for the scan (default: "scan_001")
#' @param data_type Type of data: "full" for voxel-level, "summary" for cluster means
#' @param scan_metadata Optional list of metadata for the scan
#' @param cluster_metadata Optional data.frame with cluster information
#' @param compress Logical, whether to compress the data (default: TRUE)
#' @param verbose Logical, whether to print progress messages (default: FALSE)
#'
#' @return An \code{H5ParcellatedScanSummary} or \code{H5ParcellatedScan} object
#'
#' @export
write_parcellated_scan_h5 <- function(filepath, 
                                      mask, 
                                      clusters,
                                      scan_data,
                                      scan_name = "scan_001",
                                      data_type = c("summary", "full"),
                                      scan_metadata = list(),
                                      cluster_metadata = NULL,
                                      compress = TRUE,
                                      verbose = FALSE) {
  
  # Validate inputs
  data_type <- match.arg(data_type)
  
  if (!is(mask, "LogicalNeuroVol")) {
    stop("'mask' must be a LogicalNeuroVol object")
  }
  
  if (!is(clusters, "ClusteredNeuroVol")) {
    stop("'clusters' must be a ClusteredNeuroVol object")
  }
  
  # Check dimensions match
  if (!identical(dim(mask), dim(clusters))) {
    stop("Dimensions of mask and clusters must match")
  }
  
  # Validate scan_data based on type
  if (data_type == "summary") {
    if (!is.matrix(scan_data)) {
      stop("For summary data, scan_data must be a time x clusters matrix")
    }
    n_time <- nrow(scan_data)
    n_clusters <- ncol(scan_data)
  } else {
    if (!is.list(scan_data)) {
      stop("For full data, scan_data must be a list of matrices")
    }
    # Check that all elements are matrices
    if (!all(sapply(scan_data, is.matrix))) {
      stop("All elements of scan_data list must be matrices")
    }
    n_time <- ncol(scan_data[[1]])  # Assuming voxels x time
    n_clusters <- length(scan_data)
  }
  
  # Get cluster information
  unique_clusters <- sort(unique(clusters@clusters[clusters@clusters > 0]))
  
  if (length(unique_clusters) != n_clusters) {
    stop(sprintf("Number of clusters in data (%d) doesn't match clusters in ClusteredNeuroVol (%d)",
                 n_clusters, length(unique_clusters)))
  }
  
  # Generate cluster metadata if not provided
  if (is.null(cluster_metadata)) {
    cluster_metadata <- data.frame(
      cluster_id = unique_clusters,
      n_voxels = sapply(unique_clusters, function(id) sum(clusters@clusters == id))
    )
  }
  
  # Create HDF5 file
  h5f <- NULL
  gzip_level <- if (compress) 4L else 0L
  
  tryCatch({
    if (verbose) message("Creating HDF5 file: ", filepath)
    h5f <- hdf5r::H5File$new(filepath, mode = "w")
    
    # Add class attribute for type detection
    if (data_type == "summary") {
      hdf5r::h5attr(h5f, "fmristore_class") <- "H5ParcellatedScanSummary"
    } else {
      hdf5r::h5attr(h5f, "fmristore_class") <- "H5ParcellatedScan"
    }
    hdf5r::h5attr(h5f, "fmristore_version") <- as.character(packageVersion("fmristore"))
    hdf5r::h5attr(h5f, "scan_name") <- scan_name
    hdf5r::h5attr(h5f, "n_time") <- as.integer(n_time)
    
    # Write global structures
    if (verbose) message("Writing global structures (mask, clusters, header)...")
    
    # Mask
    h5_write(h5f, "/mask", as.array(mask), 
             dtype = hdf5r::h5types$H5T_NATIVE_UCHAR, overwrite = TRUE)
    
    # Cluster map
    h5_write(h5f, "/cluster_map", clusters@clusters, 
             dtype = hdf5r::h5types$H5T_NATIVE_INT32, overwrite = TRUE)
    
    # Voxel coordinates
    h5_write(h5f, "/voxel_coords", which(as.array(mask), arr.ind = TRUE),
             dtype = hdf5r::h5types$H5T_NATIVE_INT32, overwrite = TRUE)
    
    # Header information (NIfTI-like)
    sp <- neuroim2::space(mask)
    dims_vol <- dim(sp)
    hdr_dim <- c(4L, dims_vol[1], dims_vol[2], dims_vol[3], n_time, 1L, 1L, 1L)
    hdr_pixdim <- c(0.0, neuroim2::spacing(sp)[1], neuroim2::spacing(sp)[2], 
                    neuroim2::spacing(sp)[3], 0.0, 0.0, 0.0, 0.0)
    
    q_info <- tryCatch(neuroim2::matrixToQuatern(sp@trans), error = function(e) NULL)
    h5_write(h5f, "/header/dim", hdr_dim, overwrite = TRUE)
    h5_write(h5f, "/header/pixdim", hdr_pixdim, overwrite = TRUE)
    h5_write(h5f, "/header/quatern_b", q_info$qb %||% 0.0, overwrite = TRUE)
    h5_write(h5f, "/header/quatern_c", q_info$qc %||% 0.0, overwrite = TRUE)
    h5_write(h5f, "/header/quatern_d", q_info$qd %||% 0.0, overwrite = TRUE)
    h5_write(h5f, "/header/qoffset_x", q_info$qx %||% 0.0, overwrite = TRUE)
    h5_write(h5f, "/header/qoffset_y", q_info$qy %||% 0.0, overwrite = TRUE)
    h5_write(h5f, "/header/qoffset_z", q_info$qz %||% 0.0, overwrite = TRUE)
    h5_write(h5f, "/header/qfac", q_info$qfac %||% 1.0, overwrite = TRUE)
    
    # Write cluster metadata
    if (verbose) message("Writing cluster metadata...")
    if (!is.null(cluster_metadata) && nrow(cluster_metadata) > 0) {
      h5_write(h5f, "/clusters/cluster_ids", cluster_metadata$cluster_id,
               dtype = hdf5r::h5types$H5T_NATIVE_INT32, overwrite = TRUE)
      
      # Write other columns as separate datasets
      for (col_name in setdiff(names(cluster_metadata), "cluster_id")) {
        col_data <- cluster_metadata[[col_name]]
        if (is.numeric(col_data)) {
          h5_write(h5f, paste0("/clusters/", col_name), col_data, overwrite = TRUE)
        } else if (is.character(col_data) || is.factor(col_data)) {
          h5_write(h5f, paste0("/clusters/", col_name), as.character(col_data), overwrite = TRUE)
        }
      }
    }
    
    # Write scan data under /scans/<scan_name> for consistency with multi-scan format
    if (verbose) message("Writing scan data (type: ", data_type, ")...")
    
    # Create scans group
    scans_grp <- h5f$create_group("scans")
    scan_grp <- scans_grp$create_group(scan_name)
    
    if (data_type == "summary") {
      # Write summary data under scan group
      summary_path <- paste0("/scans/", scan_name, "/clusters_summary/summary_data")
      h5_write(h5f, summary_path, scan_data,
               dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE,
               chunk_dims = c(min(n_time, 100L), n_clusters),
               compression = gzip_level,
               overwrite = TRUE)
      
    } else {
      # Write full voxel data per cluster under scan group
      for (i in seq_along(scan_data)) {
        cluster_id <- unique_clusters[i]
        dset_path <- paste0("/scans/", scan_name, "/clusters/cluster_", cluster_id)
        
        h5_write(h5f, dset_path, scan_data[[i]],
                 dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE,
                 chunk_dims = c(nrow(scan_data[[i]]), min(n_time, 100L)),
                 compression = gzip_level,
                 overwrite = TRUE)
      }
    }
    
    # Write scan metadata if provided
    if (length(scan_metadata) > 0) {
      if (verbose) message("Writing scan metadata...")
      for (key in names(scan_metadata)) {
        h5_write(h5f, paste0("/scans/", scan_name, "/metadata/", key), 
                 scan_metadata[[key]], overwrite = TRUE)
      }
    }
    
    # Close file
    h5f$close_all()
    
    # Return appropriate object based on data type
    if (data_type == "summary") {
      H5ParcellatedScanSummary(
        file = filepath,
        scan_name = scan_name,
        mask = mask,
        clusters = clusters,
        cluster_names = paste0("cluster_", unique_clusters),
        cluster_ids = unique_clusters,
        summary_dset = "summary_data"
      )
    } else {
      H5ParcellatedScan(
        file = filepath,
        scan_name = scan_name,
        mask = mask,
        clusters = clusters
      )
    }
    
  }, error = function(e) {
    if (!is.null(h5f)) {
      try(h5f$close_all(), silent = TRUE)
    }
    stop("Error writing HDF5 file: ", e$message)
  })
}