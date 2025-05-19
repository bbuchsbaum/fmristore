#' @include all_class.R io_h5_helpers.R
#' @import hdf5r
#' @import neuroim2
#' @importFrom methods is
#' @importFrom utils write.table
#' @importFrom hdf5r H5T_STRING H5S h5types
NULL

# Contains functions for writing fmristore HDF5 structures.

# Helper function to validate the runs_data list structure
validate_runs_data <- function(rd) {
  for (i in seq_along(rd)) {
    el <- rd[[i]]
    stopifnot(
      is.list(el),
      is.character(el$scan_name), length(el$scan_name) == 1, nzchar(el$scan_name),
      el$type %in% c("full", "summary"),
      is.list(el$data) || is.matrix(el$data) # Allow list for full, matrix for summary
      # Add check: if type == "full", data must be list; if type == "summary", data must be matrix
      # Add check: if type == "full", names(data) should be cluster_XXX
    )
    if (el$type == "full" && !is.list(el$data)) stop(sprintf("Run %d ('%s'): type is 'full' but data is not a list.", i, el$scan_name))
    if (el$type == "summary" && !is.matrix(el$data)) stop(sprintf("Run %d ('%s'): type is 'summary' but data is not a matrix.", i, el$scan_name))
    # TODO: Could add validation that full data list names match cluster IDs expected
  }
}

#' Write Clustered Experiment Data to HDF5
#'
#' @description
#' Writes neuroimaging data structured according to the H5ClusterExperiment
#' specification into an HDF5 file.
#'
#' This function takes R objects representing the mask, cluster definitions,
#' run-specific data (either full voxel-level or summary time series), and
#' associated metadata, and creates the corresponding HDF5 groups and datasets.
#'
#' @param filepath Character string: the path to the HDF5 file to create.
#'   If the file exists, it will be overwritten.
#' @param mask A `LogicalNeuroVol` object representing the brain mask.
#' @param clusters A `ClusteredNeuroVol` object containing cluster assignments
#'   for voxels within the mask.
#' @param runs_data A list where each element represents a single run/scan.
#'   Each element must be a list containing:
#'   \itemize{
#'     \item `scan_name`: (character) Unique identifier for the scan.
#'     \item `type`: (character) Either "full" or "summary".
#'     \item `data`: 
#'       \itemize{
#'         \item If `type` is "full", `data` must be a list where names are `cluster_<cid>`
#'               (e.g., `cluster_1`, `cluster_2`) and values are matrices 
#'               `[nVoxelsInCluster, nTime]` containing the time series for that cluster.
#'         \item If `type` is "summary", `data` must be a single matrix
#'               `[nTime, nClusters]` containing the summary time series.
#'       }
#'     \item `metadata`: (Optional) A list of key-value pairs for scan-specific metadata.
#'                     Can include `n_time` explicitly, otherwise it's inferred from data.
#'   }
#' @param cluster_metadata (Optional) A `data.frame` containing metadata for the clusters.
#'   Must contain at least a column named `cluster_id` matching the unique IDs in `clusters`.
#'   Other columns will be written as part of a compound dataset.
#' @param overwrite Logical: If `TRUE`, overwrite the file if it exists. Default `FALSE`.
#' @param compress Logical: If `TRUE`, apply GZIP compression to data arrays. Default `TRUE`.
#' @param verbose Logical: Print progress messages? Default `TRUE`.
#'
#' @return Invisibly returns `NULL`. Called for its side effect of creating the HDF5 file.
#' @export
#' @family H5ClusteredIO
#' @examples
#' if (requireNamespace("neuroim2", quietly = TRUE) &&
#'     requireNamespace("hdf5r", quietly = TRUE) &&
#'     exists("write_clustered_experiment_h5", where = "package:fmristore") &&
#'     !is.null(fmristore:::create_minimal_LogicalNeuroVol) &&
#'     !is.null(fmristore:::create_minimal_ClusteredNeuroVol)) {
#'
#'   temp_h5_file <- NULL
#'   
#'   tryCatch({
#'     # 1. Create a temporary file path
#'     temp_h5_file <- tempfile(fileext = ".h5")
#'     
#'     # 2. Create minimal mask and clusters using helpers
#'     mask_vol <- fmristore:::create_minimal_LogicalNeuroVol(dims = c(5L, 5L, 2L))
#'     # Ensure clusters are within the mask and have some content
#'     # Create clusters that align with the mask's space
#'     clust_vol <- fmristore:::create_minimal_ClusteredNeuroVol(
#'       space = neuroim2::space(mask_vol), # Use mask's space
#'       mask = mask_vol@.Data,            # Use mask's data
#'       num_clusters = 2L
#'     )
#'     
#'     # 3. Prepare minimal runs_data
#'     # Get cluster IDs and number of voxels per cluster from clust_vol
#'     unique_cids <- sort(unique(clust_vol@clusters[clust_vol@clusters > 0]))
#'     n_time_run1 <- 10L
#'     n_time_run2 <- 8L
#'     
#'     # Run 1: Full data type
#'     run1_data_list <- list()
#'     if (length(unique_cids) > 0) {
#'       for (cid in unique_cids) {
#'         n_vox_in_cluster <- sum(clust_vol@clusters == cid)
#'         if (n_vox_in_cluster > 0) {
#'            # Data: nVoxInCluster x nTime
#'           run1_data_list[[paste0("cluster_", cid)]] <- matrix(
#'             rnorm(n_vox_in_cluster * n_time_run1), 
#'             nrow = n_vox_in_cluster, 
#'             ncol = n_time_run1
#'           )
#'         }
#'       }
#'     }
#'     
#'     run1 <- list(
#'       scan_name = "ScanA_Full",
#'       type = "full",
#'       data = run1_data_list,
#'       metadata = list(subject_id = "sub-01", task = "rest", n_time = n_time_run1)
#'     )
#'     
#'     # Run 2: Summary data type
#'     # Data: nTime x nClusters
#'     run2_summary_matrix <- matrix(
#'       rnorm(n_time_run2 * length(unique_cids)), 
#'       nrow = n_time_run2, 
#'       ncol = length(unique_cids)
#'     )
#'     colnames(run2_summary_matrix) <- paste0("cluster_", unique_cids) # Optional: for clarity
#'     
#'     run2 <- list(
#'       scan_name = "ScanB_Summary",
#'       type = "summary",
#'       data = run2_summary_matrix,
#'       metadata = list(subject_id = "sub-01", task = "task", n_time = n_time_run2)
#'     )
#'     
#'     runs_data_list <- list(run1, run2)
#'     
#'     # 4. Prepare minimal cluster_metadata (optional)
#'     cluster_meta_df <- NULL
#'     if (length(unique_cids) > 0) {
#'       cluster_meta_df <- data.frame(
#'         cluster_id = unique_cids,
#'         name = paste0("Region_", LETTERS[1:length(unique_cids)]),
#'         size_vox = sapply(unique_cids, function(id) sum(clust_vol@clusters == id))
#'       )
#'     }
#'     
#'     # 5. Call the function
#'     write_clustered_experiment_h5(
#'       filepath = temp_h5_file,
#'       mask = mask_vol,
#'       clusters = clust_vol,
#'       runs_data = runs_data_list,
#'       cluster_metadata = cluster_meta_df,
#'       overwrite = TRUE,
#'       verbose = FALSE 
#'     )
#'     
#'     # Verify file was created
#'     if (file.exists(temp_h5_file)) {
#'       cat("Successfully wrote clustered experiment to:", temp_h5_file, "\\n")
#'       # Optional: Basic check of the HDF5 file structure
#'       # h5f <- hdf5r::H5File$new(temp_h5_file, mode="r")
#'       # print(h5f$ls(recursive=TRUE))
#'       # h5f$close_all()
#'     }
#'     
#'   }, error = function(e) {
#'     message("write_clustered_experiment_h5 example failed: ", e$message)
#'     if (!is.null(temp_h5_file)) message("Temporary file was: ", temp_h5_file)
#'   }, finally = {
#'     # Clean up temporary file
#'     if (!is.null(temp_h5_file) && file.exists(temp_h5_file)) {
#'       unlink(temp_h5_file)
#'     }
#'   })
#' } else {
#'   message("Skipping write_clustered_experiment_h5 example: dependencies or helpers not available.")
#' }
write_clustered_experiment_h5 <- function(filepath,
                                          mask,
                                          clusters,
                                          runs_data,
                                          cluster_metadata = NULL,
                                          overwrite = FALSE,
                                          compress = TRUE, # Added compress argument
                                          verbose = TRUE) {

  # --- Input Validation --- 
  if (!is(mask, "LogicalNeuroVol")) stop("`mask` must be a LogicalNeuroVol object.")
  if (!is(clusters, "ClusteredNeuroVol")) stop("`clusters` must be a ClusteredNeuroVol object.")
  # Use helper for dimension check
  check_same_dims(mask, clusters, dims_to_compare = 1:3, 
                  msg = "Dimensions of mask and clusters must match.")
  if (!is.list(runs_data)) stop("`runs_data` must be a list.")
  if (file.exists(filepath) && !overwrite) stop("File exists and overwrite is FALSE: ", filepath)
  if (file.exists(filepath) && overwrite) file.remove(filepath)

  n_vox_mask <- sum(mask)
  if (length(clusters@clusters) != n_vox_mask) {
      stop(sprintf("Length of clusters vector (%d) does not match number of voxels in mask (%d).", length(clusters@clusters), n_vox_mask))
  }
  
  # Validate runs_data structure
  validate_runs_data(runs_data)
  
  # Validate cluster_metadata if provided
  if (!is.null(cluster_metadata)) {
      if (!is.data.frame(cluster_metadata)) stop("`cluster_metadata` must be a data.frame.")
      if (!("cluster_id" %in% names(cluster_metadata))) stop("`cluster_metadata` must contain a 'cluster_id' column.")
      # TODO: Check if cluster_ids in metadata match unique(clusters@clusters)?
  }

  # --- File Creation --- 
  h5f <- NULL
  gzip_level <- if (compress) 4L else 0L # Set gzip level based on argument
  
  tryCatch({
    if (verbose) message("Creating HDF5 file: ", filepath)
    h5f <- H5File$new(filepath, mode = "w")
    
    # --- Write Global Structures --- 
    if (verbose) message("Writing global structures (mask, clusters, header)... ")
    # Mask - use h5_write
    h5_write(h5f, "/mask", as.array(mask), dtype = h5types$H5T_NATIVE_UCHAR, overwrite = TRUE)
    
    # Cluster Map - use h5_write
    h5_write(h5f, "/cluster_map", clusters@clusters, dtype = h5types$H5T_NATIVE_INT32, overwrite = TRUE)
    
    # Voxel Coordinates - use h5_write
    h5_write(h5f, "/voxel_coords", which(as.array(mask), arr.ind = TRUE), 
             dtype = h5types$H5T_NATIVE_INT32, overwrite = TRUE)

    # Header - use h5_write for each field
    # hdr_grp <- h5f$create_group("header") # h5_write creates parent
    sp <- space(mask)
    dims_vol <- dim(sp)
    hdr_dim <- c(4L, dims_vol[1], dims_vol[2], dims_vol[3], length(runs_data), 1L, 1L, 1L) 
    hdr_pixdim <- c(0.0, spacing(sp)[1], spacing(sp)[2], spacing(sp)[3], 0.0, 0.0, 0.0, 0.0)
    q_info <- tryCatch(neuroim2::matrixToQuatern(sp@trans), error=function(e) NULL)
    h5_write(h5f, "/header/dim", hdr_dim, overwrite = TRUE)
    h5_write(h5f, "/header/pixdim", hdr_pixdim, overwrite = TRUE)
    h5_write(h5f, "/header/quatern_b", q_info$qb %||% 0.0, overwrite = TRUE)
    h5_write(h5f, "/header/quatern_c", q_info$qc %||% 0.0, overwrite = TRUE)
    h5_write(h5f, "/header/quatern_d", q_info$qd %||% 0.0, overwrite = TRUE)
    h5_write(h5f, "/header/qoffset_x", q_info$qx %||% 0.0, overwrite = TRUE)
    h5_write(h5f, "/header/qoffset_y", q_info$qy %||% 0.0, overwrite = TRUE)
    h5_write(h5f, "/header/qoffset_z", q_info$qz %||% 0.0, overwrite = TRUE)
    h5_write(h5f, "/header/qfac", q_info$qfac %||% 1.0, overwrite = TRUE)
    # hdr_grp$close() # No longer needed

    # Global Clusters Group
    global_clus_grp <- h5f$create_group("clusters")
    unique_cluster_ids <- sort(unique(clusters@clusters))
    # Write cluster_ids using h5_write
    h5_write(global_clus_grp, "cluster_ids", unique_cluster_ids, 
             dtype = h5types$H5T_NATIVE_INT32, overwrite = TRUE)
    
    # Write cluster_metadata using h5_write
    if (!is.null(cluster_metadata)) {
      if (verbose) message("Writing global cluster metadata...")
      # meta_grp <- global_clus_grp$create_group("cluster_meta") # h5_write creates parent
      tryCatch({
          cluster_metadata_filtered <- cluster_metadata[cluster_metadata$cluster_id %in% unique_cluster_ids, , drop = FALSE]
          if (nrow(cluster_metadata_filtered) != nrow(cluster_metadata)) {
              warning("Provided cluster_metadata contained IDs not present in the clusters object; only metadata for existing IDs was written.")
          }
          for (cn in names(cluster_metadata_filtered)) {
            vec <- cluster_metadata_filtered[[cn]]
            col_path <- file.path("/clusters/cluster_meta", cn)
            h5_write(h5f, col_path, vec, 
                     dtype = guess_h5_type(vec), # guess_h5_type handles string type creation
                     overwrite = TRUE)
          }
      }, error = function(e) {
          warning("Failed during writing of cluster metadata: ", e$message)
      }) # No finally needed for meta_grp
    }
    # global_clus_grp$close() # Close group handle

    # --- Write Scans --- 
    if (verbose) message("Writing scan data...")
    scans_grp <- h5f$create_group("scans")

    # Validate that all runs_data have the same type before proceeding
    if (length(runs_data) > 0) {
      all_run_types <- vapply(runs_data, function(run) run$type, character(1))
      unique_run_types <- unique(all_run_types)
      if (length(unique_run_types) > 1) {
        stop(sprintf("All items in 'runs_data' must have the same 'type' (either 'full' or 'summary'). Found mixed types: %s", 
                     paste(unique_run_types, collapse=", ")))
      }
      # All types are the same, so the first one is representative
      summary_only <- (unique_run_types[1] == "summary")
    } else {
      # No runs, default summary_only to FALSE (or could be NA/absent, but FALSE is safer for attribute)
      summary_only <- FALSE 
    }
    hdf5r::h5attr(scans_grp, "summary_only") <- as.logical(summary_only)
    
    for (i in seq_along(runs_data)) {
        run <- runs_data[[i]]
        sname <- run$scan_name
        stype <- run$type
        sdata <- run$data
        smeta <- run$metadata %||% list()
        
        if (verbose) message(sprintf("  Writing scan: %s (type: %s)", sname, stype))
        # scan_grp <- scans_grp$create_group(sname) # h5_write creates parent
        
        # Write metadata using h5_write
        if (length(smeta) > 0) {
           # meta_grp <- scan_grp$create_group("metadata") # h5_write creates parent
           tryCatch({
             for (mname in names(smeta)) {
                mval <- smeta[[mname]]
                meta_path <- file.path("/scans", sname, "metadata", mname)
                h5_write(h5f, meta_path, mval, 
                         dtype = guess_h5_type(mval), 
                         overwrite = TRUE)
             }
           }, error = function(e) {
               warning(sprintf("Failed during writing of metadata for scan '%s': %s", sname, e$message))
           }) # No finally needed
        }
        
        # Write data based on type
        if (stype == "full") {
           # scan_clus_grp <- scan_grp$create_group("clusters") # h5_write creates parent
           tryCatch({
               # --- Determine nTime ---
               nTime <- NA_integer_
               # 1. Prioritize from metadata if available and valid
               if (!is.null(smeta$n_time) && is.numeric(smeta$n_time) && length(smeta$n_time) == 1 && smeta$n_time > 0) {
                   nTime <- as.integer(smeta$n_time)
                   if (verbose) message(sprintf("    Scan '%s': Using nTime %d from metadata.", sname, nTime))
               }

               # 2. If not in metadata, infer from data matrices in sdata
               if (is.na(nTime)) {
                   if (verbose) message(sprintf("    Scan '%s': nTime not in metadata, attempting to infer from sdata list.", sname))
                   if (length(sdata) > 0) {
                       for (data_item_name in names(sdata)) {
                           data_item <- sdata[[data_item_name]]
                           if (is.matrix(data_item) && ncol(data_item) > 0) {
                               nTime <- ncol(data_item)
                               if (verbose) message(sprintf("    Scan '%s': Inferred nTime %d from sdata item '%s'.", sname, nTime, data_item_name))
                               break # Found nTime, exit loop
                           }
                       }
                   }
               }
               
               if (is.na(nTime)) {
                   stop(sprintf("Could not determine a valid nTime for full data write in scan '%s'. Provide in metadata or ensure sdata contains valid matrices.", sname))
               }
               
               for (nm in names(sdata)) { # e.g. "cluster_17"
                 mat <- sdata[[nm]]
                 if (!is.matrix(mat) || !is.numeric(mat))
                   stop("run ", sname, ": ", nm, " is not numeric matrix")
                 dims <- dim(mat)
                 if (any(dims < 0) || length(dims) != 2) stop("Invalid dimensions for matrix in ", nm)
                 
                 chunk_dims_full <- if(prod(dims)>0) pmin(dims, c(1024L, 128L)) else NULL
                 data_path <- file.path("/scans", sname, "clusters", nm)
                 h5_write(h5f, data_path, mat, 
                          dtype = h5types$H5T_IEEE_F32LE, 
                          chunk_dims = chunk_dims_full,
                          compression = if(prod(dims)>0) gzip_level else 0L,
                          overwrite = TRUE)
               }
           }, error = function(e) {
               stop(sprintf("Failed writing full data for scan '%s': %s", sname, e$message))
           }) # No finally needed
           
        } else if (stype == "summary") {
           # scan_summary_grp <- scan_grp$create_group("clusters_summary") # h5_write creates parent
           tryCatch({
               mat <- sdata # validated above
               if (!is.matrix(mat) || !is.numeric(mat))
                 stop("run ", sname, ": summary data is not numeric matrix")
               dims <- dim(mat)
               if (any(dims < 0) || length(dims) != 2) stop("Invalid dimensions for summary matrix")
               
               # Check matrix cols match unique cluster IDs
               if (dims[2] != length(unique_cluster_ids)) {
                   stop(sprintf("Summary matrix for scan '%s' has %d columns, but %d unique cluster IDs exist.", sname, dims[2], length(unique_cluster_ids)))
               }
               
               chunk_dims_summary <- if(prod(dims)>0) pmin(dims, c(128L, 256L)) else NULL
               summary_path <- file.path("/scans", sname, "clusters_summary", "summary_data")
               h5_write(h5f, summary_path, mat, 
                        dtype = h5types$H5T_IEEE_F32LE,
                        chunk_dims = chunk_dims_summary,
                        compression = if(prod(dims)>0) gzip_level else 0L,
                        overwrite = TRUE)
           }, error = function(e) {
               stop(sprintf("Failed writing summary data for scan '%s': %s", sname, e$message))
           }) # No finally needed
           
        } else {
           warning(sprintf("Unknown run type '%s' for scan '%s'. Skipping data write.", stype, sname))
        }
        
        # No need to close scan_grp explicitly
    }
    # No need to close scans_grp explicitly

  }, error = function(e) {
     if (!is.null(h5f) && h5f$is_valid) try(h5f$close_all(), silent=TRUE)
     if (file.exists(filepath)) file.remove(filepath) # Clean up partial file on error
     stop("Failed during HDF5 file writing: ", e$message)
  }, finally = {
     if (!is.null(h5f) && h5f$is_valid) try(h5f$close_all(), silent=TRUE)
  })
  
  invisible(NULL)
} 