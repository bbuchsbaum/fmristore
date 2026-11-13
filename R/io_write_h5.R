#' @include all_class.R io_h5_helpers.R
#' @importFrom hdf5r H5File is.h5file H5S H5D H5T_STRING h5attr h5attr_names
#' @importFrom neuroim2 NeuroSpace space spacing origin trans matrixToQuatern
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

    # Validate that full data list names match cluster_XXX pattern
    if (el$type == "full" && is.list(el$data)) {
      data_names <- names(el$data)
      if (!is.null(data_names) && length(data_names) > 0) {
        cluster_pattern <- "^cluster_[0-9]+$"
        if (!all(grepl(cluster_pattern, data_names))) {
          invalid_names <- data_names[!grepl(cluster_pattern, data_names)]
          stop(sprintf(
            "Run %d ('%s'): full data list contains invalid cluster names: %s. Expected pattern: cluster_<ID>",
            i, el$scan_name, paste(invalid_names, collapse = ", ")
          ))
        }
      }
    }
  }
}

# Helper: check each run's data against the clusters and its metadata before
# anything is written. Returns the number of time points of each run.
validate_runs_data_content <- function(rd, clusters) {
  scan_names <- vapply(rd, function(el) el$scan_name, character(1))
  if (anyDuplicated(scan_names) > 0) {
    stop(sprintf(
      "Scan names in `runs_data` must be unique. Duplicated: %s",
      paste(unique(scan_names[duplicated(scan_names)]), collapse = ", ")
    ))
  }

  cluster_ids <- sort(unique(clusters@clusters))
  cluster_sizes <- tabulate(match(clusters@clusters, cluster_ids), nbins = length(cluster_ids))
  # Full data must cover every labelled cluster; ID 0 marks unclustered voxels.
  required_ids <- cluster_ids[cluster_ids > 0]

  n_time <- integer(length(rd))
  for (i in seq_along(rd)) {
    el <- rd[[i]]
    run_label <- sprintf("Run %d ('%s')", i, el$scan_name)
    sdata <- el$data

    if (el$type == "summary") {
      if (!is.numeric(sdata)) stop(run_label, ": summary data is not a numeric matrix.")
      if (ncol(sdata) != length(cluster_ids)) {
        stop(sprintf(
          "%s: summary matrix has %d columns, but %d unique cluster IDs exist.",
          run_label, ncol(sdata), length(cluster_ids)
        ))
      }
      n_time[i] <- nrow(sdata)
      next
    }

    data_names <- names(sdata)
    if (length(sdata) > 0 && (is.null(data_names) || anyDuplicated(data_names) > 0)) {
      stop(run_label, ": full data must be a list with unique names of the form cluster_<ID>.")
    }
    data_ids <- as.integer(sub("^cluster_", "", data_names))
    unknown <- data_names[!(data_ids %in% cluster_ids)]
    if (length(unknown) > 0) {
      stop(sprintf(
        "%s: full data contains clusters not present in `clusters`: %s",
        run_label, paste(unknown, collapse = ", ")
      ))
    }
    missing_ids <- setdiff(required_ids, data_ids)
    if (length(missing_ids) > 0) {
      stop(sprintf(
        "%s: full data is missing cluster(s): %s",
        run_label, paste0("cluster_", missing_ids, collapse = ", ")
      ))
    }
    for (k in seq_along(sdata)) {
      mat <- sdata[[k]]
      if (!is.matrix(mat) || !is.numeric(mat)) {
        stop(sprintf("%s: %s is not a numeric matrix.", run_label, data_names[k]))
      }
      n_vox_cluster <- cluster_sizes[match(data_ids[k], cluster_ids)]
      if (nrow(mat) != n_vox_cluster) {
        stop(sprintf(
          "%s: %s has %d rows, but cluster %d has %d voxels in the mask.",
          run_label, data_names[k], nrow(mat), data_ids[k], n_vox_cluster
        ))
      }
    }
    n_cols <- unique(vapply(sdata, ncol, integer(1)))
    if (length(n_cols) > 1) {
      stop(sprintf(
        "%s: all full data matrices must have the same number of columns (time points); found %s.",
        run_label, paste(sort(n_cols), collapse = ", ")
      ))
    }

    meta_n_time <- el$metadata$n_time
    has_meta_n_time <- !is.null(meta_n_time) && is.numeric(meta_n_time) &&
      length(meta_n_time) == 1 && meta_n_time > 0
    if (length(n_cols) == 1 && has_meta_n_time && meta_n_time != n_cols) {
      stop(sprintf(
        "%s: metadata n_time (%s) does not match the %d columns of the full data matrices.",
        run_label, format(meta_n_time), n_cols
      ))
    }
    if (length(n_cols) == 1 && n_cols > 0) {
      n_time[i] <- n_cols
    } else if (has_meta_n_time) {
      n_time[i] <- as.integer(meta_n_time)
    } else {
      stop(sprintf(
        paste(
          "Could not determine a valid nTime for full data write in scan '%s'.",
          "Provide in metadata or ensure sdata contains valid matrices."
        ),
        el$scan_name
      ))
    }
  }
  n_time
}

#' Write Parcellated Experiment Data to HDF5
#'
#' @description
#' Writes neuroimaging data structured according to the H5ParcellatedMultiScan
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
#' \dontrun{
#' if (requireNamespace("neuroim2", quietly = TRUE) &&
#'   requireNamespace("hdf5r", quietly = TRUE) &&
#'   exists("write_parcellated_experiment_h5", where = "package:fmristore") &&
#'   !is.null(fmristore:::create_minimal_LogicalNeuroVol) &&
#'   !is.null(fmristore:::create_minimal_ClusteredNeuroVol)) {
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
#'       mask = mask_vol@.Data, # Use mask's data
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
#'           # Data: nVoxInCluster x nTime
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
#'     write_parcellated_experiment_h5(
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
#'       cat("Successfully wrote parcellated experiment to:", temp_h5_file, "\\n")
#'       # Optional: Basic check of the HDF5 file structure
#'       # h5f <- hdf5r::H5File$new(temp_h5_file, mode="r")
#'       # print(h5f$ls(recursive=TRUE))
#'       # h5f$close_all()
#'     }
#'   }, error = function(e) {
#'     message("write_parcellated_experiment_h5 example failed: ", e$message)
#'     if (!is.null(temp_h5_file)) message("Temporary file was: ", temp_h5_file)
#'   }, finally = {
#'     # Clean up temporary file
#'     if (!is.null(temp_h5_file) && file.exists(temp_h5_file)) {
#'       unlink(temp_h5_file)
#'     }
#'   })
#' } else {
#'   message("Skipping example: dependencies not available.")
#' }
#' }
#' @export
write_parcellated_experiment_h5 <- function(filepath,
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
  check_same_dims(mask, clusters,
    dims_to_compare = 1:3,
    msg = "Dimensions of mask and clusters must match."
  )
  if (!is.list(runs_data)) stop("`runs_data` must be a list.")
  if (file.exists(filepath) && !overwrite) stop("File exists and overwrite is FALSE: ", filepath)

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

    # Check if cluster_ids in metadata match unique(clusters@clusters)
    unique_cluster_ids_actual <- sort(unique(clusters@clusters))
    metadata_cluster_ids <- sort(unique(cluster_metadata$cluster_id))

    if (!identical(metadata_cluster_ids, unique_cluster_ids_actual)) {
      missing_in_metadata <- setdiff(unique_cluster_ids_actual, metadata_cluster_ids)
      extra_in_metadata <- setdiff(metadata_cluster_ids, unique_cluster_ids_actual)

      error_msg <- "Cluster IDs in metadata do not match those in clusters object."
      if (length(missing_in_metadata) > 0) {
        error_msg <- paste0(error_msg, " Missing from metadata: ", paste(missing_in_metadata, collapse = ", "), ".")
      }
      if (length(extra_in_metadata) > 0) {
        error_msg <- paste0(error_msg, " Extra in metadata: ", paste(extra_in_metadata, collapse = ", "), ".")
      }
      stop(error_msg)
    }
  }

  # Validate run data against the clusters and metadata; yields nTime per run
  run_n_time <- validate_runs_data_content(runs_data, clusters)

  # --- File Creation ---
  h5f <- NULL
  gzip_level <- if (compress) 4L else 0L # Set gzip level based on argument

  # Only replace an existing file once all input has been validated
  if (file.exists(filepath) && overwrite) file.remove(filepath)

  tryCatch({
    if (verbose) message("Creating HDF5 file: ", filepath)
    h5f <- hdf5r::H5File$new(filepath, mode = "w")

    # Add class attribute for type detection
    hdf5r::h5attr(h5f, "fmristore_class") <- "H5ParcellatedMultiScan"
    hdf5r::h5attr(h5f, "fmristore_version") <- as.character(packageVersion("fmristore"))

    # --- Write Global Structures ---
    if (verbose) message("Writing global structures (mask, clusters, header)... ")
    # Mask - use h5_write
    h5_write(h5f, "/mask", as.array(mask), dtype = h5types$H5T_NATIVE_UCHAR, overwrite = TRUE)

    # Cluster Map - use h5_write
    h5_write(h5f, "/cluster_map", clusters@clusters, dtype = h5types$H5T_NATIVE_INT32, overwrite = TRUE)

    # Voxel Coordinates - use h5_write
    h5_write(h5f, "/voxel_coords", which(as.array(mask), arr.ind = TRUE),
      dtype = h5types$H5T_NATIVE_INT32, overwrite = TRUE
    )

    # Header - use h5_write for each field
    # hdr_grp <- h5f$create_group("header") # h5_write creates parent
    sp <- space(mask)
    dims_vol <- dim(sp)
    hdr_dim <- c(4L, dims_vol[1], dims_vol[2], dims_vol[3], length(runs_data), 1L, 1L, 1L)
    hdr_pixdim <- c(0.0, spacing(sp)[1], spacing(sp)[2], spacing(sp)[3], 0.0, 0.0, 0.0, 0.0)
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
    # hdr_grp$close() # No longer needed

    # Global Clusters Group
    global_clus_grp <- h5f$create_group("clusters")
    unique_cluster_ids <- sort(unique(clusters@clusters))
    # Write cluster_ids using h5_write
    h5_write(h5f, "/clusters/cluster_ids", unique_cluster_ids,
      dtype = h5types$H5T_NATIVE_INT32, overwrite = TRUE
    )

    # Write cluster_metadata using h5_write
    if (!is.null(cluster_metadata)) {
      if (verbose) message("Writing global cluster metadata...")
      # meta_grp <- global_clus_grp$create_group("cluster_meta") # h5_write creates parent
      tryCatch(
        {
          cluster_metadata_filtered <- cluster_metadata[cluster_metadata$cluster_id %in% unique_cluster_ids, , drop = FALSE]
          if (nrow(cluster_metadata_filtered) != nrow(cluster_metadata)) {
            warning("Provided cluster_metadata contained IDs not present in the clusters object; only metadata for existing IDs was written.")
          }
          for (cn in names(cluster_metadata_filtered)) {
            vec <- cluster_metadata_filtered[[cn]]
            col_path <- file.path("/clusters/cluster_meta", cn)
            h5_write(h5f, col_path, vec,
              dtype = guess_h5_type(vec), # guess_h5_type handles string type creation
              overwrite = TRUE
            )
          }
        },
        error = function(e) {
          warning("Failed during writing of cluster metadata: ", e$message)
        }
      ) # No finally needed for meta_grp
    }
    # global_clus_grp$close() # Close group handle

    # --- Write Scans ---
    if (verbose) message("Writing scan data...")
    scans_grp <- h5f$create_group("scans")

    # Check if all runs are summary-only
    if (length(runs_data) > 0) {
      all_run_types <- vapply(runs_data, function(run) run$type, character(1))
      # Set summary_only attribute based on whether ALL runs are summary type
      summary_only <- all(all_run_types == "summary")
    } else {
      # No runs, default summary_only to FALSE
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
      scan_grp <- scans_grp$create_group(sname)

      # Set class attribute based on type
      if (stype == "full") {
        hdf5r::h5attr(scan_grp, "class") <- "H5ParcellatedScan"
      } else if (stype == "summary") {
        hdf5r::h5attr(scan_grp, "class") <- "H5ParcellatedScanSummary"
      }

      # Write metadata using h5_write
      if (length(smeta) > 0) {
        # meta_grp <- scan_grp$create_group("metadata") # h5_write creates parent
        tryCatch(
          {
            for (mname in names(smeta)) {
              mval <- smeta[[mname]]
              meta_path <- file.path("/scans", sname, "metadata", mname)
              h5_write(h5f, meta_path, mval,
                dtype = guess_h5_type(mval),
                overwrite = TRUE
              )
            }
          },
          error = function(e) {
            warning(sprintf("Failed during writing of metadata for scan '%s': %s", sname, e$message))
          }
        ) # No finally needed
      }

      # Write data based on type
      if (stype == "full") {
        # scan_clus_grp <- scan_grp$create_group("clusters") # h5_write creates parent
        tryCatch(
          {
            # nTime and the matrices were checked by validate_runs_data_content()
            if (verbose) message(sprintf("    Scan '%s': nTime %d.", sname, run_n_time[i]))

            for (nm in names(sdata)) { # e.g. "cluster_17"
              mat <- sdata[[nm]]
              if (!is.matrix(mat) || !is.numeric(mat)) {
                stop("run ", sname, ": ", nm, " is not numeric matrix")
              }
              dims <- dim(mat)
              if (any(dims < 0) || length(dims) != 2) stop("Invalid dimensions for matrix in ", nm)

              chunk_dims_full <- if (prod(dims) > 0) pmin(dims, c(1024L, 128L)) else NULL
              data_path <- file.path("/scans", sname, "clusters", nm)
              h5_write(h5f, data_path, mat,
                dtype = h5types$H5T_IEEE_F32LE,
                chunk_dims = chunk_dims_full,
                compression = if (prod(dims) > 0) gzip_level else 0L,
                overwrite = TRUE
              )
            }
          },
          error = function(e) {
            stop(sprintf("Failed writing full data for scan '%s': %s", sname, e$message))
          }
        ) # No finally needed
      } else if (stype == "summary") {
        # scan_summary_grp <- scan_grp$create_group("clusters_summary") # h5_write creates parent
        tryCatch(
          {
            mat <- sdata # validated above
            if (!is.matrix(mat) || !is.numeric(mat)) {
              stop("run ", sname, ": summary data is not numeric matrix")
            }
            dims <- dim(mat)
            if (any(dims < 0) || length(dims) != 2) stop("Invalid dimensions for summary matrix")

            # Check matrix cols match unique cluster IDs
            if (dims[2] != length(unique_cluster_ids)) {
              stop(sprintf("Summary matrix for scan '%s' has %d columns, but %d unique cluster IDs exist.", sname, dims[2], length(unique_cluster_ids)))
            }

            chunk_dims_summary <- if (prod(dims) > 0) pmin(dims, c(128L, 256L)) else NULL
            summary_path <- file.path("/scans", sname, "clusters_summary", "summary_data")
            h5_write(h5f, summary_path, mat,
              dtype = h5types$H5T_IEEE_F32LE,
              chunk_dims = chunk_dims_summary,
              compression = if (prod(dims) > 0) gzip_level else 0L,
              overwrite = TRUE
            )
          },
          error = function(e) {
            stop(sprintf("Failed writing summary data for scan '%s': %s", sname, e$message))
          }
        ) # No finally needed
      } else {
        warning(sprintf("Unknown run type '%s' for scan '%s'. Skipping data write.", stype, sname))
      }

      # No need to close scan_grp explicitly
    }
    # No need to close scans_grp explicitly
  }, error = function(e) {
    if (!is.null(h5f) && h5f$is_valid) try(h5f$close_all(), silent = TRUE)
    if (file.exists(filepath)) file.remove(filepath) # Clean up partial file on error
    stop("Failed during HDF5 file writing: ", e$message)
  }, finally = {
    if (!is.null(h5f) && h5f$is_valid) try(h5f$close_all(), silent = TRUE)
  })

  invisible(NULL)
}

#' @rdname write_parcellated_experiment_h5
#' @export
write_clustered_experiment_h5 <- function(filepath,
                                          mask,
                                          clusters,
                                          runs_data,
                                          cluster_metadata = NULL,
                                          overwrite = FALSE,
                                          compress = TRUE,
                                          verbose = TRUE) {
  .Deprecated("write_parcellated_experiment_h5",
    package = "fmristore",
    msg = "write_clustered_experiment_h5 is deprecated. Please use write_parcellated_experiment_h5 instead."
  )

  write_parcellated_experiment_h5(
    filepath = filepath,
    mask = mask,
    clusters = clusters,
    runs_data = runs_data,
    cluster_metadata = cluster_metadata,
    overwrite = overwrite,
    compress = compress,
    verbose = verbose
  )
}
