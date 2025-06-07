#' Constructor for H5ClusterRun Objects
#'
#' @description
#' Creates a new \code{H5ClusterRun} object, representing a single run of full
#' voxel-level clustered data from an HDF5 file.
#'
#' This function handles file opening/closing and reads necessary metadata like
#' \code{n_time} from the HDF5 file if not provided explicitly.
#'
#' @details
#' Diagnostic messages (such as when \code{n_time} is inferred from a dataset)
#' are printed only when the option \code{fmristore.verbose} is set to
#' \code{TRUE} via \code{options(fmristore.verbose = TRUE)}.
#'
#' @param file Character string path to the HDF5 file.
#' @param scan_name The character string name of the scan (e.g., "run1").
#' @param mask A \code{LogicalNeuroVol} object for the brain mask.
#' @param clusters A \code{ClusteredNeuroVol} object for cluster assignments.
#' @param n_time (Optional) The number of time points. If \code{NULL} (default), the function
#'   will attempt to read it from the HDF5 file attributes or metadata dataset.
#' @param compress (Optional) Logical indicating compression status (metadata).
#'
#' @return A new \code{H5ClusterRun} object with an open file handle managed by the object.
#'
#' @examples
#' \dontrun{
#' # Create temporary HDF5 file with minimal experiment structure
#' temp_file <- tempfile(fileext = ".h5")
#' exp_file <- fmristore:::create_minimal_h5_for_H5ClusterExperiment(file_path = temp_file)
#'
#' # Create mask and clusters
#' mask <- fmristore:::create_minimal_LogicalNeuroVol(dims = c(5, 5, 4))
#' clusters <- fmristore:::create_minimal_ClusteredNeuroVol(mask_vol = mask, num_clusters = 3)
#'
#' # Create H5ClusterRun object
#' run <- H5ClusterRun(exp_file, scan_name = "Run1_Full", mask = mask, clusters = clusters)
#'
#' # Access properties
#' print(run@n_time)
#' print(dim(mask))
#'
#' # Clean up
#' close(run)
#' unlink(temp_file)
#' }
#'
#' @importFrom methods new is
#' @importFrom hdf5r h5attr
#' @importFrom withr defer
#' @export
H5ClusterRun <- function(file, scan_name,
                         mask, clusters,
                         n_time = NULL, compress = FALSE) {
  # --- 1. Argument Validation (Basic) ---
  if (!((is.character(file) && length(file) == 1 && nzchar(file)) || inherits(file, "H5File"))) {
    stop("[H5ClusterRun] 'file' must be a non-empty character string path or an H5File object.")
  }
  if (!is.character(scan_name) || length(scan_name) != 1 || !nzchar(scan_name)) {
    stop("[H5ClusterRun] 'scan_name' must be a non-empty character string.")
  }
  if (!is(mask, "LogicalNeuroVol")) {
    stop("[H5ClusterRun] 'mask' must be a LogicalNeuroVol object.")
  }
  if (!is(clusters, "ClusteredNeuroVol")) {
    stop("[H5ClusterRun] 'clusters' must be a ClusteredNeuroVol object.")
  }
  # Dimension consistency check (using helper)
  check_same_dims(mask, clusters, dims_to_compare = 1:3,
    msg = "[H5ClusterRun] Dimensions of 'mask' and 'clusters' must match.")

  n_vox <- sum(mask)
  if (!identical(length(clusters@clusters), as.integer(n_vox))) {
    stop(sprintf(
      "[H5ClusterRun] Mismatch: clusters@clusters length (%d) != sum(mask) (%d).",
      length(clusters@clusters), n_vox
    ))
  }

  if (!is.logical(compress) || length(compress) != 1) {
    stop("[H5ClusterRun] 'compress' must be a single logical value.")
  }

  # --- 2. Open HDF5 file ---
  # Use open_h5 with mode 'r' by default. It handles file existence check.
  # Note: open_h5 now returns a list(h5=handle, owns=TRUE/FALSE)
  fh <- open_h5(file, mode = "r")
  h5obj <- fh$h5
  # We *don't* defer closure here. The returned object will own the handle.

  determined_n_time <- n_time # Store original or NULL

  # --- 3. Determine n_time (if NULL) ---
  # This logic is copied/adapted from make_run_full
  if (is.null(determined_n_time)) {
    scan_group_path <- paste0("/scans/", scan_name)
    scan_group <- NULL

    scan_group_exists <- tryCatch(h5obj$exists(scan_group_path), error = function(e) FALSE)

    if (scan_group_exists) {
      tryCatch(
        {
          scan_group <- h5obj[[scan_group_path]]
          on.exit(if (!is.null(scan_group) && inherits(scan_group, "H5Group") && scan_group$is_valid) try(scan_group$close()), add = TRUE)

          if ("n_time" %in% hdf5r::h5attr_names(scan_group)) {
            determined_n_time <- hdf5r::h5attr(scan_group, "n_time")
          } else if (tryCatch(scan_group$exists("metadata/n_time"), error = function(e) FALSE)) {
            md_time_ds <- NULL
            tryCatch({
              md_time_ds <- scan_group[["metadata/n_time"]]
              determined_n_time <- md_time_ds$read()
            }, finally = {
              if (!is.null(md_time_ds) && md_time_ds$is_valid) try(md_time_ds$close())
            })
          } else if (tryCatch(scan_group$exists("clusters"), error = function(e) FALSE)) {
            first_cid <- clusters@clusters[1] # Assumes at least one cluster
            if (!is.null(first_cid) && !is.na(first_cid)) {
              dset_path_cid1 <- sprintf("/scans/%s/clusters/cluster_%d", scan_name, as.integer(first_cid))
              if (tryCatch(h5obj$exists(dset_path_cid1), error = function(e) FALSE)) {
                dset_cid1 <- NULL
                tryCatch({
                  dset_cid1 <- h5obj[[dset_path_cid1]]
                  dims <- dset_cid1$dims
                  if (length(dims) == 2) {
                    determined_n_time <- dims[2]
                    if (isTRUE(getOption("fmristore.verbose"))) {
                      message(sprintf("[H5ClusterRun] Inferred n_time = %d from dataset '%s'.", determined_n_time, dset_path_cid1))
                    }
                  }
                }, finally = {
                  if (!is.null(dset_cid1) && dset_cid1$is_valid) try(dset_cid1$close())
                })
              }
            }
          }
        },
        error = function(e) {
          # Close scan_group if opened before error
          if (!is.null(scan_group) && scan_group$is_valid) try(scan_group$close())
          warning(sprintf("[H5ClusterRun] Error reading n_time metadata for scan '%s': %s. Proceeding without inferred n_time.", scan_name, e$message))
        })
    }

    if (is.null(determined_n_time)) {
      # Close the file handle we opened if we couldn't determine n_time
      if (fh$owns) try(h5obj$close_all(), silent = TRUE)
      stop(sprintf("[H5ClusterRun] Could not determine 'n_time' for scan '%s'. Provide it explicitly or ensure it exists in HDF5 attributes/metadata.", scan_name))
    }
  }

  # --- 4. Final n_time validation ---
  if (!is.numeric(determined_n_time) || length(determined_n_time) != 1 || determined_n_time <= 0 || floor(determined_n_time) != determined_n_time) {
    # Close the file handle before stopping
    if (fh$owns) try(h5obj$close_all(), silent = TRUE)
    stop(sprintf("[H5ClusterRun] Determined 'n_time' (%s) must be a single positive integer.", as.character(determined_n_time)))
  }
  final_n_time <- as.integer(determined_n_time)

  # --- 5. Create the object ---
  # The object takes ownership of the h5obj handle.
  # No defer() needed here because we are returning the object containing the handle.
  # The object's finalizer (if defined) or manual closing should handle it later.
  new_obj <- tryCatch(
    {
      new("H5ClusterRun",
        obj       = h5obj, # Pass the open handle
        scan_name = scan_name,
        mask      = mask,
        clusters  = clusters,
        n_voxels  = as.integer(n_vox),
        n_time    = final_n_time,
        compress  = as.logical(compress)
      )
    },
    error = function(e) {
      # If new() fails, close the handle we opened.
      if (fh$owns) try(h5obj$close_all(), silent = TRUE)
      stop(sprintf("[H5ClusterRun] Failed to create object: %s", e$message))
    })

  # Return the created object, which now manages the H5 handle.
  return(new_obj)
}

#' Constructor for H5ClusterRunSummary Objects
#'
#' @description
#' Creates a new \code{H5ClusterRunSummary} object, representing a single run of
#' summary cluster time-series data from an HDF5 file.
#'
#' This function handles file opening/closing, validates the summary dataset,
#' determines \code{n_time} from the dataset dimensions, and reconciles cluster names/IDs.
#'
#' @param file Character string path to the HDF5 file.
#' @param scan_name The character string name of the scan (e.g., "run1").
#' @param mask A \code{LogicalNeuroVol} object for the brain mask (used for reference/consistency).
#' @param clusters (Optional) A \code{ClusteredNeuroVol} object for cluster assignments. If \code{NULL},
#'   cluster names/IDs must be provided or derivable from the dataset.
#' @param cluster_names (Optional) Character vector of names for the clusters.
#' @param cluster_ids (Optional) Integer vector of IDs for the clusters.
#' @param summary_dset (Optional) The name of the dataset within the run's summary group
#'   (default: "summary_data").
#'
#' @return A new \code{H5ClusterRunSummary} object with an open file handle managed by the object.
#'
#' @examples
#' \dontrun{
#' # Create temporary HDF5 file with minimal experiment structure
#' temp_file <- tempfile(fileext = ".h5")
#' exp_file <- fmristore:::create_minimal_h5_for_H5ClusterExperiment(file_path = temp_file)
#'
#' # Create mask
#' mask <- fmristore:::create_minimal_LogicalNeuroVol(dims = c(5, 5, 4))
#'
#' # Create H5ClusterRunSummary object
#' run_summary <- H5ClusterRunSummary(exp_file, scan_name = "Run2_Summary", mask = mask)
#'
#' # Access properties
#' print(run_summary@cluster_names)
#' print(run_summary@n_time)
#'
#' # Clean up
#' close(run_summary)
#' unlink(temp_file)
#' }
#'
#' @importFrom methods new is
#' @importFrom hdf5r H5D
#' @export
H5ClusterRunSummary <- function(file, scan_name,
                                mask, clusters = NULL,
                                cluster_names = character(), cluster_ids = integer(),
                                summary_dset = "summary_data") {
  # --- 1. Argument Validation (Basic) ---
  if (!((is.character(file) && length(file) == 1 && nzchar(file)) || inherits(file, "H5File"))) {
    stop("[H5ClusterRunSummary] 'file' must be a non-empty character string path or an H5File object.")
  }
  if (!is.character(scan_name) || length(scan_name) != 1 || !nzchar(scan_name)) {
    stop("[H5ClusterRunSummary] 'scan_name' must be a non-empty character string.")
  }
  if (!is(mask, "LogicalNeuroVol")) {
    stop("[H5ClusterRunSummary] 'mask' must be a LogicalNeuroVol object.")
  }
  n_vox <- sum(mask)
  if (!is.null(clusters)) {
    if (!is(clusters, "ClusteredNeuroVol")) {
      stop("[H5ClusterRunSummary] 'clusters' must be a ClusteredNeuroVol object.")
    }
    check_same_dims(mask, clusters, dims_to_compare = 1:3,
      msg = "[H5ClusterRunSummary] Dimensions of 'mask' and provided 'clusters' must match.")
    if (!identical(length(clusters@clusters), as.integer(n_vox))) {
      stop(sprintf(
        "[H5ClusterRunSummary] Mismatch: provided clusters@clusters length (%d) != sum(mask) (%d).",
        length(clusters@clusters), n_vox
      ))
    }
  }
  if (!is.character(summary_dset) || length(summary_dset) != 1 || !nzchar(summary_dset)) {
    stop("[H5ClusterRunSummary] 'summary_dset' must be a non-empty character string.")
  }

  # --- 2. Open HDF5 file ---
  fh <- open_h5(file, mode = "r")
  h5obj <- fh$h5
  # The returned object will own the handle. No defer needed here.

  # --- 3. Validate HDF5 structure before reading data ---
  scan_group_path <- file.path("/scans", scan_name)
  summary_group_path <- file.path(scan_group_path, "clusters_summary")
  dset_path <- file.path(summary_group_path, summary_dset)

  tryCatch(
    {
      assert_h5_path(h5obj, scan_group_path, "scan group")
      assert_h5_path(h5obj, summary_group_path, "summary group")
      assert_h5_path(h5obj, dset_path, "summary dataset")
    },
    error = function(e) {
      if (fh$owns) try(h5obj$close_all(), silent = TRUE)
      stop(paste0("[H5ClusterRunSummary] ", e$message))
    })

  # --- 4. Read dataset dimensions and reconcile cluster info ---
  final_n_time <- NA_integer_
  final_n_clusters <- NA_integer_
  final_cluster_names <- cluster_names
  final_cluster_ids <- cluster_ids

  dataset_dims <- NULL

  info_res <- tryCatch(
    {
      with_h5_dataset(h5obj, dset_path, function(ds) {
        list(
          dims = ds$dims,
          names = if (ds$attr_exists("cluster_names")) ds$attr_read("cluster_names") else character(),
          ids   = if (ds$attr_exists("cluster_ids")) ds$attr_read("cluster_ids") else integer()
        )
      })
    },
    error = function(e) {
      if (fh$owns) try(h5obj$close_all(), silent = TRUE)
      stop(sprintf("[H5ClusterRunSummary] Error processing summary dataset '%s': %s", dset_path, e$message))
    })

  dataset_dims <- info_res$dims
  if (length(dataset_dims) != 2) {
    if (fh$owns) try(h5obj$close_all(), silent = TRUE)
    stop(sprintf("Summary dataset '%s' does not have 2 dimensions (found %d).", dset_path, length(dataset_dims)))
  }
  final_n_time <- dataset_dims[1]
  final_n_clusters <- dataset_dims[2]

  if (length(final_cluster_names) == 0 && length(final_cluster_ids) == 0) {
    final_cluster_names <- info_res$names
    final_cluster_ids   <- info_res$ids

    if (length(final_cluster_names) > 0 && length(final_cluster_names) != final_n_clusters) {
      stop(sprintf("Read cluster_names length (%d) mismatch with dataset columns (%d).", length(final_cluster_names), final_n_clusters))
    }
    if (length(final_cluster_ids) > 0 && length(final_cluster_ids) != final_n_clusters) {
      stop(sprintf("Read cluster_ids length (%d) mismatch with dataset columns (%d).", length(final_cluster_ids), final_n_clusters))
    }

    if (length(final_cluster_names) == 0 && final_n_clusters > 0) {
      final_cluster_names <- paste0("Cluster", seq_len(final_n_clusters))
    }
    if (length(final_cluster_ids) == 0 && final_n_clusters > 0) {
      if (!is.null(clusters) && length(unique(clusters@clusters)) == final_n_clusters) {
        final_cluster_ids <- sort(unique(clusters@clusters))
      } else {
        final_cluster_ids <- seq_len(final_n_clusters)
      }
    }
  }

  if (length(final_cluster_names) != final_n_clusters) {
    warning(sprintf("Final number of cluster names (%d) does not match dataset columns (%d). Resetting to default names.", length(final_cluster_names), final_n_clusters))
    final_cluster_names <- paste0("Col_", seq_len(final_n_clusters))
  }
  if (length(final_cluster_ids) != final_n_clusters) {
    stop(sprintf("Final cluster_ids length (%d) does not match dataset columns (%d).", length(final_cluster_ids), final_n_clusters))
  }
  if (!is.integer(final_cluster_ids)) {
    final_cluster_ids <- as.integer(final_cluster_ids)
  }

  # --- 5. Create the object ---
  new_obj <- tryCatch(
    {
      new("H5ClusterRunSummary",
        obj = h5obj, # Pass the open handle
        scan_name = scan_name,
        mask = mask,
        n_voxels = as.integer(n_vox), # Add n_voxels, inherited from H5ClusteredArray
        clusters = if (is.null(clusters)) new("ClusteredNeuroVol") else clusters,
        n_time = as.integer(final_n_time),
        cluster_names = final_cluster_names,
        cluster_ids = final_cluster_ids,
        summary_dset = summary_dset
      )
    },
    error = function(e) {
      if (fh$owns) try(h5obj$close_all(), silent = TRUE)
      stop(sprintf("[H5ClusterRunSummary] Failed to create H5ClusterRunSummary object: %s", e$message))
    })

  return(new_obj)
}
