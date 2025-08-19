#' @include all_class.R cluster_array.R
#' @importFrom methods is
NULL

# Contains methods and helper functions for H5ParcellatedMultiScan objects

# Helper function for stricter validation of a user-provided mask against HDF5 content
# Internal to H5ParcellatedMultiScan constructor logic
#
# @param user_mask The mask object provided by the user (already passed basic checks by ensure_mask).
# @param h5_handle An open H5File handle to the experiment file.
# @param expected_cmap_len The expected number of voxels based on /cluster_map length.
# @param master_space The NeuroSpace of the experiment (for dimension information).
# @return Invisible NULL if validation passes, otherwise stops with an error.
# @keywords internal
.validate_user_provided_mask_h5 <- function(user_mask, h5_handle, expected_cmap_len, master_space) {
  if (!is(user_mask, "LogicalNeuroVol")) {
    stop("Internal Error: .validate_user_provided_mask_h5 expects a LogicalNeuroVol for user_mask.")
  }
  if (!is(h5_handle, "H5File") || !h5_handle$is_valid) {
    stop("Internal Error: .validate_user_provided_mask_h5 requires a valid H5File handle.")
  }
  if (!is(master_space, "NeuroSpace")) {
    stop("Internal Error: .validate_user_provided_mask_h5 requires a NeuroSpace object for master_space.")
  }

  message("[H5ParcellatedMultiScan] Performing stricter HDF5-based validation for user-provided mask...")

  # 1. Check consistency with cluster map length (number of active voxels)
  if (sum(user_mask@.Data) != expected_cmap_len) {
    stop(sprintf("Provided mask's TRUE voxel count (%d) does not match length of HDF5 /cluster_map (%d). User mask override rejected.",
      sum(user_mask@.Data), expected_cmap_len))
  }

  # 2. Check voxel pattern consistency against /voxel_coords or /mask in the HDF5 file
  file_voxel_indices <- NULL
  vox_coords_dset <- NULL
  mask_dset_val <- NULL

  # Local on.exit for HDF5 handles opened within this helper
  on.exit(
    {
      if (!is.null(vox_coords_dset) && inherits(vox_coords_dset, "H5D") && vox_coords_dset$is_valid) close_h5_safely(vox_coords_dset)
      if (!is.null(mask_dset_val) && inherits(mask_dset_val, "H5D") && mask_dset_val$is_valid) close_h5_safely(mask_dset_val)
    },
    add = TRUE) # Add to existing on.exit handlers if any in calling scope (though this is top-level helper)

  tryCatch(
    {
      if (h5_handle$exists("/voxel_coords")) {
        vox_coords_dset <- h5_handle[["/voxel_coords"]]
        coords_matrix <- vox_coords_dset$read()
        # Ensure master_space@dim is valid and has 3 dimensions for this calculation
        if (length(master_space@dim) < 3) stop("Master space dimensions are less than 3.")
        # Convert array indices (1-based) to linear indices
        # coords_matrix contains 1-based indices from which(..., arr.ind = TRUE)
        file_voxel_indices <- coords_matrix[, 1] +
          (coords_matrix[, 2] - 1) * master_space@dim[1] +
          (coords_matrix[, 3] - 1) * master_space@dim[1] * master_space@dim[2]
      } else if (h5_handle$exists("/mask")) {
        warning("Dataset '/voxel_coords' not found, validating user mask against HDF5 '/mask' data. Consider adding /voxel_coords for efficiency.")
        mask_dset_val <- h5_handle[["/mask"]]
        file_mask_data <- mask_dset_val$read()
        file_voxel_indices <- which(as.logical(file_mask_data))
      } else {
        stop("Cannot validate user-provided mask pattern: neither '/voxel_coords' nor '/mask' found in HDF5 file.")
      }
    },
    error = function(e) {
      stop(sprintf("Error reading HDF5 '/voxel_coords' or '/mask' for user mask validation: %s", e$message))
    })

  provided_voxel_indices <- which(user_mask@.Data)
  # Convert to same type for comparison (file_voxel_indices might be numeric due to calculations)
  if (!identical(sort(as.integer(provided_voxel_indices)), sort(as.integer(file_voxel_indices)))) {
    # Comparing sorted indices because order might not be guaranteed identical even if sets are same
    stop("User-provided mask's pattern of TRUE voxels does not match the pattern derived from the HDF5 file ('/voxel_coords' or '/mask'). User mask override rejected.")
  }

  message("[H5ParcellatedMultiScan] User-provided mask pattern validation successful.")
  invisible(NULL)
}

#' @describeIn series_concat Concatenate voxel time series for H5ParcellatedMultiScan
#' @export
#' @family H5Parcellated
setMethod("series_concat",
  signature(experiment = "H5ParcellatedMultiScan", mask_idx = "numeric"),
  function(experiment, mask_idx, run_indices = NULL) {

    if (length(experiment@runs) == 0) {
      warning("[series_concat] Experiment contains no runs. Returning empty matrix.")
      return(matrix(numeric(0), nrow = 0, ncol = length(mask_idx %||% 0)))
    }
    # Use the assertion helper
    assert_non_empty_numeric(mask_idx, arg = "mask_idx", fn = "series_concat")

    # Validate run_indices if provided
    if (is.null(run_indices)) {
      run_indices <- seq_along(experiment@runs)
    } else {
      if (!is.numeric(run_indices) || any(run_indices < 1) || any(run_indices > length(experiment@runs))) {
        stop(sprintf("[series_concat] Invalid 'run_indices'. Must be between 1 and %d.", length(experiment@runs)))
      }
      run_indices <- as.integer(run_indices)
    }
    if (length(run_indices) == 0) {
      warning("[series_concat] No runs selected by 'run_indices'. Returning empty matrix.")
      return(matrix(numeric(0), nrow = 0, ncol = length(mask_idx)))
    }

    all_series <- list()
    total_time <- 0
    first_run_voxels <- NULL # To check consistency

    for (idx in run_indices) {
      current_run <- experiment@runs[[idx]]

      if (!is(current_run, "H5ParcellatedScan")) {
        stop(sprintf("[series_concat] Run %d (scan: '%s') is not an H5ParcellatedScan object. Voxel-level series cannot be extracted.",
          idx, current_run@scan_name))
      }

      tryCatch(
        {
          # series() returns [nTime, nVoxels] matrix
          run_series <- series(current_run, i = mask_idx)

          # Basic consistency check on number of columns returned
          if (is.null(first_run_voxels)) {
            first_run_voxels <- ncol(run_series)
            if (first_run_voxels != length(mask_idx)) {
              warning(sprintf("[series_concat] Run %d returned %d voxels, but %d mask indices were requested. Mismatch possible if some indices were outside mask.", idx, first_run_voxels, length(mask_idx)))
            }
          } else if (ncol(run_series) != first_run_voxels) {
            stop(sprintf("[series_concat] Inconsistent number of voxels returned between runs (Run 1: %d, Run %d: %d). Concatenation aborted.", first_run_voxels, idx, ncol(run_series)))
          }

          all_series[[length(all_series) + 1]] <- run_series
          total_time <- total_time + nrow(run_series)

        },
        error = function(e) {
          stop(sprintf("[series_concat] Failed to extract series for run %d (scan: '%s'). Error: %s",
            idx, current_run@scan_name, e$message))
        })
    }

    if (length(all_series) == 0) {
      return(matrix(numeric(0), nrow = 0, ncol = length(mask_idx))) # Should match first_run_voxels if non-null?
    }

    # Use do.call for efficiency
    result_matrix <- do.call(rbind, all_series)

    # Final dimension check
    if (nrow(result_matrix) != total_time || (!is.null(first_run_voxels) && ncol(result_matrix) != first_run_voxels)) {
      warning("[series_concat] Final matrix dimensions seem inconsistent with expectations after rbind. Check results carefully.")
    }

    return(result_matrix)
  })


#' @describeIn matrix_concat Concatenate summary matrices for H5ParcellatedMultiScan
#' @export
#' @family H5Parcellated
setMethod("matrix_concat",
  signature(experiment = "H5ParcellatedMultiScan"),
  function(experiment, run_indices = NULL) {

    if (length(experiment@runs) == 0) {
      warning("[matrix_concat] Experiment contains no runs. Returning empty matrix.")
      return(matrix(numeric(0), nrow = 0, ncol = 0))
    }

    # Validate run_indices if provided
    if (is.null(run_indices)) {
      run_indices <- seq_along(experiment@runs)
    } else {
      if (!is.numeric(run_indices) || any(run_indices < 1) || any(run_indices > length(experiment@runs))) {
        stop(sprintf("[matrix_concat] Invalid 'run_indices'. Must be between 1 and %d.", length(experiment@runs)))
      }
      run_indices <- as.integer(run_indices)
    }
    if (length(run_indices) == 0) {
      warning("[matrix_concat] No runs selected by 'run_indices'. Returning empty matrix.")
      return(matrix(numeric(0), nrow = 0, ncol = 0))
    }

    all_matrices <- list()
    total_time <- 0
    n_clusters <- NULL # Check for consistency

    for (idx in run_indices) {
      current_run <- experiment@runs[[idx]]

      if (!is(current_run, "H5ParcellatedScanSummary")) {
        stop(sprintf("[matrix_concat] Run %d (scan: '%s') is not an H5ParcellatedScanSummary object. Summary matrix cannot be extracted.",
          idx, current_run@scan_name))
      }

      tryCatch(
        {
          run_matrix <- as.matrix(current_run)

          # Check consistency of number of clusters (columns)
          if (is.null(n_clusters)) {
            n_clusters <- ncol(run_matrix)
          } else if (ncol(run_matrix) != n_clusters) {
            stop(sprintf("[matrix_concat] Inconsistent number of clusters found. Run %d has %d columns, expected %d. Concatenation aborted.",
              idx, ncol(run_matrix), n_clusters))
          }

          all_matrices[[length(all_matrices) + 1]] <- run_matrix
          total_time <- total_time + nrow(run_matrix)

        },
        error = function(e) {
          stop(sprintf("[matrix_concat] Failed to extract matrix for run %d (scan: '%s'). Error: %s",
            idx, current_run@scan_name, e$message))
        })
    }

    if (length(all_matrices) == 0) {
      return(matrix(numeric(0), nrow = 0, ncol = (n_clusters %||% 0)))
    }

    result_matrix <- do.call(rbind, all_matrices)

    # Final dimension check
    if (nrow(result_matrix) != total_time || (!is.null(n_clusters) && ncol(result_matrix) != n_clusters)) {
      warning("[matrix_concat] Final matrix dimensions seem inconsistent with expectations after rbind. Check results carefully.")
    }

    return(result_matrix)
  })


#' Constructor for H5ParcellatedMultiScan Objects
#'
#' @description
#' Creates a new `H5ParcellatedMultiScan` object, representing a collection of
#' clustered neuroimaging runs sharing a common HDF5 file, mask, and cluster map.
#'
#' This function handles opening the HDF5 file (if a path is provided),
#' identifying available scans, and creating the appropriate run objects
#' (`H5ParcellatedScan` or `H5ParcellatedScanSummary`) for each scan based on
#' the available data within the HDF5 file structure (following the
#' ClusteredTimeSeriesSpec).
#'
#' @param file Either a character string path to the HDF5 file or an
#'   open \code{H5File} object.
#' @param scan_names (Optional) A character vector specifying which scans under `/scans/`
#'   to include in the experiment. If `NULL` (default), the constructor attempts
#'   to discover all available scan groups under `/scans/`.
#' @param mask (Optional) A `LogicalNeuroVol` object for the brain mask. If `NULL`,
#'   the constructor attempts to load it from `/mask` in the HDF5 file.
#' @param clusters (Optional) A `ClusteredNeuroVol` object for cluster assignments.
#'   If `NULL`, the constructor attempts to load it from `/cluster_map` and
#'   potentially `/voxel_coords` in the HDF5 file.
#' @param scan_metadata (Optional) A list to override or supplement metadata read
#'   from the HDF5 file. If provided, its length should match the number of scans.
#' @param cluster_metadata (Optional) A data.frame to override or supplement
#'   cluster metadata read from `/clusters/cluster_meta` in the HDF5 file.
#' @param summary_preference (Optional) Character string controlling which run type to load.
#'   If \code{NULL} (default), the constructor reads the \code{/scans@summary_only}
#'   attribute to choose a default: \code{"require"} if \code{summary_only=TRUE},
#'   \code{"ignore"} if \code{summary_only=FALSE}, otherwise \code{"prefer"}.
#'   Explicit values can be "require" (only load summary runs, error if missing),
#'   "prefer" (load summary if available, else full), or "ignore" (load full runs only).
#'   This influences whether \code{make_run_summary} or \code{make_run_full} is called.
#'   *Note: This parameter requires careful implementation based on HDF5 content checks.*
#' @param keep_handle_open (Logical) Only relevant if \code{file_source} is a path.
#'   If \code{TRUE} (default), the HDF5 file handle is kept open within the returned
#'   object. If \code{FALSE}, the handle is closed after reading metadata.
#'   *Note:* For most operations, the handle needs to remain open.
#'
#' @return A new \code{H5ParcellatedMultiScan} object.
#'
#' @examples
#' \dontrun{
#' # Create temporary HDF5 file with minimal experiment structure
#' temp_file <- tempfile(fileext = ".h5")
#' exp_file <- fmristore:::create_minimal_h5_for_H5ParcellatedMultiScan(file_path = temp_file)
#'
#' # Create H5ParcellatedMultiScan object
#' experiment <- H5ParcellatedMultiScan(exp_file)
#'
#' # Access scan names
#' print(scan_names(experiment))
#'
#' # Get number of scans
#' print(n_scans(experiment))
#'
#' # Access runs
#' print(names(experiment@runs))
#'
#' # Extract data for specific voxels (first 5 mask indices)
#' voxel_data <- series_concat(experiment, mask_idx = 1:5)
#' print(dim(voxel_data))
#'
#' # Clean up
#' close(experiment)
#' unlink(temp_file)
#' }
#'
#' @importFrom hdf5r H5File list.groups H5A H5D h5attr h5attr_names list.datasets
#' @importFrom methods new is
#' @export
#' @family H5Parcellated
H5ParcellatedMultiScan <- function(file,
                                scan_names = NULL,
                                mask = NULL,
                                clusters = NULL,
                                scan_metadata = NULL,
                                cluster_metadata = NULL,
                                summary_preference = NULL,
                                keep_handle_open = TRUE
) {
  # Rename file_source to file internally for clarity
  file_source <- file
  rm(file) # Remove the old name from the scope

  h5obj <- NULL
  opened_here <- FALSE
  # Track objects opened within this function scope for cleanup
  opened_groups <- list()
  opened_datasets <- list()
  on.exit(
    {
      # Close datasets first
      lapply(opened_datasets, function(ds) if (!is.null(ds) && is(ds, "H5D") && ds$is_valid) close_h5_safely(ds))
      # Then close groups
      lapply(opened_groups, function(grp) if (!is.null(grp) && is(grp, "H5Group") && grp$is_valid) close_h5_safely(grp))
      # Finally, close file if opened here and not keeping open
      if (opened_here && !keep_handle_open && !is.null(h5obj) && h5obj$is_valid) {
        close_h5_safely(h5obj)
      }
      # Finalizer registration is handled later in the function
    },
    add = TRUE)

  # --- 1. Handle File Source and Read Header/Create Space ---
  fh <- open_h5(file_source, mode = "r")
  h5obj <- fh$h5
  # Defer closing the file only if we opened it and the user doesn't want to keep it open.
  # Note: keep_handle_open is TRUE by default, so defer usually won't close it here.
  # The H5ParcellatedMultiScan object's finalizer handles closing if keep_handle_open=TRUE.
  defer(
    {
      if (fh$owns && !keep_handle_open && h5obj$is_valid) {
        message("[H5ParcellatedMultiScan] Closing HDF5 file handle opened by constructor.")
        try(h5obj$close_all(), silent = TRUE)
      }
    },
    envir = parent.frame()) # Defer in the context of the calling function H5ParcellatedMultiScan

  # --- Read Header Info and Create Master NeuroSpace ---
  master_space <- NULL
  hdr_grp <- NULL
  tryCatch(
    {
      hdr_grp_path <- "/header"
      if (!h5obj$exists(hdr_grp_path)) stop("Required '/header' group not found in HDF5 file.")
      hdr_grp <- h5obj[[hdr_grp_path]]
      opened_groups[["header"]] <- hdr_grp

      .rd_hdr <- function(nm) {
        d <- NULL
        val <- NULL
        tryCatch({
          if (hdr_grp$exists(nm)) {
            d <- hdr_grp[[nm]]
            val <- d[]
          }
        }, finally = {
          if (!is.null(d) && d$is_valid) d$close()
        })
        return(val)
      }

      dims   <- .rd_hdr("dim")
      pixdim <- .rd_hdr("pixdim")
      qb     <- .rd_hdr("quatern_b")
      qc <- .rd_hdr("quatern_c")
      qd <- .rd_hdr("quatern_d")
      qx     <- .rd_hdr("qoffset_x")
      qy <- .rd_hdr("qoffset_y")
      qz <- .rd_hdr("qoffset_z")
      qfac   <- .rd_hdr("qfac")

      if (is.null(dims) || length(dims) < 4 || dims[1] != 4) {
        stop("Invalid or missing 'dim' in /header. Must start with 4 and have at least X,Y,Z dims.")
      }
      XYZ_dims <- dims[2:4]

      if (!is.null(pixdim) && length(pixdim) >= 4) {
        spacing_dims <- pixdim[2:4]
      } else {
        warning("Missing or incomplete 'pixdim' in header. Using default spacing (1,1,1).")
        spacing_dims <- c(1, 1, 1)
      }

      qfac_val <- qfac %||% 1.0
      if (!all(sapply(list(qb, qc, qd, qx, qy, qz), function(x) !is.null(x) && is.numeric(x)))) {
        warning("Missing or non-numeric quaternion parameters in header. Using default transform.")
        transform_mat <- diag(4)
        transform_mat[1, 1] <- spacing_dims[1]
        transform_mat[2, 2] <- spacing_dims[2]
        transform_mat[3, 3] <- spacing_dims[3]
        origin_vec <- c(qx %||% 0, qy %||% 0, qz %||% 0)
        transform_mat[1:3, 4] <- origin_vec # NeuroSpace uses corner origin convention
      } else {
        transform_mat <- tryCatch(
          neuroim2::quaternToMatrix(
            quat     = c(qb, qc, qd),
            origin   = c(qx, qy, qz),
            stepSize = spacing_dims,
            qfac     = qfac_val
          ),
          error = function(e) {
            warning("Error calling quaternToMatrix: ", e$message, ". Using default transform.")
            mat_fallback <- diag(4)
            mat_fallback[1, 1] <- spacing_dims[1]
            mat_fallback[2, 2] <- spacing_dims[2]
            mat_fallback[3, 3] <- spacing_dims[3]
            mat_fallback
          }
        )
      }
      master_space <- NeuroSpace(dim = XYZ_dims, spacing = spacing_dims, trans = transform_mat)

    },
    error = function(e) {
      stop(sprintf("[H5ParcellatedMultiScan] Failed to read header and create NeuroSpace: %s", e$message))
    })
  # Header group is automatically closed by on.exit if added to opened_groups

  # --- Auto-load/Validate mask and clusters, using the master_space ---
  # Cache cluster map length for validation
  cmap_len <- NULL
  cmap_dset <- NULL
  tryCatch(
    {
      cmap_path <- "/cluster_map"
      if (!h5obj$exists(cmap_path)) stop("Required dataset '/cluster_map' not found.")
      cmap_dset <- h5obj[[cmap_path]]
      opened_datasets[["cmap"]] <- cmap_dset
      # Get length without reading the whole vector if possible (use H5S V1 API)
      cmap_space <- cmap_dset$get_space()
      cmap_len <- cmap_space$get_simple_extent_dims()$dims
    },
    error = function(e) {
      stop(sprintf("Failed to get length of /cluster_map: %s", e$message))
    })
  # cmap_dset is closed via on.exit
  if (is.null(cmap_len) || length(cmap_len) != 1) {
    stop("/cluster_map must be a 1D dataset.")
  }

  # --- Validate/Load MASK ---
  initial_mask_arg <- mask # Capture the original mask argument from the function call

  # ensure_mask will validate a provided mask or load one from HDF5 if initial_mask_arg is NULL.
  # The 'mask' variable is updated with the result.
  mask <- ensure_mask(initial_mask_arg, h5obj, master_space, path = "/mask")

  # If the user originally provided a mask (initial_mask_arg was not NULL),
  # then perform stricter HDF5-based validation on the 'mask' object
  # (which is now the user-provided mask, assuming it passed ensure_mask's checks).
  if (!is.null(initial_mask_arg)) {
    .validate_user_provided_mask_h5(user_mask = mask,
      h5_handle = h5obj,
      expected_cmap_len = cmap_len,
      master_space = master_space)
  }
  # 'mask' now holds the final, validated mask for subsequent use.

  # --- Validate/Load CLUSTERS ---
  if (is.null(clusters)) {
    message("Clusters argument is NULL, attempting to read from HDF5 file (/cluster_map).")
    # The helper `read_h5_clusters_to_ClusteredNeuroVol` already performs length check against sum(mask)
    # and uses space(mask). Mask is now the final, validated mask object.
    clusters <- read_h5_clusters_to_ClusteredNeuroVol(h5obj, mask, "/cluster_map")
  } else {
    # CLUSTERS PROVIDED - Perform validation
    message("Clusters argument provided. Validating against mask...")
    # 1. Space check
    if (!identical(space(clusters), space(mask))) {
      stop("Provided clusters object's NeuroSpace does not match the mask's NeuroSpace.")
    }
    # 2. Length check
    if (length(clusters@clusters) != sum(mask)) {
      stop(sprintf("Provided clusters vector length (%d) does not match the number of TRUE voxels in the mask (%d).",
        length(clusters@clusters), sum(mask)))
    }
    message("Provided clusters validation successful.")
    # Use the validated provided clusters
  }
  # --- End Auto-load/Validation ---

  # --- 2. Determine Scan Names ---
  scans_group_path <- H5_PATHS$SCANS_GROUP
  if (!h5obj$exists(scans_group_path)) {
    stop(sprintf("[H5ParcellatedMultiScan] Scans group not found at '%s'.", scans_group_path))
  }
  scans_group <- h5obj[[scans_group_path]]
  opened_groups[["scans"]] <- scans_group

  # Read summary_only attribute if present
  summary_only_attr <- NULL
  if ("summary_only" %in% h5attr_names(scans_group)) {
    summary_only_attr <- tryCatch(h5attr(scans_group, "summary_only"),
      error = function(e) NULL)
  }

  if (is.null(summary_preference)) {
    if (isTRUE(summary_only_attr)) {
      summary_preference <- "require"
    } else if (identical(summary_only_attr, FALSE)) {
      summary_preference <- "ignore"
    } else {
      summary_preference <- "prefer"
    }
  }

  # Use names() to list only direct members of /scans
  available_scans <- tryCatch(names(scans_group), error = function(e) character(0))

  if (is.null(scan_names)) {
    if (length(available_scans) == 0) {
      warning("[H5ParcellatedMultiScan] No scan groups found under '/scans'.")
      scan_names <- character(0)
    } else {
      scan_names <- available_scans
    }
  } else {
    if (!is.character(scan_names) || length(scan_names) == 0) {
      stop("[H5ParcellatedMultiScan] Provided 'scan_names' must be a non-empty character vector.")
    }
    missing_scans <- setdiff(scan_names, available_scans)
    if (length(missing_scans) > 0) {
      stop(sprintf("[H5ParcellatedMultiScan] Specified scan names not found under '/scans': %s",
        paste(missing_scans, collapse = ", ")))
    }
    # Keep only requested scans in the specified order
  }

  # --- 5. Create Run Objects and Load Scan Metadata ---
  runs_list <- vector("list", length(scan_names))
  names(runs_list) <- scan_names
  final_scan_metadata <- vector("list", length(scan_names))
  names(final_scan_metadata) <- scan_names

  found_full <- FALSE
  found_summary <- FALSE

  for (sname in scan_names) {
    scan_path <- file.path(scans_group_path, sname)

    # Read class attribute to determine scan type
    scan_class <- NULL
    scan_grp <- NULL
    tryCatch({
      scan_grp <- h5obj[[scan_path]]
      if ("class" %in% h5attr_names(scan_grp)) {
        scan_class <- h5attr(scan_grp, "class")
      }
    }, error = function(e) {
      # Ignore errors, will fall back to data existence checks
    }, finally = {
      if (!is.null(scan_grp) && scan_grp$is_valid) try(scan_grp$close(), silent = TRUE)
    })

    has_full_data <- h5obj$exists(file.path(scan_path, "clusters"))
    has_summary_data <- h5obj$exists(file.path(scan_path, "clusters_summary")) # Check group first
    summary_dset_name <- "summary_data" # Default, could be made configurable
    summary_dset_exists <- FALSE
    summary_group_path <- file.path(scan_path, "clusters_summary") # Define path here

    if (has_summary_data) {
      summary_group <- NULL
      tryCatch({
        # Only need to check existence, H5ParcellatedScanSummary will open group/dataset
        # summary_group <- h5obj[[summary_group_path]]; opened_groups[[paste0(sname,"_summary")]] <- summary_group
        summary_group <- h5obj[[summary_group_path]] # Open just to check dataset existence
        summary_dset_exists <- summary_group$exists(summary_dset_name)
      }, error = function(e) {
        warning(sprintf("Error checking summary dataset for scan %s: %s", sname, e$message))
      }, finally = {
        # Close the group if we opened it just for the check
        if (!is.null(summary_group) && inherits(summary_group, "H5Group") && summary_group$is_valid) {
          close_h5_safely(summary_group)
        }
      })
    }

    found_full <- found_full || has_full_data
    found_summary <- found_summary || summary_dset_exists

    # Determine scan type based on class attribute first, then fall back to data existence
    create_summary <- FALSE
    if (!is.null(scan_class)) {
      # Use class attribute if available
      if (scan_class == "H5ParcellatedScanSummary") {
        create_summary <- TRUE
        if (!summary_dset_exists) {
          stop(sprintf("Scan '%s' has class='H5ParcellatedScanSummary' but summary data not found at %s",
            sname, file.path(summary_group_path, summary_dset_name)))
        }
      } else if (scan_class == "H5ParcellatedScan") {
        create_summary <- FALSE
        if (!has_full_data) {
          stop(sprintf("Scan '%s' has class='H5ParcellatedScan' but full data not found under %s",
            sname, file.path(scan_path, "clusters")))
        }
      } else {
        # Unknown class, fall back to preference logic
        if (summary_preference == "require") {
          if (!summary_dset_exists) stop(sprintf("Summary data required but not found for scan '%s' at %s", sname, file.path(summary_group_path, summary_dset_name)))
          create_summary <- TRUE
        } else if (summary_preference == "prefer") {
          create_summary <- summary_dset_exists
        } else { # ignore
          create_summary <- FALSE
        }
      }
    } else {
      # No class attribute, use preference logic
      if (summary_preference == "require") {
        if (!summary_dset_exists) stop(sprintf("Summary data required but not found for scan '%s' at %s", sname, file.path(summary_group_path, summary_dset_name)))
        create_summary <- TRUE
      } else if (summary_preference == "prefer") {
        create_summary <- summary_dset_exists
      } else { # ignore
        create_summary <- FALSE
      }
    }

    # --- Load Scan-Specific Metadata ---
    current_scan_meta <- list()
    meta_path <- file.path(scan_path, "metadata")
    meta_grp <- NULL
    if (h5obj$exists(meta_path)) {
      tryCatch(
        {
          meta_grp <- h5obj[[meta_path]]
          opened_groups[[paste0(sname, "_meta")]] <- meta_grp
          meta_items <- list.datasets(meta_grp) # Use list.datasets
          for (item_name in meta_items) {
            dset <- NULL
            tryCatch({
              dset <- meta_grp[[item_name]]
              # Handle scalar vs array appropriately
              val <- if (length(dset$dims) == 0 || is.null(dset$dims)) dset[1] else dset[]
              current_scan_meta[[item_name]] <- val
            }, error = function(e_read) {
              warning(sprintf("Failed to read metadata item '%s' for scan '%s': %s", item_name, sname, e_read$message))
            }, finally = {
              if (!is.null(dset) && dset$is_valid) try(dset$close(), silent = TRUE)
            })
          }
        },
        error = function(e_grp) {
          warning(sprintf("Failed to access metadata group for scan '%s': %s", sname, e_grp$message))
        })
      # Meta group closed via on.exit
    }
    final_scan_metadata[[sname]] <- current_scan_meta
    # --- End Scan Metadata Loading ---

    # The new constructors handle cluster name/ID loading internally if needed.
    # scan_specific_cluster_names <- character()
    # scan_specific_cluster_ids <- integer()

    if (create_summary && summary_dset_exists) {
      if (!has_summary_data) stop("Internal logic error: create_summary is TRUE but has_summary_data is FALSE") # Should not happen
      tryCatch(
        {
          # Use the canonical constructor
          runs_list[[sname]] <- H5ParcellatedScanSummary(
            file = h5obj, # MODIFIED: Pass H5File handle
            scan_name = sname,
            mask = mask,       # Pass validated mask
            clusters = clusters, # Pass validated clusters (can be NULL)
            # cluster_names = character(), # Let constructor handle defaults/loading
            # cluster_ids = integer(),     # Let constructor handle defaults/loading
            summary_dset = summary_dset_name
          )
        },
        error = function(e) stop(sprintf("Failed to create H5ParcellatedScanSummary for scan '%s': %s", sname, e$message)))

    } else {

      if (!has_full_data) stop(sprintf("Full voxel data required but not found for scan '%s' under group %s", sname, file.path(scan_path, "clusters")))

      # Get n_time from metadata if available, else pass NULL (constructor will try to find it)
      scan_specific_n_time <- current_scan_meta$n_time %||% NULL

      tryCatch(
        {
          # Use the canonical constructor
          runs_list[[sname]] <- H5ParcellatedScan(
            file = h5obj, # MODIFIED: Pass H5File handle
            scan_name = sname,
            mask = mask,       # Pass validated mask
            clusters = clusters, # Pass validated clusters
            n_time = scan_specific_n_time # Pass potential n_time from metadata
            # compress = ? # Let constructor handle reading compress attr if needed
          )
        },
        error = function(e) stop(sprintf("Failed to create H5ParcellatedScan for scan '%s': %s", sname, e$message)))
    }
  }

  if (!is.null(summary_only_attr)) {
    if (isTRUE(summary_only_attr) && found_full) {
      warning("[H5ParcellatedMultiScan] '/scans@summary_only' is TRUE but full data was found in at least one scan.")
    }
  }

  # --- 6. Load Global Cluster Metadata ---
  if (is.null(cluster_metadata)) {
    clusters_group_path <- "/clusters"
    cluster_meta_group_path <- file.path(clusters_group_path, "cluster_meta")
    global_meta_df <- data.frame()
    if (tryCatch(h5obj$exists(cluster_meta_group_path), error = function(e) FALSE)) {
      meta_obj <- NULL
      tryCatch(
        {
          meta_obj <- h5obj[[cluster_meta_group_path]]
          opened_groups[["global_cluster_meta"]] <- meta_obj

          if (is(meta_obj, "H5Group")) {
            meta_datasets <- list.datasets(meta_obj)
            if (length(meta_datasets) > 0) {
              meta_list <- list()
              for (dname in meta_datasets) {
                dset <- NULL
                tryCatch({
                  dset <- meta_obj[[dname]]
                  meta_list[[dname]] <- if (length(dset$dims) == 0 || is.null(dset$dims)) dset[1] else dset[]
                }, error = function(e_read) {
                  warning(sprintf("Failed to read global cluster metadata dataset '%s': %s", dname, e_read$message))
                }, finally = {
                  if (!is.null(dset) && dset$is_valid) try(dset$close(), silent = TRUE)
                })
              }
              lens <- vapply(meta_list, length, integer(1))
              if (length(unique(lens)) == 1) {
                global_meta_df <- as.data.frame(meta_list)
              } else {
                warning("Global cluster metadata datasets have inconsistent lengths. Returning as a list.")
                global_meta_df <- meta_list
              }
            }
          } else if (is(meta_obj, "H5D")) {
            meta_data <- NULL
            tryCatch(
              {
                meta_data <- meta_obj$read()
              },
              error = function(e_read) {
                warning(sprintf("Failed to read global cluster metadata dataset '%s': %s", cluster_meta_group_path, e_read$message))
              })
            if (!is.null(meta_data)) {
              if (is.data.frame(meta_data)) global_meta_df <- meta_data else global_meta_df <- as.data.frame(meta_data)
            }
            if (!is.null(meta_obj) && meta_obj$is_valid) try(meta_obj$close(), silent = TRUE)
            opened_groups[["global_cluster_meta"]] <- NULL
          } else {
            warning(sprintf("Object '%s' is neither an H5Group nor H5D (class '%s').", cluster_meta_group_path, class(meta_obj)[1]))
            if (!is.null(meta_obj) && inherits(meta_obj, "H5RefClass") && meta_obj$is_valid) {
              try(meta_obj$close(), silent = TRUE)
              opened_groups[["global_cluster_meta"]] <- NULL
            }
          }
        },
        error = function(e_grp) {
          warning(sprintf("Error reading global cluster metadata from '%s': %s", cluster_meta_group_path, e_grp$message))
        })
    }
    cluster_metadata <- global_meta_df
  } else {
    if (!is.data.frame(cluster_metadata)) stop("[H5ParcellatedMultiScan] Provided 'cluster_metadata' must be a data.frame.")
    # User provided metadata overrides loaded
  }

  # --- 7. Override scan metadata if provided ---
  if (!is.null(scan_metadata)) {
    if (!is.list(scan_metadata) || length(scan_metadata) != length(runs_list)) {
      stop(sprintf("[H5ParcellatedMultiScan] Provided 'scan_metadata' must be a list of length %d.", length(runs_list)))
    }
    # Merge/replace carefully. For now, just replace.
    # Ensure names match if replacing
    if (!identical(sort(names(scan_metadata)), sort(names(final_scan_metadata)))) {
      warning("Names in provided 'scan_metadata' do not match the loaded scan names. Using provided names.")
    }
    final_scan_metadata <- scan_metadata
    names(final_scan_metadata) <- names(scan_metadata) # Use names from provided metadata
  }

  # --- 8. Create Experiment Object ---
  exp_obj <- new("H5ParcellatedMultiScan",
    runs = runs_list,
    scan_metadata = final_scan_metadata,
    cluster_metadata = cluster_metadata)

  # --- Add Finalizer if needed (Checklist Item 4.3) ---
  if (opened_here && keep_handle_open) {
    if (isTRUE(getOption("fmristore.verbose"))) {
      message("[H5ParcellatedMultiScan] Registering finalizer to close HDF5 handle when experiment object is garbage collected.")
    }
    reg.finalizer(exp_obj, function(e) {
      # Access the handle from the object's structure
      h5_handle_to_close <- NULL
      if (length(e@runs) > 0 && inherits(e@runs[[1]], "H5ParcellatedArray")) {
         h5_handle_to_close <- e@runs[[1]]@obj
      }
      if (!is.null(h5_handle_to_close) && inherits(h5_handle_to_close, "H5File") && h5_handle_to_close$is_valid) {
         try(h5_handle_to_close$close_all(), silent = TRUE)
      }
    }, onexit = TRUE)
  } else if (opened_here && !keep_handle_open) {
      # Close handle immediately if opened here and not keeping open
      try(h5obj$close_all(), silent = TRUE)
  } # else: handle was passed in open, user manages its lifecycle
  # --- End Finalizer ---

  return(exp_obj)
}


#' Accessor Methods for H5ParcellatedMultiScan

#' Access Slots/Properties using `$`
#'
#' Provides convenient access to shared properties like `mask`, `clusters`,
#' and the underlying HDF5 file object (`obj`) by retrieving them from the
#' first run object stored within the experiment.
#' Also provides access to the experiment's own slots (`runs`, `scan_metadata`, `cluster_metadata`).
#'
#' @param x An `H5ParcellatedMultiScan` object.
#' @param name The name of the property or slot to access (`mask`, `clusters`, `obj`, `runs`, `scan_metadata`, `cluster_metadata`).
#'
#' @return The requested object or value.
#' @export
#' @family H5Parcellated
setMethod("$", "H5ParcellatedMultiScan", function(x, name) {
  # Check for experiment's own slots first
  if (name %in% slotNames(x)) {
    return(slot(x, name))
  }

  # Check for convenient accessors from the first run
  if (name %in% c("mask", "clusters", "obj", "h5file", "n_voxels")) {
    if (length(x@runs) == 0) {
      stop(sprintf("Cannot access property '%s': Experiment contains no runs.", name))
    }
    first_run <- x@runs[[1]]
    if (name == "h5file") name <- "obj" # Map generic to slot name
    if (!hasMethod("slot", signature = class(first_run)) || !(name %in% slotNames(first_run))) {
      stop(sprintf("Cannot access property '%s': Slot not found in the first run object of class '%s'.", name, class(first_run)[1]))
    }
    return(slot(first_run, name))
  }

  # Default behavior or error for unknown names
  # stop(sprintf("'%s' is not a valid slot or accessible property for H5ParcellatedMultiScan", name))
  # Returning NULL might be safer than erroring? Or follow default $ behavior if possible.
  return(NULL)
})

#' Get the HDF5 file object via generic
#' @param x H5ParcellatedMultiScan object
#' @return The HDF5 file object from the first run.
#' @rdname h5file-methods
#' @export
#' @family H5Parcellated
setMethod("h5file", "H5ParcellatedMultiScan", function(x) {
  if (length(x@runs) == 0) {
    stop("Cannot get HDF5 file: Experiment contains no runs.")
  }
  # Assumes all runs share the same handle (enforced by validity check)
  if (is.null(x@runs[[1]]@obj)) {
    stop("H5File object is NULL")
  }
  x@runs[[1]]@obj$filename
})

#' Get the H5File object for H5ParcellatedScan objects
#' @param x H5ParcellatedScan object
#' @return The H5File object handle
#' @rdname h5file-methods
#' @export
#' @family H5Parcellated
setMethod("h5file", "H5ParcellatedScan", function(x) {
  if (is.null(x@obj)) {
    stop("H5File object is NULL")
  }
  x@obj$filename
})

#' Get the H5File object for H5ParcellatedScanSummary objects
#' @param x H5ParcellatedScanSummary object
#' @return The H5File object handle
#' @rdname h5file-methods
#' @export
#' @family H5Parcellated
setMethod("h5file", "H5ParcellatedScanSummary", function(x) {
  if (is.null(x@obj)) {
    stop("H5File object is NULL")
  }
  x@obj$filename
})

#' Get the mask object via generic
#' @param x H5ParcellatedMultiScan object
#' @return The mask object from the first run.
#' @rdname mask-methods
#' @export
#' @family H5Parcellated
setMethod("mask", "H5ParcellatedMultiScan", function(x) {
  if (length(x@runs) == 0) {
    stop("Cannot get mask object: Experiment contains no runs.")
  }
  x@runs[[1]]@mask
})

#' Get the clusters object via generic
#' @param x H5ParcellatedMultiScan object
#' @return The clusters object from the first run.
#' @rdname clusters-methods
#' @export
#' @family H5Parcellated
setMethod("clusters", "H5ParcellatedMultiScan", function(x) {
  if (length(x@runs) == 0) {
    stop("Cannot get clusters object: Experiment contains no runs.")
  }
  x@runs[[1]]@clusters
})

#' Get scan names
#' @param x H5ParcellatedMultiScan object
#' @return Character vector of scan names.
#' @rdname scan_names-methods
#' @export
#' @family H5Parcellated
setMethod("scan_names", "H5ParcellatedMultiScan", function(x) {
  names(x@runs) %||% character(0)
})

#' Get number of scans
#' @param x H5ParcellatedMultiScan object
#' @return Integer number of scans.
#' @rdname n_scans-methods
#' @export
#' @family H5Parcellated
setMethod("n_scans", "H5ParcellatedMultiScan", function(x) {
  length(x@runs)
})

#' Get scan metadata
#' @param x H5ParcellatedMultiScan object
#' @return List of scan metadata.
#' @rdname scan_metadata-methods
#' @export
#' @family H5Parcellated
setMethod("scan_metadata", "H5ParcellatedMultiScan", function(x) {
  x@scan_metadata
})

#' Get cluster metadata
#' @param x H5ParcellatedMultiScan object
#' @return Data frame of cluster metadata.
#' @rdname cluster_metadata-methods
#' @export
#' @family H5Parcellated
setMethod("cluster_metadata", "H5ParcellatedMultiScan", function(x) {
  x@cluster_metadata
})

#' @rdname show-methods
#' @importFrom methods show
setMethod("show", "H5ParcellatedMultiScan", function(object) {
  cat("\nH5ParcellatedMultiScan\n")
  cat("  # runs        :", length(object@runs), "\n")
  # Use %||% which should be defined in io_h5_helpers.R now
  num_clusters <- length(object@cluster_metadata$cluster_id) %||%
    (if (length(object@runs) > 0) length(object@runs[[1]]@cluster_ids) else NA)
  cat("  # clusters    :", num_clusters %||% "Unknown", "\n")
  if (length(object@runs) > 0) {
    cat("  mask dims     :", paste(dim(object@runs[[1]]@mask), collapse = "x"), "\n")
    h5_handle <- object@runs[[1]]@obj
    file_path <- if (is_h5file(h5_handle) && h5_handle$is_valid) h5_handle$get_filename() else "Invalid/Closed Handle"
    cat("  HDF5 file     :", file_path, "\n\n")
  } else {
    cat("  mask dims     : Unknown (no runs loaded)\n")
    cat("  HDF5 file     : Unknown (no runs loaded)\n\n")
  }
})

#' Close the HDF5 file handle
#' @param con H5ParcellatedMultiScan object
#' @param ... Additional arguments (ignored)
#' @return NULL invisibly
#' @rdname close
#' @export
#' @family H5Parcellated
setMethod("close", "H5ParcellatedMultiScan", function(con, ...) {
  # All runs share the same H5File handle, so we only need to close it once
  if (length(con@runs) > 0 && !is.null(con@runs[[1]]@obj)) {
    safe_h5_close(con@runs[[1]]@obj)
  }
  invisible(NULL)
})
