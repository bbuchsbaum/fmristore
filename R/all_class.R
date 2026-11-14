#' @import methods
#' @importFrom Matrix Matrix sparseMatrix t
#' @importFrom hdf5r H5File
NULL

# S4 needs to know 'H5File' is an S3 class from hdf5r
setOldClass("H5File")

#' LabeledVolumeSet Class
#'
#' @description
#' A class representing a multi-volume dataset stored in HDF5, where each volume
#' is labeled (similar to a named 4th dimension). We store:
#' \itemize{
#'   \item a 3D mask
#'   \item a set of label strings
#'   \item for each label, data only at the mask's nonzero entries
#' }
#'
#' This extends \code{NeuroVec}, so it is logically a 4D object with dimension
#' \code{[X, Y, Z, #labels]}.
#'
#' @slot obj An \code{H5File} reference (the file handle).
#' @slot mask A \code{LogicalNeuroVol} of shape [X, Y, Z].
#' @slot labels A \code{character} vector for the volume labels.
#' @slot load_env An \code{environment} storing references for lazy loading.
#'
#' @seealso \code{\link[neuroim2]{NeuroVec-class}}
#'
#' @examples
#' \dontrun{
#' if (requireNamespace("neuroim2", quietly = TRUE) &&
#'   requireNamespace("hdf5r", quietly = TRUE) &&
#'   exists("read_labeled_vec", where = "package:fmristore") &&
#'   !is.null(fmristore:::create_minimal_h5_for_LabeledVolumeSet)) {
#'   # Setup: Create a temporary HDF5 file suitable for read_labeled_vec
#'   temp_h5_path <- NULL
#'   lvs <- NULL
#'   tryCatch({
#'     temp_h5_path <- fmristore:::create_minimal_h5_for_LabeledVolumeSet(
#'       vol_dims = c(5L, 4L, 3L),
#'       labels = c("CondA", "CondB"),
#'       num_vols_per_label = 2L
#'     )
#'
#'     # Create a LabeledVolumeSet object using the constructor
#'     lvs <- fmristore::read_labeled_vec(file_name = temp_h5_path)
#'
#'     print(lvs)
#'     print(labels(lvs)) # Access labels
#'
#'     # Access data for a specific label (returns a NeuroVol or similar)
#'     # data_cond_a <- lvs[["CondA"]]
#'     # print(dim(data_cond_a)) # Should be 3D x num_vols_per_label
#'   }, error = function(e) {
#'     message("LabeledVolumeSet example failed: ", e$message)
#'   }, finally = {
#'     # Close HDF5 handle owned by LabeledVolumeSet object
#'     if (!is.null(lvs) && inherits(lvs, "LabeledVolumeSet")) {
#'       try(close(lvs), silent = TRUE)
#'     }
#'     # Cleanup temporary file
#'     if (!is.null(temp_h5_path) && file.exists(temp_h5_path)) {
#'       unlink(temp_h5_path)
#'     }
#'   })
#' } else {
#'   message("Skipping LabeledVolumeSet example: fmristore, neuroim2, hdf5r, or helper not available.")
#' }
#' }
#'
#' @export
setClass("LabeledVolumeSet",
  slots = c(
    obj      = "H5File", # pointer to open file
    mask     = "LogicalNeuroVol", # 3D mask
    labels   = "character", # volume (sub-volume) names
    load_env = "environment" # environment for lazy loading
  ),
  contains = c("NeuroVec")
) # extends NeuroVec

#' H5ParcellatedArray (Virtual Base Class)
#'
#' @description
#' A **virtual** base class for representing parcellated neuroimaging data stored in HDF5.
#' It holds the common elements shared across different representations (e.g., full voxel data, summary data).
#' This class is not intended to be instantiated directly.
#'
#' @slot obj An \code{H5File} object representing the open HDF5 file.
#' @slot mask A \code{LogicalNeuroVol} defining the brain mask (shared across runs).
#' @slot clusters A \code{ClusteredNeuroVol} containing cluster assignments (shared across runs).
#' @slot n_voxels An \code{integer} caching the number of voxels in the mask (sum(mask)).
#'
#' @keywords internal
#' @importFrom hdf5r H5File
#' @importFrom methods validObject
#' @importClassesFrom neuroim2 LogicalNeuroVol ClusteredNeuroVol
#' @family H5Parcellated
#'
#' @examples
#' # H5ParcellatedArray is a virtual class and cannot be directly instantiated.
#' # See its subclasses H5ParcellatedScan and H5ParcellatedScanSummary for examples.
#'
#' # You can check if an object inherits from it:
#' # run_object <- H5ParcellatedScan(...) # Assuming run_object is created
#' # inherits(run_object, "H5ParcellatedArray") # Should return TRUE
#'
#' @export
setClass("H5ParcellatedArray",
  slots = c(
    obj      = "H5File",
    mask     = "LogicalNeuroVol",
    clusters = "ClusteredNeuroVol",
    n_voxels = "integer"
  ),
  prototype = list(
    obj      = NULL,
    mask     = new("LogicalNeuroVol"),
    clusters = new("ClusteredNeuroVol"),
    n_voxels = NA_integer_
  ),
  validity = function(object) {
    errors <- character()
    # Check 1: n_voxels should match sum(mask), if mask is valid
    mask_space_valid <- !is.null(object@mask) &&
      is(object@mask, "LogicalNeuroVol") &&
      validObject(object@mask@space, test = TRUE)

    if (mask_space_valid) {
      expected_nvox <- sum(object@mask)
      if (!is.na(object@n_voxels) && object@n_voxels != expected_nvox) {
        errors <- c(
          errors,
          sprintf(
            "Slot 'n_voxels' (%d) does not match sum(mask) (%d).",
            object@n_voxels, expected_nvox
          )
        )
      }
      # Check 2: length of cluster vector should match n_voxels (which should match sum(mask))
      if (!is.null(object@clusters) && is(object@clusters, "ClusteredNeuroVol") && length(object@clusters@clusters) > 0) {
        if (length(object@clusters@clusters) != expected_nvox) {
          errors <- c(
            errors,
            sprintf(
              "Length of clusters@clusters (%d) does not match sum(mask) (%d).",
              length(object@clusters@clusters), expected_nvox
            )
          )
        }
      }
    }
    # Check 3: Basic type check for H5File handle
    if (!is.null(object@obj) && !inherits(object@obj, "H5File")) {
      errors <- c(errors, "Slot 'obj' must be an H5File object or NULL.")
    }

    if (length(errors) == 0) TRUE else errors
  },
  contains = "VIRTUAL"
)

#' H5ParcellatedScan Class
#'
#' @description
#' Represents a single "run" or "scan" of full voxel-level clustered time-series data
#' stored in an HDF5 file. It inherits common properties (mask, clusters, file handle)
#' from `H5ParcellatedArray` and adds run-specific information.
#'
#' @slot scan_name A \code{character} string identifying the scan (e.g., corresponds to a group under `/scans/`).
#' @slot n_time An \code{integer} specifying the number of time points in this run.
#' @slot compress A \code{logical} indicating if compression was intended or used (metadata).
#' Inherits slots `obj`, `mask`, `clusters`, `n_voxels` from `H5ParcellatedArray`.
#'
#' @seealso \code{\link{H5ParcellatedArray-class}}, \code{\link{H5ParcellatedScanSummary-class}}
#' @family H5Parcellated
#' @export
setClass("H5ParcellatedScan",
  slots = c(
    scan_name = "character",
    n_time    = "integer",
    compress  = "logical" # Primarily for metadata/write path
  ),
  prototype = list(
    scan_name = character(),
    n_time = NA_integer_,
    compress = FALSE
    # Inherits prototype for obj, mask, clusters, n_voxels
  ),
  contains = "H5ParcellatedArray"
)

#' H5ParcellatedScanSummary Class
#'
#' @description
#' Represents a single "run" or "scan" containing only *summary* time-series data
#' for clusters (e.g., mean signal per cluster) stored in an HDF5 file.
#' It inherits common properties from `H5ParcellatedArray` but provides methods suited
#' for accessing summary data (like `as.matrix`) rather than full voxel data.
#'
#' @slot scan_name A \code{character} string identifying the scan.
#' @slot n_time An \code{integer} specifying the number of time points in this run.
#' @slot cluster_names A \code{character} vector providing names for the clusters (columns in summary matrix).
#' @slot cluster_ids An \code{integer} vector of cluster IDs corresponding to `cluster_names`.
#' @slot summary_dset A \code{character} string giving the name of the summary dataset within the run's HDF5 group (e.g., "summary_data").
#' Inherits slots `obj`, `mask`, `n_voxels` from `H5ParcellatedArray`.
#' Note: The `clusters` slot inherited from `H5ParcellatedArray` might be NULL or contain the map for reference, but voxel-level access methods are typically disabled/error.
#'
#' @seealso \code{\link{H5ParcellatedArray-class}}, \code{\link{H5ParcellatedScan-class}}
#' @family H5Parcellated
#' @export
setClass("H5ParcellatedScanSummary",
  slots = c(
    scan_name     = "character",
    n_time        = "integer",
    cluster_names = "character",
    cluster_ids   = "integer",
    summary_dset  = "character"
  ),
  prototype = list(
    scan_name     = character(),
    n_time        = NA_integer_,
    cluster_names = character(),
    cluster_ids   = integer(),
    summary_dset  = "summary_data"
  ),
  contains = "H5ParcellatedArray"
)

#' H5ParcellatedMultiScan Class
#'
#' @description
#' Represents a collection of clustered neuroimaging runs (scans) stored within a single HDF5 file.
#' It acts as a container for `H5ParcellatedScan` and/or `H5ParcellatedScanSummary` objects,
#' along with associated metadata.
#'
#' This class facilitates managing multiple runs that share the same mask and cluster definitions.
#'
#' @slot runs A `list` where each element is an object inheriting from `H5ParcellatedArray`
#'   (typically `H5ParcellatedScan` or `H5ParcellatedScanSummary`).
#' @slot scan_metadata A `list` containing metadata for each scan in the `runs` list.
#' @slot cluster_metadata A `data.frame` containing metadata associated with the clusters.
#'
#' @seealso \code{\link{H5ParcellatedScan-class}}, \code{\link{H5ParcellatedScanSummary-class}}
#' @family H5Cluster
#' @export
setClass("H5ParcellatedMultiScan",
  slots = c(
    runs             = "list",
    scan_metadata    = "list",
    cluster_metadata = "data.frame"
  ),
  prototype = list(
    runs             = list(),
    scan_metadata    = list(),
    cluster_metadata = data.frame()
  ),
  validity = function(object) {
    errors <- character()
    n_runs <- length(object@runs)

    # Check 1: All elements in 'runs' list must inherit from H5ParcellatedArray
    if (n_runs > 0) {
      is_valid_run <- vapply(object@runs, inherits, logical(1),
        what = "H5ParcellatedArray"
      )
      if (!all(is_valid_run)) {
        invalid_indices <- which(!is_valid_run)
        errors <- c(
          errors,
          paste0(
            "Elements at indices ",
            paste(invalid_indices, collapse = ", "),
            " in the 'runs' list do not inherit from H5ParcellatedArray."
          )
        )
        # Stop further checks if basic type is wrong
        return(errors)
      }
    }

    # Check 2: Ensure scan_metadata has same length as runs if not empty
    if (length(object@scan_metadata) > 0 && length(object@scan_metadata) != n_runs) {
      errors <- c(
        errors,
        sprintf(
          "Length of 'scan_metadata' (%d) does not match length of 'runs' list (%d).",
          length(object@scan_metadata), n_runs
        )
      )
    }

    # Check 3: If multiple runs, verify they share the same H5File, mask, and clusters
    if (n_runs > 1) {
      first_run <- object@runs[[1]]
      # Check H5 File (using filename as a robust check)
      first_filename <- tryCatch(first_run@obj$get_filename(), error = function(e) NA_character_)
      if (is.na(first_filename)) {
        errors <- c(errors, "Could not get HDF5 filename from the first run object.")
      }

      for (i in 2:n_runs) {
        current_run <- object@runs[[i]]
        # Check H5 File consistency
        current_filename <- tryCatch(current_run@obj$get_filename(), error = function(e) NA_character_)
        if (is.na(current_filename) || !identical(first_filename, current_filename)) {
          errors <- c(
            errors,
            sprintf(
              "Run %d uses a different HDF5 file ('%s') than Run 1 ('%s').",
              i, current_filename, first_filename
            )
          )
        }
        # Check mask consistency (using identical to check for same object in memory)
        if (!identical(first_run@mask, current_run@mask)) {
          errors <- c(errors, sprintf("Run %d has a different mask object than Run 1.", i))
        }
        # Check clusters consistency (using identical)
        # Need to handle NULL case for H5ParcellatedScanSummary
        if (!identical(first_run@clusters, current_run@clusters)) {
          # Allow comparison if both are NULL (summary runs might have NULL clusters)
          if (!(is.null(first_run@clusters) && is.null(current_run@clusters))) {
            errors <- c(errors, sprintf("Run %d has a different clusters object than Run 1.", i))
          }
        }
      }
    }

    if (length(errors) == 0) TRUE else errors
  }
  # Does not contain H5ParcellatedArray
)

