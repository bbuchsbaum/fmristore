#' @import methods
#' @importFrom hdf5r H5File
NULL

# S4 needs to know 'H5File' is an S3 class from hdf5r
setOldClass("H5File")

#' H5NeuroVol Class
#'
#' @description
#' A class representing a three-dimensional brain image backed by an HDF5 dataset.
#' \code{H5NeuroVol} objects provide efficient storage and access for large-scale
#' brain imaging data in the HDF5 format.
#'
#' @slot h5obj An instance of class \code{H5File} from the \pkg{hdf5r} package,
#'   representing the underlying HDF5 file containing the brain image data.
#'
#' @details
#' The \code{H5NeuroVol} class inherits a \code{space} slot from
#' \code{\link[neuroim2]{NeuroVol-class}}, which specifies the spatial domain
#' (dimensions, orientation, etc.). Data I/O is performed by reading/writing
#' subsets of the HDF5 dataset, allowing efficient handling of large 3D volumes.
#'
#' @section Inheritance:
#' \code{H5NeuroVol} inherits from:
#' \itemize{
#'   \item \code{\link[neuroim2]{NeuroVol-class}}: Base class for 3D brain volumes
#'   \item \code{\link[neuroim2]{ArrayLike3D-class}}: Interface for 3D array-like operations
#' }
#'
#' @seealso
#' \code{\link[neuroim2]{NeuroVol-class}} for the base 3D brain volume class.
#' \code{\link[hdf5r]{H5File}} for details on HDF5 file objects.
#'
#' @examples
#' if (requireNamespace("neuroim2", quietly = TRUE) &&
#'     requireNamespace("hdf5r", quietly = TRUE) &&
#'     exists("H5NeuroVol", where = "package:fmristore") && # Check if constructor is available
#'     !is.null(fmristore:::create_minimal_h5_for_H5NeuroVol)) { # Check helper
#'
#'   # Setup: Create a temporary HDF5 file using a helper
#'   # The helper creates a dataset named "data/elements" by default.
#'   temp_h5_path <- NULL
#'   h5_vol <- NULL
#'   tryCatch({
#'     temp_h5_path <- fmristore:::create_minimal_h5_for_H5NeuroVol(dims = c(5L, 5L, 3L))
#'     
#'     # Create an H5NeuroVol object using the constructor
#'     # The constructor defaults to dataset_name = "data/elements" if not specified
#'     # and if a NeuroSpace is not given, it reads it from /space in the HDF5 file.
#'     h5_vol <- fmristore::H5NeuroVol(file_name = temp_h5_path) 
#'     
#'     print(h5_vol)
#'     
#'     # Access a subset of the data
#'     subset_data <- h5_vol[1:2, 1:2, 1]
#'     print(subset_data)
#'     
#'   }, error = function(e) {
#'     message("H5NeuroVol example failed: ", e$message)
#'   }, finally = {
#'     # Close HDF5 handle owned by H5NeuroVol object
#'     if (!is.null(h5_vol) && inherits(h5_vol, "H5NeuroVol")) {
#'       try(close(h5_vol), silent = TRUE)
#'     }
#'     # Cleanup temporary file
#'     if (!is.null(temp_h5_path) && file.exists(temp_h5_path)) {
#'       unlink(temp_h5_path)
#'     }
#'   })
#' } else {
#'   message("Skipping H5NeuroVol example: fmristore, neuroim2, hdf5r, or helper not available.")
#' }
#'
#' @export
#' @rdname H5NeuroVol-class
#' @importClassesFrom neuroim2 ArrayLike3D
#' @importMethodsFrom neuroim2 [
#' @importMethodsFrom neuroim2 linear_access
setClass("H5NeuroVol",
         slots = c(
           h5obj = "H5File"  # underlying HDF5 file handle
         ),
         prototype = list(
           h5obj = NULL # Default to NULL, needs explicit file handle
         ),
         contains = c("NeuroVol", "ArrayLike3D"))

#' H5NeuroVec Class
#'
#' @description
#' A class representing a four-dimensional brain image backed by an HDF5 file.
#' \code{H5NeuroVec} objects provide efficient storage and access for large-scale
#' 4D neuroimaging data using the HDF5 file format.
#'
#' @slot obj An instance of class \code{H5File} from the \pkg{hdf5r} package,
#'   representing the underlying HDF5 file containing the 4D brain image data.
#'
#' @details
#' \code{H5NeuroVec} inherits a \code{space} slot from
#' \code{\link[neuroim2]{NeuroVec-class}}, defining the 4D dimensions
#' (e.g., \code{x}, \code{y}, \code{z}, \code{time}). Data are stored in an HDF5
#' dataset and accessed on demand.
#'
#' @section Inheritance:
#' \code{H5NeuroVec} inherits from:
#' \itemize{
#'   \item \code{\link[neuroim2]{NeuroVec-class}}: Base class for 4D brain images
#'   \item \code{\link[neuroim2]{ArrayLike4D-class}}: Interface for 4D array-like operations
#' }
#'
#' @seealso
#' \code{\link[neuroim2]{NeuroVec-class}} for the base 4D brain image class.
#' \code{\link{H5NeuroVol-class}} for the 3D counterpart.
#' \code{\link[hdf5r]{H5File}} for details on HDF5 file objects.
#'
#' @examples
#' if (requireNamespace("neuroim2", quietly = TRUE) &&
#'     requireNamespace("hdf5r", quietly = TRUE) &&
#'     exists("H5NeuroVec", where = "package:fmristore") && 
#'     !is.null(fmristore:::create_minimal_h5_for_H5NeuroVec)) { 
#'
#'   # Setup: Create a temporary HDF5 file using a helper
#'   temp_h5_path <- NULL
#'   h5_vec <- NULL
#'   tryCatch({
#'     temp_h5_path <- fmristore:::create_minimal_h5_for_H5NeuroVec(dims = c(5L, 5L, 3L, 4L))
#'     
#'     # Create an H5NeuroVec object using the constructor
#'     # Constructor defaults: dataset_name = "data/elements", space from HDF5 /space group
#'     h5_vec <- fmristore::H5NeuroVec(file_name = temp_h5_path) 
#'     
#'     print(h5_vec)
#'     
#'     # Access a subset of the data
#'     subset_data <- h5_vec[1:2, 1:2, 1, 1:2]
#'     print(subset_data)
#'     
#'   }, error = function(e) {
#'     message("H5NeuroVec example failed: ", e$message)
#'   }, finally = {
#'     # Close HDF5 handle owned by H5NeuroVec object
#'     if (!is.null(h5_vec) && inherits(h5_vec, "H5NeuroVec")) {
#'       try(close(h5_vec), silent = TRUE)
#'     }
#'     # Cleanup temporary file
#'     if (!is.null(temp_h5_path) && file.exists(temp_h5_path)) {
#'       unlink(temp_h5_path)
#'     }
#'   })
#' } else {
#'   message("Skipping H5NeuroVec example: fmristore, neuroim2, hdf5r, or helper not available.")
#' }
#'
#' @export
#' @rdname H5NeuroVec-class
setClass("H5NeuroVec",
         slots = c(
           obj = "H5File"  # underlying HDF5 file handle
         ),
         prototype = list(
           obj = NULL # Default to NULL
         ),
         contains = c("NeuroVec", "ArrayLike4D"))

#' H5Format Class
#'
#' @description
#' This class represents the HDF5 (Hierarchical Data Format version 5) file format.
#' It extends \code{FileFormat-class} with HDF5-specific attributes.
#'
#' @seealso \code{\link[neuroim2]{FileFormat-class}}
#'
#' @keywords internal
#' @noRd
setClass("H5Format",
         contains = c("FileFormat"))

#' H5NeuroVecSource Class
#'
#' @description
#' A class used internally to produce an \code{\linkS4class{H5NeuroVec}} from an HDF5 file.
#'
#' @slot file_name A \code{character} string specifying the name of the HDF5 file.
#'
#' @seealso \code{\link{H5NeuroVec-class}}
#'
#' @keywords internal
#' @noRd
setClass("H5NeuroVecSource",
         representation(file_name = "character"),
         prototype = list(
           file_name = character() # Default to empty character vector
         ))

#' LatentNeuroVec Class
#'
#' @description
#' A class that represents a 4-dimensional neuroimaging array using a latent space
#' decomposition. It stores the data as a set of basis functions (dictionary) and
#' a corresponding set of loadings (coefficients), enabling efficient representation
#' and manipulation of high-dimensional data.
#'
#' @slot basis A \code{Matrix} object where each column represents a basis vector
#'   in the latent space.
#' @slot loadings A \code{Matrix} object (often sparse) containing the coefficients
#'   for each basis vector across the spatial dimensions.
#' @slot offset A \code{numeric} vector representing a constant offset term for
#'   each voxel or spatial location.
#' @slot map A \code{IndexLookupVol} object representing the mapping from basis to loadings.
#' @slot label A \code{character} string representing the label for the latent vector.
#'
#' @details
#' \code{LatentNeuroVec} inherits from \code{\link[neuroim2]{NeuroVec-class}}
#' and \code{\link[neuroim2]{AbstractSparseNeuroVec-class}}. The original 4D data
#' can be reconstructed as:
#' \deqn{data[v,t] = \sum_k \bigl(basis[t,k] \times loadings[v,k]\bigr) + offset[v]}.
#' (Note: `v` indexes voxels within the mask).
#' 
#' **Important Naming Note:**
#' * In this R object: `@basis` stores temporal components (`nTime x k`), `@loadings` stores spatial components (`nVox x k`).
#' * In the HDF5 spec: `/scans/.../embedding` stores temporal (`nTime x k`), `/basis/basis_matrix` stores spatial (`k x nVox`).
#' The I/O functions handle the mapping and transposition.
#'
#' This approach is especially useful for large datasets where storing the full
#' 4D array is expensive.
#'
#' @section Inheritance:
#' \code{LatentNeuroVec} inherits from:
#' \itemize{
#'   \item \code{\link[neuroim2]{NeuroVec-class}}
#'   \item \code{\link[neuroim2]{AbstractSparseNeuroVec-class}}
#' }
#'
#' @seealso
#' \code{\link[neuroim2]{NeuroVec-class}},
#' \code{\link[neuroim2]{AbstractSparseNeuroVec-class}}.
#'
#' @examples
#' if (requireNamespace("neuroim2", quietly = TRUE) && 
#'     requireNamespace("Matrix", quietly = TRUE) &&
#'     !is.null(fmristore:::create_minimal_LatentNeuroVec)) {
#'   
#'   # Create a LatentNeuroVec object using the helper
#'   # The helper creates a mask, basis, and loadings internally.
#'   # It uses new("LatentNeuroVec", ...) after creating constituent parts if not directly calling
#'   # a neuroim2::LatentNeuroVec constructor, or directly calls a constructor.
#'   # Our helper fmristore:::create_minimal_LatentNeuroVec returns a neuroim2::LatentNeuroVec.
#' 
#'   latent_vec <- NULL
#'   tryCatch({
#'     latent_vec <- fmristore:::create_minimal_LatentNeuroVec(
#'       space_dims = c(5L, 5L, 3L), 
#'       n_time = 8L, 
#'       n_comp = 2L
#'     )
#'     
#'     print(latent_vec)
#'     
#'     # Access slots (example)
#'     # print(dim(latent_vec@basis))
#'     # print(dim(latent_vec@loadings))
#'     
#'     # Example of accessing data (reconstruction for a voxel would be more complex)
#'     # This class is more about representation; direct element access is usually via methods.
#'     # For example, a method might be `series(latent_vec, vox_indices = c(1,2,3))`
#'     # For a simple demonstration, we can show its dimensions:
#'     print(dim(latent_vec)) # from NeuroVec inheritance
#'     
#'   }, error = function(e) {
#'     message("LatentNeuroVec example failed: ", e$message)
#'   })
#'   
#' } else {
#'   message("Skipping LatentNeuroVec example: neuroim2, Matrix, or helper not available.")
#' }
#'
#' @export
#' @rdname LatentNeuroVec-class
setClass("LatentNeuroVec",
         slots = c(
           basis = "Matrix",
           loadings = "Matrix",
           offset = "numeric",
           map = "IndexLookupVol",
           label = "character"
         ),
        
         contains = c("NeuroVec", "AbstractSparseNeuroVec"))

#' LatentNeuroVecSource Class
#'
#' @description
#' A class used intprotoernally to produce a \code{\linkS4class{LatentNeuroVec}} instance.
#'
#' @slot file_name A \code{character} string specifying the file name.
#'
#' @seealso \code{\link{LatentNeuroVec-class}}
#'
#' @keywords internal
#' @noRd
setClass("LatentNeuroVecSource",
         representation(file_name="character"),
         prototype = list(
           file_name = character() # Default to empty character vector
         ))


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
#' if (requireNamespace("neuroim2", quietly = TRUE) &&
#'     requireNamespace("hdf5r", quietly = TRUE) &&
#'     exists("read_labeled_vec", where = "package:fmristore") && 
#'     !is.null(fmristore:::create_minimal_h5_for_LabeledVolumeSet)) {
#'
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
#'     
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
#'
#' @export
setClass("LabeledVolumeSet",
         slots = c(
           obj      = "H5File",          # pointer to open file
           mask     = "LogicalNeuroVol", # 3D mask
           labels   = "character",       # volume (sub-volume) names
           load_env = "environment"     # environment for lazy loading
         ),
         contains = c("NeuroVec"))  # extends NeuroVec

#' H5ClusteredArray (Virtual Base Class)
#'
#' @description
#' A **virtual** base class for representing clustered neuroimaging data stored in HDF5.
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
#' @family H5Cluster
#'
#' @examples
#' # H5ClusteredArray is a virtual class and cannot be directly instantiated.
#' # See its subclasses H5ClusterRun and H5ClusterRunSummary for examples.
#' 
#' # You can check if an object inherits from it:
#' # run_object <- H5ClusterRun(...) # Assuming run_object is created
#' # inherits(run_object, "H5ClusteredArray") # Should return TRUE
#'
#' @export
setClass("H5ClusteredArray",
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
                                 validObject(object@mask@space, test=TRUE)
                                
            if (mask_space_valid) {
                expected_nvox <- sum(object@mask)
                if (!is.na(object@n_voxels) && object@n_voxels != expected_nvox) {
                   errors <- c(errors, 
                               sprintf("Slot 'n_voxels' (%d) does not match sum(mask) (%d).", 
                                       object@n_voxels, expected_nvox))
                }
                # Check 2: length of cluster vector should match n_voxels (which should match sum(mask))
                if (!is.null(object@clusters) && is(object@clusters, "ClusteredNeuroVol") && length(object@clusters@clusters) > 0) {
                   if (length(object@clusters@clusters) != expected_nvox) {
                      errors <- c(errors,
                                  sprintf("Length of clusters@clusters (%d) does not match sum(mask) (%d).",
                                          length(object@clusters@clusters), expected_nvox))
                   }
                }
            }
            # Check 3: Basic type check for H5File handle
            if (!is.null(object@obj) && !inherits(object@obj, "H5File")) {
                errors <- c(errors, "Slot 'obj' must be an H5File object or NULL.")
            }
            
            if (length(errors) == 0) TRUE else errors
         },
         contains = "VIRTUAL")

#' H5ClusterRun Class
#'
#' @description
#' Represents a single "run" or "scan" of full voxel-level clustered time-series data
#' stored in an HDF5 file. It inherits common properties (mask, clusters, file handle)
#' from `H5ClusteredArray` and adds run-specific information.
#'
#' @slot scan_name A \code{character} string identifying the scan (e.g., corresponds to a group under `/scans/`).
#' @slot n_time An \code{integer} specifying the number of time points in this run.
#' @slot compress A \code{logical} indicating if compression was intended or used (metadata).
#' Inherits slots `obj`, `mask`, `clusters`, `n_voxels` from `H5ClusteredArray`.
#'
#' @seealso \code{\link{H5ClusteredArray-class}}, \code{\link{H5ClusterRunSummary-class}}
#' @family H5Cluster
#' @export
setClass("H5ClusterRun",
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
         contains = "H5ClusteredArray")

#' H5ClusterRunSummary Class
#'
#' @description
#' Represents a single "run" or "scan" containing only *summary* time-series data
#' for clusters (e.g., mean signal per cluster) stored in an HDF5 file.
#' It inherits common properties from `H5ClusteredArray` but provides methods suited
#' for accessing summary data (like `as.matrix`) rather than full voxel data.
#'
#' @slot scan_name A \code{character} string identifying the scan.
#' @slot n_time An \code{integer} specifying the number of time points in this run.
#' @slot cluster_names A \code{character} vector providing names for the clusters (columns in summary matrix).
#' @slot cluster_ids An \code{integer} vector of cluster IDs corresponding to `cluster_names`.
#' @slot summary_dset A \code{character} string giving the name of the summary dataset within the run's HDF5 group (e.g., "summary_data").
#' Inherits slots `obj`, `mask`, `n_voxels` from `H5ClusteredArray`.
#' Note: The `clusters` slot inherited from `H5ClusteredArray` might be NULL or contain the map for reference, but voxel-level access methods are typically disabled/error.
#'
#' @seealso \code{\link{H5ClusteredArray-class}}, \code{\link{H5ClusterRun-class}}
#' @family H5Cluster
#' @export
setClass("H5ClusterRunSummary",
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
         contains = "H5ClusteredArray")

#' H5ClusterExperiment Class
#'
#' @description
#' Represents a collection of clustered neuroimaging runs (scans) stored within a single HDF5 file.
#' It acts as a container for `H5ClusterRun` and/or `H5ClusterRunSummary` objects,
#' along with associated metadata.
#'
#' This class facilitates managing multiple runs that share the same mask and cluster definitions.
#'
#' @slot runs A `list` where each element is an object inheriting from `H5ClusteredArray`
#'   (typically `H5ClusterRun` or `H5ClusterRunSummary`).
#' @slot scan_metadata A `list` containing metadata for each scan in the `runs` list.
#' @slot cluster_metadata A `data.frame` containing metadata associated with the clusters.
#'
#' @seealso \code{\link{H5ClusterRun-class}}, \code{\link{H5ClusterRunSummary-class}}
#' @family H5Cluster
#' @export
setClass("H5ClusterExperiment",
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

            # Check 1: All elements in 'runs' list must inherit from H5ClusteredArray
            if (n_runs > 0) {
                is_valid_run <- vapply(object@runs, inherits, logical(1),
                                       what = "H5ClusteredArray")
                if (!all(is_valid_run)) {
                    invalid_indices <- which(!is_valid_run)
                    errors <- c(errors,
                                paste0("Elements at indices ",
                                       paste(invalid_indices, collapse=", "),
                                       " in the 'runs' list do not inherit from H5ClusteredArray."))
                    # Stop further checks if basic type is wrong
                    return(errors)
                }
            }

            # Check 2: Ensure scan_metadata has same length as runs if not empty
            if (length(object@scan_metadata) > 0 && length(object@scan_metadata) != n_runs) {
                 errors <- c(errors,
                           sprintf("Length of 'scan_metadata' (%d) does not match length of 'runs' list (%d).",
                                   length(object@scan_metadata), n_runs))
            }

            # Check 3: If multiple runs, verify they share the same H5File, mask, and clusters
            if (n_runs > 1) {
                first_run <- object@runs[[1]]
                # Check H5 File (using filename as a robust check)
                first_filename <- tryCatch(first_run@obj$get_filename(), error=function(e) NA_character_)
                if(is.na(first_filename)) {
                    errors <- c(errors, "Could not get HDF5 filename from the first run object.")
                }

                for (i in 2:n_runs) {
                    current_run <- object@runs[[i]]
                    # Check H5 File consistency
                    current_filename <- tryCatch(current_run@obj$get_filename(), error=function(e) NA_character_)
                    if (is.na(current_filename) || !identical(first_filename, current_filename)) {
                        errors <- c(errors,
                                   sprintf("Run %d uses a different HDF5 file ('%s') than Run 1 ('%s').",
                                           i, current_filename, first_filename))
                    }
                    # Check mask consistency (using identical to check for same object in memory)
                    if (!identical(first_run@mask, current_run@mask)) {
                        errors <- c(errors, sprintf("Run %d has a different mask object than Run 1.", i))
                    }
                    # Check clusters consistency (using identical)
                    # Need to handle NULL case for H5ClusterRunSummary
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
         # Does not contain H5ClusteredArray
)

#' @seealso \code{\link{H5ClusterRun-class}}, \code{\link{H5ClusterRunSummary-class}}
#' @family H5Cluster
#'
#' @examples
#' if (requireNamespace("neuroim2", quietly = TRUE) &&
#'     requireNamespace("hdf5r", quietly = TRUE) &&
#'     exists("H5ClusterExperiment", where = "package:fmristore") &&
#'     !is.null(fmristore:::create_minimal_h5_for_H5ClusterExperiment)) {
#'
#'   temp_exp_file <- NULL
#'   exp <- NULL
#'   tryCatch({
#'     # 1. Create an HDF5 file for an experiment containing runs
#'     temp_exp_file <- fmristore:::create_minimal_h5_for_H5ClusterExperiment()
#'     
#'     # 2. Load the experiment using the constructor
#'     exp <- fmristore::H5ClusterExperiment(file_path = temp_exp_file)
#'     
#'     # 3. Show the experiment object
#'     print(exp)
#'     print(names(runs(exp))) # Show the names of the runs it loaded
#'     
#'   }, error = function(e) {
#'     message("H5ClusterExperiment example failed: ", e$message)
#'   }, finally = {
#'     # Close H5ClusterExperiment handle (closes underlying file)
#'     if (!is.null(exp)) try(close(exp), silent = TRUE)
#'     # Cleanup temporary file
#'     if (!is.null(temp_exp_file) && file.exists(temp_exp_file)) {
#'       unlink(temp_exp_file)
#'     }
#'   })
#' } else {
#'   message("Skipping H5ClusterExperiment example: dependencies or helper not available.")
#' }
#'
#' @export
setClass("H5ClusterExperiment",
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

            # Check 1: All elements in 'runs' list must inherit from H5ClusteredArray
            if (n_runs > 0) {
                is_valid_run <- vapply(object@runs, inherits, logical(1),
                                       what = "H5ClusteredArray")
                if (!all(is_valid_run)) {
                    invalid_indices <- which(!is_valid_run)
                    errors <- c(errors,
                                paste0("Elements at indices ",
                                       paste(invalid_indices, collapse=", "),
                                       " in the 'runs' list do not inherit from H5ClusteredArray."))
                    # Stop further checks if basic type is wrong
                    return(errors)
                }
            }

            # Check 2: Ensure scan_metadata has same length as runs if not empty
            if (length(object@scan_metadata) > 0 && length(object@scan_metadata) != n_runs) {
                 errors <- c(errors,
                           sprintf("Length of 'scan_metadata' (%d) does not match length of 'runs' list (%d).",
                                   length(object@scan_metadata), n_runs))
            }

            # Check 3: If multiple runs, verify they share the same H5File, mask, and clusters
            if (n_runs > 1) {
                first_run <- object@runs[[1]]
                # Check H5 File (using filename as a robust check)
                first_filename <- tryCatch(first_run@obj$get_filename(), error=function(e) NA_character_)
                if(is.na(first_filename)) {
                    errors <- c(errors, "Could not get HDF5 filename from the first run object.")
                }

                for (i in 2:n_runs) {
                    current_run <- object@runs[[i]]
                    # Check H5 File consistency
                    current_filename <- tryCatch(current_run@obj$get_filename(), error=function(e) NA_character_)
                    if (is.na(current_filename) || !identical(first_filename, current_filename)) {
                        errors <- c(errors,
                                   sprintf("Run %d uses a different HDF5 file ('%s') than Run 1 ('%s').",
                                           i, current_filename, first_filename))
                    }
                    # Check mask consistency (using identical to check for same object in memory)
                    if (!identical(first_run@mask, current_run@mask)) {
                        errors <- c(errors, sprintf("Run %d has a different mask object than Run 1.", i))
                    }
                    # Check clusters consistency (using identical)
                    # Need to handle NULL case for H5ClusterRunSummary
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
         # Does not contain H5ClusteredArray
)
