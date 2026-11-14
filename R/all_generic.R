#' Get the number of scans
#'
#' This generic returns the number of scans in an object (e.g. a sequence of fMRI scans).
#'
#' @param x The object from which to retrieve the number of scans
#' @return An integer representing the number of scans
#'
#' @examples
#' \dontrun{
#' if (requireNamespace("neuroim2", quietly = TRUE) &&
#'     requireNamespace("hdf5r", quietly = TRUE) &&
#'     exists("H5ParcellatedMultiScan", where = "package:fmristore") &&
#'     exists("n_scans", where = "package:fmristore") &&
#'     !is.null(fmristore:::create_minimal_h5_for_H5ParcellatedMultiScan)) {
#'
#'   temp_exp_file <- NULL
#'   exp_obj <- NULL
#'   tryCatch({
#'     # Create a minimal H5ParcellatedMultiScan which contains runs (scans)
#'     temp_exp_file <- fmristore:::create_minimal_h5_for_H5ParcellatedMultiScan()
#'     exp_obj <- fmristore::H5ParcellatedMultiScan(file_path = temp_exp_file)
#'
#'     # Get the number of scans
#'     num_scans <- n_scans(exp_obj)
#'     print(num_scans) # Should be 2 based on the helper
#'
#'   }, error = function(e) {
#'     message("n_scans example failed: ", e$message)
#'   }, finally = {
#'     if (!is.null(exp_obj)) try(close(exp_obj), silent = TRUE)
#'     if (!is.null(temp_exp_file) && file.exists(temp_exp_file)) {
#'       unlink(temp_exp_file)
#'     }
#'   })
#' } else {
#'   message("Skipping n_scans example: dependencies or helper not available.")
#' }
#' }
#'
#' @export
setGeneric("n_scans", function(x) standardGeneric("n_scans"))

#' Get the scan names
#'
#' This generic returns a character vector of scan names or labels.
#'
#' @param x The object from which to retrieve the scan names
#' @return A character vector of scan names
#'
#' @examples
#' \dontrun{
#' if (requireNamespace("neuroim2", quietly = TRUE) &&
#'     requireNamespace("hdf5r", quietly = TRUE) &&
#'     exists("H5ParcellatedMultiScan", where = "package:fmristore") &&
#'     exists("scan_names", where = "package:fmristore") &&
#'     !is.null(fmristore:::create_minimal_h5_for_H5ParcellatedMultiScan)) {
#'
#'   temp_exp_file <- NULL
#'   exp_obj <- NULL
#'   tryCatch({
#'     # Create a minimal H5ParcellatedMultiScan
#'     temp_exp_file <- fmristore:::create_minimal_h5_for_H5ParcellatedMultiScan()
#'     exp_obj <- fmristore::H5ParcellatedMultiScan(file_path = temp_exp_file)
#'
#'     # Get the scan names
#'     s_names <- scan_names(exp_obj)
#'     print(s_names) # Should be c("Run1_Full", "Run2_Summary") or similar
#'
#'   }, error = function(e) {
#'     message("scan_names example failed: ", e$message)
#'   }, finally = {
#'     if (!is.null(exp_obj)) try(close(exp_obj), silent = TRUE)
#'     if (!is.null(temp_exp_file) && file.exists(temp_exp_file)) {
#'       unlink(temp_exp_file)
#'     }
#'   })
#' } else {
#'   message("Skipping scan_names example: dependencies or helper not available.")
#' }
#' }
#'
#' @export
setGeneric("scan_names", function(x) standardGeneric("scan_names"))

#' Get scan metadata
#'
#' This generic returns any available metadata associated with each scan.
#'
#' @param x The object from which to retrieve the scan metadata
#' @return A list (or other structure) containing metadata for each scan
#'
#' @examples
#' \dontrun{
#' if (requireNamespace("neuroim2", quietly = TRUE) &&
#'     requireNamespace("hdf5r", quietly = TRUE) &&
#'     exists("H5ParcellatedMultiScan", where = "package:fmristore") &&
#'     exists("scan_metadata", where = "package:fmristore") &&
#'     !is.null(fmristore:::create_minimal_h5_for_H5ParcellatedMultiScan)) {
#'
#'   temp_exp_file <- NULL
#'   exp_obj <- NULL
#'   tryCatch({
#'     # Create a minimal H5ParcellatedMultiScan
#'     temp_exp_file <- fmristore:::create_minimal_h5_for_H5ParcellatedMultiScan()
#'     exp_obj <- fmristore::H5ParcellatedMultiScan(file_path = temp_exp_file)
#'
#'     # Get the scan metadata
#'     s_meta <- scan_metadata(exp_obj)
#'     print(s_meta)
#'     # The helper currently doesn't add rich scan_metadata,
#'     # so this might be an empty list or list of NULLs by default.
#'     # length(s_meta) == n_scans(exp_obj) # This should hold TRUE
#'
#'   }, error = function(e) {
#'     message("scan_metadata example failed: ", e$message)
#'   }, finally = {
#'     if (!is.null(exp_obj)) try(close(exp_obj), silent = TRUE)
#'     if (!is.null(temp_exp_file) && file.exists(temp_exp_file)) {
#'       unlink(temp_exp_file)
#'     }
#'   })
#' } else {
#'   message("Skipping scan_metadata example: dependencies or helper not available.")
#' }
#' }
#'
#' @export
setGeneric("scan_metadata", function(x) standardGeneric("scan_metadata"))

#' Get cluster metadata
#'
#' This generic returns any available metadata associated with clusters in an object.
#'
#' @param x The object from which to retrieve the cluster metadata
#' @return A data frame or other structure containing metadata for each cluster
#'
#' @examples
#' \dontrun{
#' if (requireNamespace("neuroim2", quietly = TRUE) &&
#'     requireNamespace("hdf5r", quietly = TRUE) &&
#'     exists("H5ParcellatedMultiScan", where = "package:fmristore") &&
#'     exists("cluster_metadata", where = "package:fmristore") &&
#'     !is.null(fmristore:::create_minimal_h5_for_H5ParcellatedMultiScan)) {
#'
#'   temp_exp_file <- NULL
#'   exp_obj <- NULL
#'   tryCatch({
#'     # Create a minimal H5ParcellatedMultiScan
#'     temp_exp_file <- fmristore:::create_minimal_h5_for_H5ParcellatedMultiScan()
#'     exp_obj <- fmristore::H5ParcellatedMultiScan(file_path = temp_exp_file)
#'
#'     # Get the cluster metadata
#'     c_meta <- cluster_metadata(exp_obj)
#'     print(c_meta)
#'     # The helper currently doesn't add rich cluster_metadata,
#'     # so this is likely an empty data.frame or one with default cluster names/IDs.
#'
#'   }, error = function(e) {
#'     message("cluster_metadata example failed: ", e$message)
#'   }, finally = {
#'     if (!is.null(exp_obj)) try(close(exp_obj), silent = TRUE)
#'     if (!is.null(temp_exp_file) && file.exists(temp_exp_file)) {
#'       unlink(temp_exp_file)
#'     }
#'   })
#' } else {
#'   message("Skipping cluster_metadata example: dependencies or helper not available.")
#' }
#' }
#'
#' @export
setGeneric("cluster_metadata", function(x) standardGeneric("cluster_metadata"))

#' Get the HDF5 file object
#'
#' This generic returns the HDF5 file object associated with the object.
#'
#' @param x The object from which to retrieve the HDF5 file object
#' @return The HDF5 file object (from package hdf5r)
#'
#' @examples
#' \dontrun{
#' if (!is.null(fmristore:::create_minimal_h5_for_H5ParcellatedMultiScan)) {
#'   temp_exp_file <- NULL
#'   exp_obj <- NULL
#'   tryCatch({
#'     temp_exp_file <- fmristore:::create_minimal_h5_for_H5ParcellatedMultiScan()
#'     exp_obj <- fmristore::H5ParcellatedMultiScan(file_path = temp_exp_file)
#'
#'     # Get the H5File object
#'     h5f <- h5file(exp_obj)
#'     print(h5f)
#'     # if (requireNamespace("hdf5r", quietly = TRUE)) print(h5f$is_valid)
#'
#'   }, error = function(e) {
#'     message("h5file example failed: ", e$message)
#'   }, finally = {
#'     # Closing exp_obj will close the h5file handle it owns
#'     if (!is.null(exp_obj)) try(close(exp_obj), silent = TRUE)
#'     if (!is.null(temp_exp_file) && file.exists(temp_exp_file)) {
#'       unlink(temp_exp_file)
#'     }
#'   })
#' } else {
#'   message("Skipping h5file example: helper not available.")
#' }
#' }
#'
#' @export
setGeneric("h5file", function(x) standardGeneric("h5file"))

#' Read Dataset from HDF5
#'
#' Generic function for reading neuroimaging datasets from HDF5 format
#'
#' @param x The file path or H5File object to read from
#' @param type Optional character string specifying the expected type.
#'   If NULL (default), the type is auto-detected from the file structure.
#' @param ... Additional arguments passed to methods
#' @return The appropriate fmristore object
#' @export
#' @examples
#' \dontrun{
#' # Auto-detect type from HDF5 file
#' obj <- read_dataset("data.h5")
#'
#' # Specify type explicitly
#' exp <- read_dataset("data.h5", type = "H5ParcellatedMultiScan")
#' }
#' @rdname read_dataset
setGeneric("read_dataset", function(x, ...) standardGeneric("read_dataset"))

#' Get the mask volume
#' @param x An object with a mask, like H5ParcellatedMultiScan
#' @param ... Additional arguments
#' @return The mask object (e.g., a LogicalNeuroVol)
#'
#' @examples
#' \dontrun{
#' # For H5ParcellatedMultiScan:
#' if (!is.null(fmristore:::create_minimal_h5_for_H5ParcellatedMultiScan)) {
#'   temp_exp_file <- NULL
#'   exp_obj <- NULL
#'   tryCatch({
#'     temp_exp_file <- fmristore:::create_minimal_h5_for_H5ParcellatedMultiScan()
#'     exp_obj <- fmristore::H5ParcellatedMultiScan(file_path = temp_exp_file)
#'     mask_vol_exp <- mask(exp_obj)
#'     print(mask_vol_exp)
#'   }, error = function(e) {
#'     message("mask example for H5ParcellatedMultiScan failed: ", e$message)
#'   }, finally = {
#'     if (!is.null(exp_obj)) try(close(exp_obj), silent = TRUE)
#'     if (!is.null(temp_exp_file) && file.exists(temp_exp_file)) {
#'       unlink(temp_exp_file)
#'     }
#'   })
#' } else {
#'   message("Skipping mask example for H5ParcellatedMultiScan: helper not available.")
#' }
#' }
#'
#' @name mask
#' @importFrom neuroim2 mask
#' @export
mask

#' Get the cluster map object
#' @param x An object with cluster assignments (e.g., H5ParcellatedMultiScan)
#' @param ... Additional arguments
#' @return The clusters object (e.g., a ClusteredNeuroVol)
#'
#' @examples
#' \dontrun{
#' # For H5ParcellatedMultiScan:
#' if (!is.null(fmristore:::create_minimal_h5_for_H5ParcellatedMultiScan)) {
#'   temp_exp_file <- NULL
#'   exp_obj <- NULL
#'   tryCatch({
#'     temp_exp_file <- fmristore:::create_minimal_h5_for_H5ParcellatedMultiScan(
#'       master_mask_dims = c(4L,4L,3L),
#'       num_master_clusters = 2L
#'     )
#'     exp_obj <- fmristore::H5ParcellatedMultiScan(file_path = temp_exp_file)
#'
#'     # Get the master cluster map from the experiment
#'     cluster_vol <- clusters(exp_obj)
#'     print(cluster_vol)
#'     # if (requireNamespace("neuroim2", quietly=TRUE)) print(is(cluster_vol, "ClusteredNeuroVol"))
#'
#'     # Individual runs also have cluster information, potentially accessible via their own methods
#'     # run1 <- runs(exp_obj)[["Run1_Full"]]
#'     # run1_clusters <- clusters(run1) # Assuming a method for H5ParcellatedScan
#'     # print(run1_clusters)
#'
#'   }, error = function(e) {
#'     message("clusters example for H5ParcellatedMultiScan failed: ", e$message)
#'   }, finally = {
#'     if (!is.null(exp_obj)) try(close(exp_obj), silent = TRUE)
#'     if (!is.null(temp_exp_file) && file.exists(temp_exp_file)) {
#'       unlink(temp_exp_file)
#'     }
#'   })
#' } else {
#'   message("Skipping clusters example for H5ParcellatedMultiScan: helper not available.")
#' }
#' }
#'
#' @export
#' @rdname clusters-methods
setGeneric("clusters", function(x, ...) standardGeneric("clusters"))


# --- Generics for H5ParcellatedMultiScan Helpers ---

#' Concatenate Voxel Time Series Across Runs (Generic)
#'
#' @param experiment The experiment object (typically \code{\link{H5ParcellatedMultiScan-class}}).
#' @param mask_idx Indices of voxels within the mask of the experiment.
#' @param run_indices Optional: A numeric or character vector specifying which runs to include.
#' @param ... Additional arguments for methods.
#' @return A concatenated matrix (typically time x voxels).
#'
#' @examples
#' \dontrun{
#' if (!is.null(fmristore:::create_minimal_h5_for_H5ParcellatedMultiScan)) {
#'   temp_exp_file <- NULL
#'   exp_obj <- NULL
#'   tryCatch({
#'     temp_exp_file <- fmristore:::create_minimal_h5_for_H5ParcellatedMultiScan(
#'       master_mask_dims = c(4L,4L,2L), # Smaller mask for example
#'       n_time_run1 = 5L # Shorter time series for Run1_Full
#'     )
#'     exp_obj <- fmristore::H5ParcellatedMultiScan(file_path = temp_exp_file)
#'
#'     # Get some valid voxel indices from the mask
#'     # The master mask is created by the helper.
#'     # exp_mask <- mask(exp_obj)
#'     # valid_mask_indices <- which(exp_mask@.Data)
#'     # For simplicity, let's assume first few indices if mask is not empty.
#'     # If mask has at least 2 voxels:
#'     # selected_vox_indices <- valid_mask_indices[1:min(2, length(valid_mask_indices))]
#'     # A more robust way for an example without directly loading mask here:
#'     # The H5ParcellatedMultiScan has a mask, and its linear indices are 1:sum(mask_data)
#'     # For a 4x4x2 mask, if 50% are true, there are 16 true voxels. Indices are 1 to 16.
#'     # Let's pick first 2 voxels in the mask's internal indexing (1-based).
#'
#'     # Note: series_concat typically works on runs with full data (H5ParcellatedScan)
#'     # The helper creates "Run1_Full".
#'     # concatenated_series <- series_concat(exp_obj, mask_idx = c(1, 2), 
#'     #                                       run_indices = "Run1_Full")
#'     # print(dim(concatenated_series)) # Should be n_time_run1 x 2
#'     # print(head(concatenated_series))
#'
#'     # For a fully runnable example, need to ensure mask_idx is valid for the created object.
#'     # Since the helper creates a mask, we can try to use mask_idx = 1 (first voxel in mask).
#'     # The `series_concat` method for H5ParcellatedMultiScan handles this.
#'     if (n_voxels(exp_obj) > 0) { # n_voxels from H5ParcellatedArray slot
#'        conc_series <- series_concat(exp_obj, mask_idx = 1, run_indices = "Run1_Full")
#'        print(paste("Dimensions of concatenated series for voxel 1 from Run1_Full:",
#'                    paste(dim(conc_series), collapse="x")))
#'     } else {
#'        message("Skipping series_concat demonstration as experiment mask is empty.")
#'     }
#'
#'   }, error = function(e) {
#'     message("series_concat example failed: ", e$message)
#'   }, finally = {
#'     if (!is.null(exp_obj)) try(close(exp_obj), silent = TRUE)
#'     if (!is.null(temp_exp_file) && file.exists(temp_exp_file)) {
#'       unlink(temp_exp_file)
#'     }
#'   })
#' } else {
#'   message("Skipping series_concat example: helper not available.")
#' }
#' }
#'
#' @export
setGeneric("series_concat",
           function(experiment, mask_idx, run_indices = NULL, ...) standardGeneric("series_concat"))

#' Concatenate Cluster Summary Matrices Across Runs (Generic)
#'
#' @param experiment The experiment object (typically \code{\link{H5ParcellatedMultiScan-class}}).
#' @param run_indices Optional: A numeric or character vector specifying which runs to include.
#'   If NULL, uses all runs that are of summary type or can produce a summary matrix.
#' @param ... Additional arguments for methods.
#' @return A concatenated matrix (typically time x clusters).
#'
#' @examples
#' \dontrun{
#' if (!is.null(fmristore:::create_minimal_h5_for_H5ParcellatedMultiScan)) {
#'   temp_exp_file <- NULL
#'   exp_obj <- NULL
#'   tryCatch({
#'     temp_exp_file <- fmristore:::create_minimal_h5_for_H5ParcellatedMultiScan(
#'       n_time_run2 = 6L # Shorter time series for Run2_Summary
#'     )
#'     exp_obj <- fmristore::H5ParcellatedMultiScan(file_path = temp_exp_file)
#'
#'     # Note: matrix_concat typically works on runs with summary data (H5ParcellatedScanSummary)
#'     # The helper creates "Run2_Summary".
#'     # concatenated_matrix <- matrix_concat(exp_obj, run_indices = "Run2_Summary")
#'     # print(dim(concatenated_matrix))
#'     # Should be n_time_run2 x n_master_clusters (e.g., 6x3 if num_master_clusters is 3)
#'     # print(head(concatenated_matrix))
#'
#'     # The method should correctly identify and use the summary run.
#'     # If only one summary run is present, run_indices can often be omitted.
#'     conc_matrix <- matrix_concat(exp_obj, run_indices = "Run2_Summary")
#'     print(paste("Dimensions of concatenated matrix from Run2_Summary:",
#'                 paste(dim(conc_matrix), collapse="x")))
#'
#'     # Example with all compatible runs (should pick up Run2_Summary)
#'     # conc_matrix_all <- matrix_concat(exp_obj)
#'     # print(paste("Dimensions of concatenated matrix from all compatible runs:",
#'     #             paste(dim(conc_matrix_all), collapse="x")))
#'
#'   }, error = function(e) {
#'     message("matrix_concat example failed: ", e$message)
#'   }, finally = {
#'     if (!is.null(exp_obj)) try(close(exp_obj), silent = TRUE)
#'     if (!is.null(temp_exp_file) && file.exists(temp_exp_file)) {
#'       unlink(temp_exp_file)
#'     }
#'   })
#' } else {
#'   message("Skipping matrix_concat example: helper not available.")
#' }
#' }
#'
#' @export
setGeneric("matrix_concat",
           function(experiment, run_indices = NULL, ...) standardGeneric("matrix_concat"))


#' Get HDF5 Dataset Path (Internal Generic)
#'
#' @description
#' An internal generic function used to determine the HDF5 path to the dataset
#' corresponding to a specific cluster ID within a given run object.
#' Concrete subclasses (like `H5ParcellatedScan`) must implement a method for this generic.
#'
#' @param x An object inheriting from `H5ParcellatedArray`.
#' @param cid The cluster ID (integer).
#' @param ... Additional arguments (not typically used).
#'
#' @return A character string representing the HDF5 dataset path.
#' @keywords internal
setGeneric(".dataset_path",
           function(x, cid, ...) standardGeneric(".dataset_path"))

#' Extract Time Series Data
#'
#' @description
#' A generic function to extract time series data from an object.
#' Specific methods will determine how data is extracted based on the object's class
#' and the indices provided.
#'
#' @param x The object from which to extract series data.
#' @param i An index, typically specifying voxels or elements.
#'          The interpretation of `i` depends on the specific method.
#'          It can be numeric indices, a matrix of coordinates, etc.
#' @param j Optional y-coordinate or further index specification used by some methods.
#' @param k Optional z-coordinate or further index specification used by some methods.
#' @param ... Additional arguments passed to specific methods (e.g., `drop`).
#'
#' @return A matrix or vector containing the time series data. The exact format
#'   (e.g., time x voxels, or a simple vector for a single voxel) depends on the
#'   method and arguments.
#'
#' @seealso Methods for this generic are available for classes like
#'   \code{\link{H5ParcellatedScan}}.
#'
#' @examples
#' # This is a generic function; examples are provided with its specific methods.
#' # For example, see help("series,H5ParcellatedScan-method")
#'
#' @export
#' @rdname series-methods
#' @name series
series <- neuroim2::series

#' Linear Access to Neuroimaging Data (Methods for neuroim2 Generic)
#'
#' @description
#' These methods provide 4D linear access to data for specific `fmristore` classes,
#' implementing the \code{\link[neuroim2:linear_access]{linear_access}} generic function
#' from the \code{neuroim2} package.
#'
#' The `linear_access` generic allows direct access to data elements using a single
#' numeric index that spans the entire 4D space of the object (X, Y, Z, Time).
#' Refer to the documentation for \code{linear_access} in the \code{neuroim2} package
#' for general details about the generic concept.
#'
#' @param x An object for which a `linear_access` method is defined (e.g., `H5ParcellatedScan`).
#' @param i A numeric vector of 4D linear indices.
#' @param ... Additional arguments, not typically used by `fmristore` methods for `linear_access`
#'   but may be relevant for other methods of this generic.
#'
#' @return A numeric vector of values corresponding to the provided linear indices.
#'   The order of values in the returned vector matches the order of indices in `i`.
#'
#' @seealso \code{neuroim2::\link[neuroim2]{linear_access}}, specific methods like
#'   \code{\link{linear_access,H5ParcellatedScan,numeric-method}}.
#'
#' @name linear_access-methods
#' @aliases linear_access
#' @importFrom neuroim2 linear_access
#' @keywords internal
linear_access <- neuroim2::linear_access



#' Convert to Matrix
#'
#' @description
#' These methods convert `fmristore` specific objects to matrices,
#' implementing the S4 generic function \code{\link[base]{as.matrix}} from the `base` package.
#'
#' @param x An object for which a `as.matrix` method is defined.
#' @param ... Additional arguments passed to specific as.matrix methods.
#' @rdname as.matrix-methods
#' @name as.matrix
as.matrix <- neuroim2::as.matrix



#' Convert to Data Frame
#'
#' @description
#' These methods convert `fmristore` specific objects to data frames,
#' implementing the S4 generic function \code{\link[base]{as.data.frame}} from the `base` package.
#'
#' @param x An object for which a `as.data.frame` method is defined.
#' @param row.names A character vector giving the row names for the data frame, or \code{NULL}.
#' @param optional Logical. If \code{TRUE}, setting row names and converting column names 
#'   is optional.
#' @param ... Additional arguments passed to the underlying as.data.frame methods.
#' @rdname as.data.frame-methods
#' @name as.data.frame
as.data.frame <- base::as.data.frame

#' Get Object Names (Methods for base R Generic)
#'
#' @description
#' These methods retrieve the names from `fmristore` specific objects,
#' implementing the S4 generic function \code{\link[base]{names}} from the `base` package.
#'
#' @param x An object for which a `names` method is defined.
#'
#' @return A character vector of names associated with the object.
#'
#' @seealso \code{base::\link[base]{names}}
#'
#' @rdname names-methods
#' @name names
#' @keywords internal
NULL

#' Close an Object
#'
#' @description
#' Close connections or release resources associated with an object.
#' This generic is used for objects that manage external resources like
#' HDF5 file handles.
#'
#' @param con The object to close (e.g., an object containing an HDF5 file handle)
#' @param ... Additional arguments passed to specific methods
#' @return Usually returns NULL invisibly
#'
#' @examples
#' \dontrun{
#' # Example with LabeledVolumeSet
#' if (!is.null(fmristore:::create_minimal_h5_for_LabeledVolumeSet)) {
#'   temp_file <- fmristore:::create_minimal_h5_for_LabeledVolumeSet()
#'   lvs <- read_labeled_vec(temp_file)
#'   # Use the object...
#'   close(lvs) # Close the HDF5 file handle
#'   unlink(temp_file)
#' }
#' }
#'
#' @export
#' @rdname close
setGeneric("close", function(con, ...) standardGeneric("close"),
           useAsDefault = base::close)
