#' Get the number of scans
#'
#' This generic returns the number of scans in an object (e.g. a sequence of fMRI scans).
#'
#' @param x The object from which to retrieve the number of scans
#' @return An integer representing the number of scans
#'
#' @examples
#' if (requireNamespace("neuroim2", quietly = TRUE) &&
#'     requireNamespace("hdf5r", quietly = TRUE) &&
#'     exists("H5ClusterExperiment", where = "package:fmristore") &&
#'     exists("n_scans", where = "package:fmristore") &&
#'     !is.null(fmristore:::create_minimal_h5_for_H5ClusterExperiment)) {
#'
#'   temp_exp_file <- NULL
#'   exp_obj <- NULL
#'   tryCatch({
#'     # Create a minimal H5ClusterExperiment which contains runs (scans)
#'     temp_exp_file <- fmristore:::create_minimal_h5_for_H5ClusterExperiment()
#'     exp_obj <- fmristore::H5ClusterExperiment(file_path = temp_exp_file)
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
#' if (requireNamespace("neuroim2", quietly = TRUE) &&
#'     requireNamespace("hdf5r", quietly = TRUE) &&
#'     exists("H5ClusterExperiment", where = "package:fmristore") &&
#'     exists("scan_names", where = "package:fmristore") &&
#'     !is.null(fmristore:::create_minimal_h5_for_H5ClusterExperiment)) {
#'
#'   temp_exp_file <- NULL
#'   exp_obj <- NULL
#'   tryCatch({
#'     # Create a minimal H5ClusterExperiment
#'     temp_exp_file <- fmristore:::create_minimal_h5_for_H5ClusterExperiment()
#'     exp_obj <- fmristore::H5ClusterExperiment(file_path = temp_exp_file)
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
#' if (requireNamespace("neuroim2", quietly = TRUE) &&
#'     requireNamespace("hdf5r", quietly = TRUE) &&
#'     exists("H5ClusterExperiment", where = "package:fmristore") &&
#'     exists("scan_metadata", where = "package:fmristore") &&
#'     !is.null(fmristore:::create_minimal_h5_for_H5ClusterExperiment)) {
#'
#'   temp_exp_file <- NULL
#'   exp_obj <- NULL
#'   tryCatch({
#'     # Create a minimal H5ClusterExperiment
#'     temp_exp_file <- fmristore:::create_minimal_h5_for_H5ClusterExperiment()
#'     exp_obj <- fmristore::H5ClusterExperiment(file_path = temp_exp_file)
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
#' if (requireNamespace("neuroim2", quietly = TRUE) &&
#'     requireNamespace("hdf5r", quietly = TRUE) &&
#'     exists("H5ClusterExperiment", where = "package:fmristore") &&
#'     exists("cluster_metadata", where = "package:fmristore") &&
#'     !is.null(fmristore:::create_minimal_h5_for_H5ClusterExperiment)) {
#'
#'   temp_exp_file <- NULL
#'   exp_obj <- NULL
#'   tryCatch({
#'     # Create a minimal H5ClusterExperiment
#'     temp_exp_file <- fmristore:::create_minimal_h5_for_H5ClusterExperiment()
#'     exp_obj <- fmristore::H5ClusterExperiment(file_path = temp_exp_file)
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
#' if (!is.null(fmristore:::create_minimal_h5_for_H5ClusterExperiment)) {
#'   temp_exp_file <- NULL
#'   exp_obj <- NULL
#'   tryCatch({
#'     temp_exp_file <- fmristore:::create_minimal_h5_for_H5ClusterExperiment()
#'     exp_obj <- fmristore::H5ClusterExperiment(file_path = temp_exp_file)
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
#'
#' @export
setGeneric("h5file", function(x) standardGeneric("h5file"))

#' @keywords internal
setGeneric("write_dataset", function(x, name, data) standardGeneric("write_dataset"))

#' @keywords internal
setGeneric("read_dataset", function(x, name) standardGeneric("read_dataset"))

#' @keywords internal
setGeneric("has_dataset", function(x, name) standardGeneric("has_dataset"))


# --- Generics for LatentNeuroVec Accessors --- 

#' Get the basis matrix (temporal components)
#' @param x An object, likely a LatentNeuroVec or similar
#' @param ... Additional arguments
#' @return The basis matrix (typically time x components)
#' 
#' @examples
#' # For LatentNeuroVec:
#' if (!is.null(fmristore:::create_minimal_LatentNeuroVec)) {
#'   lnv <- NULL
#'   tryCatch({
#'     lnv <- fmristore:::create_minimal_LatentNeuroVec(
#'       space_dims = c(4L, 4L, 2L), 
#'       n_time = 10L, 
#'       n_comp = 3L
#'     )
#'     b_matrix <- basis(lnv)
#'     print(dim(b_matrix)) # Should be n_time x n_comp (e.g., 10x3)
#'   }, error = function(e) {
#'     message("basis example for LatentNeuroVec failed: ", e$message)
#'   })
#' } else {
#'   message("Skipping basis example for LatentNeuroVec: helper not available.")
#' }
#' 
#' @export
#' @rdname basis-methods
setGeneric("basis", function(x, ...) standardGeneric("basis"))

#' Get the loadings matrix (spatial components)
#' @param x An object, likely a LatentNeuroVec or similar
#' @param ... Additional arguments
#' @return The loadings matrix (typically voxels x components)
#' 
#' @examples
#' # For LatentNeuroVec:
#' if (!is.null(fmristore:::create_minimal_LatentNeuroVec)) {
#'   lnv <- NULL
#'   tryCatch({
#'     # Helper creates a mask, n_mask_voxels determined internally or by arg
#'     lnv <- fmristore:::create_minimal_LatentNeuroVec(
#'       space_dims = c(4L, 4L, 2L), 
#'       n_time = 10L, 
#'       n_comp = 3L
#'     )
#'     l_matrix <- loadings(lnv)
#'     # Dimensions should be n_voxels_in_mask x n_comp
#'     print(dim(l_matrix)) 
#'   }, error = function(e) {
#'     message("loadings example for LatentNeuroVec failed: ", e$message)
#'   })
#' } else {
#'   message("Skipping loadings example for LatentNeuroVec: helper not available.")
#' }
#' 
#' @export
#' @rdname loadings-methods
setGeneric("loadings", function(x, ...) standardGeneric("loadings"))

#' Get the offset vector
#' @param x An object, likely a LatentNeuroVec or similar
#' @param ... Additional arguments
#' @return The offset vector
#' 
#' @examples
#' # For LatentNeuroVec:
#' if (!is.null(fmristore:::create_minimal_LatentNeuroVec)) {
#'   lnv <- NULL
#'   tryCatch({
#'     lnv <- fmristore:::create_minimal_LatentNeuroVec(
#'       space_dims = c(3L, 3L, 2L), 
#'       n_mask_voxels = 4L # Specify a small number of voxels in mask
#'     )
#'     off_vector <- offset(lnv)
#'     print(head(off_vector))
#'     # The neuroim2::LatentNeuroVec constructor (used by helper) defaults to zero offset.
#'     # Length should be n_voxels_in_mask.
#'     if (inherits(lnv, "LatentNeuroVec")) {
#'        # Assuming lnv@mask is a LogicalNeuroVol created by the helper
#'        # The offset vector length should match the number of TRUE voxels in the mask.
#'        # print(length(off_vector) == sum(lnv@mask@.Data))
#'     }
#'   }, error = function(e) {
#'     message("offset example for LatentNeuroVec failed: ", e$message)
#'   })
#' } else {
#'   message("Skipping offset example for LatentNeuroVec: helper not available.")
#' }
#' 
#' @export
#' @rdname offset-methods
setGeneric("offset", function(x, ...) standardGeneric("offset"))

#' Get the mask volume
#' @param x An object with a mask, like LatentNeuroVec or H5ClusterExperiment
#' @param ... Additional arguments
#' @return The mask object (e.g., a LogicalNeuroVol)
#' 
#' @examples
#' # For LatentNeuroVec:
#' if (!is.null(fmristore:::create_minimal_LatentNeuroVec)) {
#'   lnv <- NULL
#'   tryCatch({
#'     lnv <- fmristore:::create_minimal_LatentNeuroVec(space_dims = c(3L,3L,2L))
#'     mask_vol <- mask(lnv)
#'     print(mask_vol)
#'     # if (requireNamespace("neuroim2", quietly=TRUE)) print(is(mask_vol, "LogicalNeuroVol"))
#'   }, error = function(e) {
#'     message("mask example for LatentNeuroVec failed: ", e$message)
#'   })
#' } else {
#'   message("Skipping mask example for LatentNeuroVec: helper not available.")
#' }
#' 
#' # For H5ClusterExperiment:
#' if (!is.null(fmristore:::create_minimal_h5_for_H5ClusterExperiment)) {
#'   temp_exp_file <- NULL
#'   exp_obj <- NULL
#'   tryCatch({
#'     temp_exp_file <- fmristore:::create_minimal_h5_for_H5ClusterExperiment()
#'     exp_obj <- fmristore::H5ClusterExperiment(file_path = temp_exp_file)
#'     mask_vol_exp <- mask(exp_obj)
#'     print(mask_vol_exp)
#'   }, error = function(e) {
#'     message("mask example for H5ClusterExperiment failed: ", e$message)
#'   }, finally = {
#'     if (!is.null(exp_obj)) try(close(exp_obj), silent = TRUE)
#'     if (!is.null(temp_exp_file) && file.exists(temp_exp_file)) {
#'       unlink(temp_exp_file)
#'     }
#'   })
#' } else {
#'   message("Skipping mask example for H5ClusterExperiment: helper not available.")
#' }
#' 
#' @export
#' @rdname mask-methods
setGeneric("mask", function(x, ...) standardGeneric("mask"))

#' Get the index map volume
#' @param x An object with an index map, like LatentNeuroVec
#' @param ... Additional arguments
#' @return The index map object (e.g., an IndexLookupVol from neuroim2)
#' 
#' @examples
#' # For LatentNeuroVec:
#' if (!is.null(fmristore:::create_minimal_LatentNeuroVec)) {
#'   lnv <- NULL
#'   tryCatch({
#'     lnv <- fmristore:::create_minimal_LatentNeuroVec(space_dims = c(3L,3L,2L))
#'     map_vol <- map(lnv)
#'     print(map_vol)
#'     # if (requireNamespace("neuroim2", quietly=TRUE)) print(is(map_vol, "IndexLookupVol"))
#'   }, error = function(e) {
#'     message("map example for LatentNeuroVec failed: ", e$message)
#'   })
#' } else {
#'   message("Skipping map example for LatentNeuroVec: helper not available.")
#' }
#' 
#' @export
#' @rdname map-methods
setGeneric("map", function(x, ...) standardGeneric("map"))

#' Get the cluster map object
#' @param x An object with cluster assignments (e.g., H5ClusterExperiment)
#' @param ... Additional arguments
#' @return The clusters object (e.g., a ClusteredNeuroVol)
#' 
#' @examples
#' # For H5ClusterExperiment:
#' if (!is.null(fmristore:::create_minimal_h5_for_H5ClusterExperiment)) {
#'   temp_exp_file <- NULL
#'   exp_obj <- NULL
#'   tryCatch({
#'     temp_exp_file <- fmristore:::create_minimal_h5_for_H5ClusterExperiment(
#'       master_mask_dims = c(4L,4L,3L), 
#'       num_master_clusters = 2L
#'     )
#'     exp_obj <- fmristore::H5ClusterExperiment(file_path = temp_exp_file)
#'     
#'     # Get the master cluster map from the experiment
#'     cluster_vol <- clusters(exp_obj)
#'     print(cluster_vol)
#'     # if (requireNamespace("neuroim2", quietly=TRUE)) print(is(cluster_vol, "ClusteredNeuroVol"))
#'     
#'     # Individual runs also have cluster information, potentially accessible via their own methods
#'     # run1 <- runs(exp_obj)[["Run1_Full"]]
#'     # run1_clusters <- clusters(run1) # Assuming a method for H5ClusterRun
#'     # print(run1_clusters)
#'     
#'   }, error = function(e) {
#'     message("clusters example for H5ClusterExperiment failed: ", e$message)
#'   }, finally = {
#'     if (!is.null(exp_obj)) try(close(exp_obj), silent = TRUE)
#'     if (!is.null(temp_exp_file) && file.exists(temp_exp_file)) {
#'       unlink(temp_exp_file)
#'     }
#'   })
#' } else {
#'   message("Skipping clusters example for H5ClusterExperiment: helper not available.")
#' }
#' 
#' @export
#' @rdname clusters-methods
setGeneric("clusters", function(x, ...) standardGeneric("clusters"))


# --- Generics for H5ClusterExperiment Helpers --- 

#' Concatenate Voxel Time Series Across Runs (Generic)
#' 
#' @param experiment The experiment object (typically \code{\link{H5ClusterExperiment-class}}).
#' @param mask_idx Indices of voxels within the mask of the experiment.
#' @param run_indices Optional: A numeric or character vector specifying which runs to include.
#' @param ... Additional arguments for methods.
#' @return A concatenated matrix (typically time x voxels).
#' 
#' @examples
#' if (!is.null(fmristore:::create_minimal_h5_for_H5ClusterExperiment)) {
#'   temp_exp_file <- NULL
#'   exp_obj <- NULL
#'   tryCatch({
#'     temp_exp_file <- fmristore:::create_minimal_h5_for_H5ClusterExperiment(
#'       master_mask_dims = c(4L,4L,2L), # Smaller mask for example
#'       n_time_run1 = 5L # Shorter time series for Run1_Full
#'     )
#'     exp_obj <- fmristore::H5ClusterExperiment(file_path = temp_exp_file)
#'     
#'     # Get some valid voxel indices from the mask
#'     # The master mask is created by the helper.
#'     # exp_mask <- mask(exp_obj) 
#'     # valid_mask_indices <- which(exp_mask@.Data) 
#'     # For simplicity, let's assume first few indices if mask is not empty.
#'     # If mask has at least 2 voxels: 
#'     # selected_vox_indices <- valid_mask_indices[1:min(2, length(valid_mask_indices))]
#'     # A more robust way for an example without directly loading mask here:
#'     # The H5ClusterExperiment has a mask, and its linear indices are 1:sum(mask_data)
#'     # For a 4x4x2 mask, if 50% are true, there are 16 true voxels. Indices are 1 to 16.
#'     # Let's pick first 2 voxels in the mask's internal indexing (1-based).
#'     
#'     # Note: series_concat typically works on runs with full data (H5ClusterRun)
#'     # The helper creates "Run1_Full".
#'     # concatenated_series <- series_concat(exp_obj, mask_idx = c(1, 2), run_indices = "Run1_Full")
#'     # print(dim(concatenated_series)) # Should be n_time_run1 x 2
#'     # print(head(concatenated_series))
#'     
#'     # For a fully runnable example, need to ensure mask_idx is valid for the created object.
#'     # Since the helper creates a mask, we can try to use mask_idx = 1 (first voxel in mask).
#'     # The `series_concat` method for H5ClusterExperiment handles this.
#'     if (n_voxels(exp_obj) > 0) { # n_voxels from H5ClusteredArray slot
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
#' 
#' @export
setGeneric("series_concat", 
           function(experiment, mask_idx, run_indices = NULL, ...) standardGeneric("series_concat"))

#' Concatenate Cluster Summary Matrices Across Runs (Generic)
#' 
#' @param experiment The experiment object (typically \code{\link{H5ClusterExperiment-class}}).
#' @param run_indices Optional: A numeric or character vector specifying which runs to include.
#'   If NULL, uses all runs that are of summary type or can produce a summary matrix.
#' @param ... Additional arguments for methods.
#' @return A concatenated matrix (typically time x clusters).
#' 
#' @examples
#' if (!is.null(fmristore:::create_minimal_h5_for_H5ClusterExperiment)) {
#'   temp_exp_file <- NULL
#'   exp_obj <- NULL
#'   tryCatch({
#'     temp_exp_file <- fmristore:::create_minimal_h5_for_H5ClusterExperiment(
#'       n_time_run2 = 6L # Shorter time series for Run2_Summary
#'     )
#'     exp_obj <- fmristore::H5ClusterExperiment(file_path = temp_exp_file)
#'     
#'     # Note: matrix_concat typically works on runs with summary data (H5ClusterRunSummary)
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
#' 
#' @export
setGeneric("matrix_concat", 
           function(experiment, run_indices = NULL, ...) standardGeneric("matrix_concat"))


#' Get HDF5 Dataset Path (Internal Generic)
#'
#' @description
#' An internal generic function used to determine the HDF5 path to the dataset
#' corresponding to a specific cluster ID within a given run object.
#' Concrete subclasses (like `H5ClusterRun`) must implement a method for this generic.
#'
#' @param x An object inheriting from `H5ClusteredArray`.
#' @param cid The cluster ID (integer).
#' @param ... Additional arguments (not typically used).
#'
#' @return A character string representing the HDF5 dataset path.
#' @keywords internal
setGeneric(".dataset_path",
           function(x, cid, ...) standardGeneric(".dataset_path"))


#' Generic function to convert R objects to HDF5 format
#'
#' @description
#' A generic function for converting various types of R objects to HDF5 format,
#' providing a standardized interface for serialization to HDF5.
#'
#' @param object The R object to convert to HDF5 (e.g., a \code{NeuroVol} or \code{NeuroVec}).
#' @param file The path to the HDF5 file to create or modify.
#' @param ... Additional arguments specific to the particular method (see Details).
#' @param data_type For NeuroVec/LatentNeuroVec methods: Storage type (e.g., "FLOAT"). Default "FLOAT"
#' @param chunk_dim For NeuroVec method: Chunk dimensions. Default depends on input dimensions
#' @param compression For all methods: Integer compression level [0..9]. Default varies by method
#' @param mask For LabeledVolumeSet method: The mask to use (LogicalNeuroVol)
#' @param labels For LabeledVolumeSet method: Character vector of labels
#' @param dtype For LabeledVolumeSet method: HDF5 data type for values. Default H5T_NATIVE_DOUBLE
#' @param chunk_size For LabeledVolumeSet/list methods: Integer chunk size for HDF5. Default 1024
#' @param header_values For LabeledVolumeSet method: List of additional header values
#' @param scan_names For list method: Character vector of scan names
#' @param clusters For list method: ClusteredNeuroVol with cluster IDs
#' @param scan_metadata For list method: List of metadata lists, one per scan
#' @param cluster_metadata For list method: Optional data.frame with cluster descriptions
#' @param summary_only For list method: If TRUE, save summary data only
#'
#' @section Methods:
#' \describe{
#'   \item{\code{signature(object = "NeuroVec")}}{
#'     Creates an HDF5 file from a 4D NeuroVec object.
#'     Additional parameters:
#'     \describe{
#'       \item{\code{data_type}}{Storage type (e.g., "FLOAT"). Default "FLOAT"}
#'       \item{\code{chunk_dim}}{Chunk dimensions. Default depends on input dimensions}
#'       \item{\code{compression}}{Integer [0..9], default 6}
#'     }
#'     Returns an H5NeuroVec referencing the new HDF5 file.
#'   }
#'   \item{\code{signature(object = "LatentNeuroVec")}}{
#'     Saves a LatentNeuroVec to an HDF5 file in BasisEmbeddingSpec format.
#'     Additional parameters:
#'     \describe{
#'       \item{\code{data_type}}{Storage type (e.g., "FLOAT"). Default "FLOAT"}
#'       \item{\code{compression}}{Integer [1..9], default 6}
#'     }
#'     Returns an HDF5 file object.
#'   }
#'   \item{\code{signature(object = "LabeledVolume")}}{
#'     Saves a LabeledVolume to an HDF5 file.
#'     Additional parameters:
#'     \describe{
#'       \item{\code{mask}}{The mask to use (LogicalNeuroVol)}
#'       \item{\code{labels}}{Character vector of labels}
#'       \item{\code{compression}}{Integer [0..9], default 4}
#'       \item{\code{dtype}}{HDF5 data type for values. Default H5T_NATIVE_DOUBLE}
#'       \item{\code{chunk_size}}{Integer chunk size for HDF5, default 1024}
#'       \item{\code{header_values}}{List of additional header values}
#'     }
#'     Returns an HDF5 file object.
#'   }
#'   \item{\code{signature(object = "list")}}{
#'     Writes a cluster-based time-series dataset to an HDF5 file.
#'     Additional parameters:
#'     \describe{
#'       \item{\code{scan_names}}{Character vector of scan names}
#'       \item{\code{mask}}{LogicalNeuroVol for 3D geometry}
#'       \item{\code{clusters}}{ClusteredNeuroVol with cluster IDs}
#'       \item{\code{scan_metadata}}{List of metadata lists, one per scan}
#'       \item{\code{cluster_metadata}}{Optional data.frame with cluster descriptions}
#'       \item{\code{summary_only}}{Logical; if TRUE, store only summary data}
#'       \item{\code{compression}}{Integer [0..9], default 4}
#'       \item{\code{chunk_size}}{Chunk dimension for 2D writes, default 1024}
#'     }
#'     Returns an HDF5-backed object representing the clustered dataset (e.g., H5ClusterExperiment).
#'   }
#' }
#'
#' @return An object representing the HDF5 storage, typically of a class
#'   corresponding to the input type (e.g., \code{H5NeuroVol} for \code{NeuroVol} input,
#'   \code{H5NeuroVec} for \code{NeuroVec} input).
#'
#' @examples
#' # Example 1: NeuroVec (DenseNeuroVec) to HDF5
#' # Ensure helper function is available and as_h5 exists
#' # if (!is.null(fmristore:::create_minimal_DenseNeuroVec) &&
#' #     exists("as_h5", where = "package:fmristore")) {
#'   
#' dvec <- fmristore:::create_minimal_DenseNeuroVec(dims = c(3L,3L,2L,4L))
#' temp_h5_file <- tempfile(fileext = ".h5")
#' h5_obj <- NULL
#'   
#' tryCatch({
#'     # Convert DenseNeuroVec to an HDF5 file and get an H5NeuroVec object back
#'     h5_obj <- as_h5(dvec, file = temp_h5_file, 
#'                     data_type = "FLOAT",
#'                     chunk_dim = c(2, 2, 2, 4),
#'                     compression = 4)
#'     
#'     print(h5_obj) # Should be an H5NeuroVec
#'     
#' }, error = function(e) {
#'     message("as_h5 NeuroVec example failed: ", e$message)
#' }, finally = {
#'     if (!is.null(h5_obj)) try(close(h5_obj), silent = TRUE)
#'     if (file.exists(temp_h5_file)) {
#'       unlink(temp_h5_file)
#'     }
#' })
#' # }
#' 
#' # Example 2: LatentNeuroVec to HDF5
#' # if (!is.null(fmristore:::create_minimal_LatentNeuroVec) &&
#' #     exists("as_h5", where = "package:fmristore")) {
#'   
#' lnv <- fmristore:::create_minimal_LatentNeuroVec(
#'     space_dims = c(4L, 4L, 2L),
#'     n_time = 6L,
#'     n_comp = 2L
#' )
#' temp_h5_file_lnv <- tempfile(fileext = ".h5") # Use a different temp file name
#' h5_obj_lnv <- NULL # Use a different object name
#'   
#' tryCatch({
#'     # Convert LatentNeuroVec to HDF5
#'     h5_obj_lnv <- as_h5(lnv, file = temp_h5_file_lnv, compression = 4)
#'     
#'     # File should exist and h5_obj_lnv should be a valid H5File object
#'     if (file.exists(temp_h5_file_lnv)) {
#'       print("LatentNeuroVec HDF5 file created successfully: ", temp_h5_file_lnv)
#'     }
#'     
#' }, error = function(e) {
#'     message("as_h5 LatentNeuroVec example failed: ", e$message)
#' }, finally = {
#'     if (!is.null(h5_obj_lnv) && inherits(h5_obj_lnv, "H5File") && h5_obj_lnv$is_valid) {
#'       try(h5_obj_lnv$close_all(), silent = TRUE)
#'     }
#'     if (file.exists(temp_h5_file_lnv)) {
#'       unlink(temp_h5_file_lnv)
#'     }
#' })
#' # }
#'
#' # Example 3: Clustered dataset (list) to HDF5
#' # (This requires more complex setup, so we use a simplified theoretical example)
#' \donttest{
#'   # Removed: if (requireNamespace("neuroim2", quietly = TRUE) && ... 
#'   
#'   # In practice, you would:
#'   # 1. Create a list of NeuroVec/DenseNeuroVec objects (scan data)
#'   # 2. Create a LogicalNeuroVol mask
#'   # 3. Create a ClusteredNeuroVol defining clusters
#'   # 4. Define scan names and metadata
#'   # 5. Call as_h5() with these components
#'   
#'   # Note: This example is simplified for documentation and won't execute
#'   # as the actual clustered dataset structure is more complex
#'   
#'   message("Example usage for clustered dataset (as_h5 for list method):")
#'   message("  # vecs_list <- ... list of NeuroVec objects ...")
#'   message("  # mask_obj <- fmristore:::create_minimal_LogicalNeuroVol(...)")
#'   message("  # clusters_obj <- fmristore:::create_minimal_ClusteredNeuroVol(...)")
#'   message("  # scan_meta <- list(list(TR=2), list(TR=2))")
#'   message("  # h5_clust_obj <- as_h5(vecs_list, file = tempfile(fileext=\".h5\"), ")
#'   message("  #                       scan_names = c('run1', 'run2'), ")
#'   message("  #                       mask = mask_obj, clusters = clusters_obj, ")
#'   message("  #                       scan_metadata = scan_meta)")
#'   message("  # print(h5_clust_obj)")
#'   message("  # if (!is.null(h5_clust_obj)) try(close(h5_clust_obj), silent=TRUE)")
#'   message("  # if (file.exists(attr(h5_clust_obj, \"filepath\"))) unlink(attr(h5_clust_obj, \"filepath\"))")
#' 
#' }
#'
#' @export
#' @rdname as_h5-methods
setGeneric("as_h5", function(object, file, ...) {
  standardGeneric("as_h5")
})

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
#'   \code{\link{H5ClusterRun}}.
#'
#' @examples
#' # This is a generic function; examples are provided with its specific methods.
#' # For example, see help("series,H5ClusterRun-method")
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
#' @param x An object for which a `linear_access` method is defined (e.g., `H5ClusterRun`).
#' @param i A numeric vector of 4D linear indices.
#' @param ... Additional arguments, not typically used by `fmristore` methods for `linear_access`
#'   but may be relevant for other methods of this generic.
#'
#' @return A numeric vector of values corresponding to the provided linear indices.
#'   The order of values in the returned vector matches the order of indices in `i`.
#'
#' @seealso \code{neuroim2::\link[neuroim2]{linear_access}}, specific methods like
#'   \code{\link{linear_access,H5ClusterRun,numeric-method}}.
#'
#' @name linear_access-methods
#' @aliases linear_access
#' @importFrom neuroim2 linear_access
#' @keywords internal
linear_access <- neuroim2::linear_access



#' Get Dimensions of an Object (Methods for base R Generic)
#'
#' @description
#' These methods retrieve the dimensions of `fmristore` specific objects,
#' implementing the S4 generic function \code{\link[base]{dim}} from the `base` package.
#'
#' For objects that represent neuroimaging data, this typically returns a numeric
#' vector indicating the size of each dimension (e.g., X, Y, Z, Time).
#'
#' @param x An object for which a `dim` method is defined.
#'
#' @return A numeric vector of dimensions. The length and interpretation of the
#'   vector depend on the specific class of `x`.
#'
#' @seealso \code{base::\link[base]{dim}}, specific methods like
#'   \code{dim,H5ClusterRun-method} (if defined and exported).
#'
#' @rdname dim-methods
#' @name dim
#' @aliases dim
#' @keywords internal
NULL


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
#' @param optional Logical. If \code{TRUE}, setting row names and converting column names is optional.
#' @param ... Additional arguments passed to the underlying as.data.frame methods.
#' @rdname as.data.frame-methods
#' @name as.data.frame
as.data.frame <- base::as.data.frame