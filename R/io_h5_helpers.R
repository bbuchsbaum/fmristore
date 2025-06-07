# Helper functions for HDF5 I/O

#' Check if an object is an H5File handle
#' @param x Object to check
#' @return Logical TRUE if inherits from H5File, FALSE otherwise.
#' @keywords internal
is_h5file <- function(x) inherits(x, "H5File")



#' Read /mask into a LogicalNeuroVol
#' @param h5 An open `H5File` handle.
#' @param dset The path to the mask dataset (default: "/mask").
#' @param ref_space A `NeuroSpace` object representing the expected space of the mask.
#' @return A `LogicalNeuroVol` object corresponding to the mask dataset.
#' @keywords internal
#' @importFrom neuroim2 LogicalNeuroVol NeuroSpace space
#' @importFrom assertthat assert_that
read_h5_mask_to_LogicalNeuroVol <- function(h5, dset = "/mask", ref_space) {
  assert_that(inherits(ref_space, "NeuroSpace"), msg = "read_h5_mask_to_LogicalNeuroVol: ref_space must be a NeuroSpace object")
  assert_that(h5$exists(dset), msg = sprintf("read_h5_mask_to_LogicalNeuroVol: mask dataset not found at %s", dset))

  mask_dset <- h5[[dset]]
  mask_data_raw <- NULL
  dims <- NULL
  tryCatch({
    # Get dimensions FIRST from the dataset object
    dims <- mask_dset$dims
    if (is.null(dims) || length(dims) != 3) { # Expecting 3D mask
      stop(sprintf("Dataset '%s' does not have expected 3 dimensions. Found dimensions: %s",
        dset, paste(dims, collapse = "x")))
    }

    # Validate dataset dimensions against the reference space
    ref_dims <- dim(ref_space)
    if (!identical(dims, ref_dims)) {
      stop(sprintf(
        "Mask dimensions in HDF5 (%s) do not match reference space dimensions (%s) for dataset '%s'",
        paste(dims, collapse = "x"),
        paste(ref_dims, collapse = "x"),
        dset
      ))
    }

    # Now read the data (might be flattened)
    mask_data_raw <- mask_dset$read()
  }, finally = {
    if (!is.null(mask_dset) && mask_dset$is_valid) try(mask_dset$close(), silent = TRUE)
  })

  if (is.null(mask_data_raw)) stop("Failed to read mask data from ", dset)

  # Reshape potentially flattened data using stored dimensions, coerce to logical
  mask_data_logical <- array(as.logical(mask_data_raw), dim = dims)

  # Construct and return LogicalNeuroVol using the reference space
  LogicalNeuroVol(mask_data_logical, ref_space)
}

#' Read /cluster_map (+ coords) into ClusteredNeuroVol
#' @keywords internal
#' @importFrom neuroim2 ClusteredNeuroVol space
read_h5_clusters_to_ClusteredNeuroVol <- function(h5, mask,
                                                  map_dset = "/cluster_map") {


  stopifnot(h5$exists(map_dset))
  cmap_dset <- h5[[map_dset]]
  vec <- NULL
  tryCatch({
    # Use $read() here too for consistency, although it's expected to be 1D
    vec <- as.integer(cmap_dset$read())
  }, finally = {
    if (!is.null(cmap_dset) && cmap_dset$is_valid) try(cmap_dset$close(), silent = TRUE)
  })
  if (is.null(vec)) stop("Failed to read cluster map from ", map_dset)

  if (length(vec) != sum(mask)) {
    stop(sprintf(
      "Length of %s (%d) does not equal number of TRUE voxels in mask (%d)",
      map_dset, length(vec), sum(mask)
    ))
  }
  # Pass vector and space (extracted from mask) to constructor
  ClusteredNeuroVol(mask, vec)
}

#' Safely close an HDF5 object (file, group, dataset, attribute, dataspace, property list)
#'
#' Checks if the object is valid and not NULL before attempting to close.
#'
#' @keywords internal
close_h5_safely <- function(h5obj) {
  if (!is.null(h5obj) && h5obj$is_valid) {
    try(h5obj$close(), silent = TRUE)
  }
}
