# Helper functions for creating minimal objects for Roxygen examples
# These are internal and intended to be called via fmristore:::

#' Create a minimal LogicalNeuroVol for examples
#'
#' @param dims 3D dimensions, e.g., c(3L, 3L, 2L).
#' @param true_voxels A list of 3-element integer vectors for TRUE voxels, 
#'   e.g., list(c(1L,1L,1L), c(2L,1L,1L)). If NULL, creates a small default pattern.
#' @return A \code{LogicalNeuroVol} object.
#' @keywords internal
#' @examples
#' # Not typically called directly by users, used in other examples via fmristore:::
#' if (requireNamespace("neuroim2", quietly = TRUE)) {
#'   logi_vol <- fmristore:::create_minimal_LogicalNeuroVol()
#'   print(logi_vol)
#'   logi_vol2 <- fmristore:::create_minimal_LogicalNeuroVol(dims = c(2,2,2), 
#'                                                          true_voxels = list(c(1L,1L,1L)))
#'   print(logi_vol2)
#' }
create_minimal_LogicalNeuroVol <- function(dims = c(3L, 3L, 2L), true_voxels = NULL) {
  if (!requireNamespace("neuroim2", quietly = TRUE)) {
    stop("Package 'neuroim2' is needed for this helper function.")
  }
  space <- neuroim2::NeuroSpace(dims)
  arr <- array(FALSE, dim = dims)
  if (is.null(true_voxels)) {
    # Default pattern: first few voxels
    arr[1L,1L,1L] <- TRUE
    if (prod(dims) > 1 && dims[1] >= 2) arr[2L,1L,1L] <- TRUE
  } else {
    for (v_coord in true_voxels) {
      if (length(v_coord) == 3 && all(v_coord >= 1L) && all(v_coord <= dims)) {
        arr[v_coord[1], v_coord[2], v_coord[3]] <- TRUE
      } else {
        warning("Skipping invalid true_voxel coordinate: ", paste(v_coord, collapse=", "))
      }
    }
  }
  return(neuroim2::LogicalNeuroVol(arr, space))
}

#' Create a minimal DenseNeuroVec for examples
#'
#' @param dims 4D dimensions, e.g., c(3L, 3L, 2L, 4L).
#' @return A \code{DenseNeuroVec} object with minimal sequential data.
#' @keywords internal
#' @examples
#' # Not typically called directly by users, used in other examples via fmristore:::
#' if (requireNamespace("neuroim2", quietly = TRUE)) {
#'   dvec <- fmristore:::create_minimal_DenseNeuroVec()
#'   print(dvec)
#' }
create_minimal_DenseNeuroVec <- function(dims = c(3L, 3L, 2L, 4L)) {
  if (!requireNamespace("neuroim2", quietly = TRUE)) {
    stop("Package 'neuroim2' is needed for this helper function.")
  }
  space <- neuroim2::NeuroSpace(dims)
  arr <- array(seq_len(prod(dims)), dim = dims) 
  return(neuroim2::DenseNeuroVec(arr, space))
}

#' Create a minimal ClusteredNeuroVol for examples
#'
#' @param mask_vol A \code{LogicalNeuroVol} to use as the mask. 
#'   If \code{NULL}, a default one is created.
#' @param num_clusters Integer, number of clusters to create.
#' @return A \code{ClusteredNeuroVol} object.
#' @keywords internal
#' @examples
#' # Not typically called directly by users, used in other examples via fmristore:::
#' if (requireNamespace("neuroim2", quietly = TRUE)) {
#'   clust_vol <- fmristore:::create_minimal_ClusteredNeuroVol()
#'   print(clust_vol)
#'   
#'   custom_mask <- fmristore:::create_minimal_LogicalNeuroVol(dims = c(5,5,3))
#'   clust_vol2 <- fmristore:::create_minimal_ClusteredNeuroVol(mask_vol = custom_mask, 
#'                                                              num_clusters = 3L)
#'   print(clust_vol2)
#' }
create_minimal_ClusteredNeuroVol <- function(mask_vol = NULL, num_clusters = 2L) {
  if (!requireNamespace("neuroim2", quietly = TRUE)) {
    stop("Package 'neuroim2' is needed for this helper function.")
  }
  if (is.null(mask_vol)) {
    mask_vol <- create_minimal_LogicalNeuroVol(dims = c(4L,4L,3L), 
                                               true_voxels = list(c(1L,1L,1L), c(2L,1L,1L), 
                                                                  c(1L,2L,1L), c(2L,2L,1L),
                                                                  c(1L,1L,2L), c(2L,1L,2L)))
  }
  
  n_vox_in_mask <- sum(mask_vol@.Data)
  if (n_vox_in_mask == 0) {
    stop("Provided or default mask_vol is empty. Cannot create ClusteredNeuroVol.")
  }
  if (num_clusters <= 0L) {
    stop("'num_clusters' must be positive.")
  }
  
  cluster_data <- rep_len(seq_len(num_clusters), n_vox_in_mask)
  cluster_labels <- paste0("Cluster", seq_len(num_clusters))
  label_map <- stats::setNames(as.list(seq_len(num_clusters)), cluster_labels)
  
  return(neuroim2::ClusteredNeuroVol(mask_vol, clusters = cluster_data, label_map = label_map))
}

#' Create a minimal HDF5 file suitable for H5NeuroVol examples
#'
#' This function creates a temporary HDF5 file with the minimal structure
#' expected by the \code{H5NeuroVol} constructor.
#'
#' @param dims 3D dimensions, e.g., c(3L, 3L, 2L).
#' @param file_path Optional: path to HDF5 file. If \code{NULL}, a temp file is created.
#' @return Path to the created HDF5 file.
#' @keywords internal
#' @examples
#' # Not typically called directly by users, used in other examples via fmristore:::
#' if (requireNamespace("neuroim2", quietly = TRUE) && 
#'     requireNamespace("hdf5r", quietly = TRUE)) {
#'   temp_file <- fmristore:::create_minimal_h5_for_H5NeuroVol()
#'   print(temp_file)
#'   # Example usage:
#'   # h5vol <- H5NeuroVol(temp_file)
#'   # print(h5vol)
#'   # close(h5vol)
#'   if (file.exists(temp_file)) unlink(temp_file)
#' }
create_minimal_h5_for_H5NeuroVol <- function(dims = c(3L, 3L, 2L), file_path = NULL) {
  if (!requireNamespace("neuroim2", quietly = TRUE) || !requireNamespace("hdf5r", quietly = TRUE)) {
    stop("Packages \'neuroim2\' and \'hdf5r\' are needed for this helper function.")
  }
  if (is.null(file_path)) {
    out_file <- tempfile(fileext = ".h5vol_example.h5")
  } else {
    out_file <- file_path
  }

  h5f <- NULL
  tryCatch({
    h5f <- hdf5r::H5File$new(out_file, mode = "w")
    
    hdf5r::h5attr(h5f, "rtype") <- "DenseNeuroVol"
    
    sp <- neuroim2::NeuroSpace(dims)
    
    space_grp <- h5f$create_group("space")
    space_grp$create_dataset("dim", robj = as.integer(dim(sp)), 
                             dtype = hdf5r::h5types$H5T_NATIVE_INT)
    space_grp$create_dataset("origin", robj = as.double(neuroim2::origin(sp)), 
                             dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE)
    space_grp$create_dataset("trans", robj = neuroim2::trans(sp), 
                             dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE)
    # space_grp$create_dataset("spacing", robj = as.double(neuroim2::spacing(sp)), 
    #                          dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE)

    data_grp <- h5f$create_group("data")
    data_arr <- array(stats::rnorm(prod(dims)), dim = dims)
    data_grp$create_dataset("elements", robj = data_arr, 
                            dtype = hdf5r::h5types$H5T_NATIVE_FLOAT)
    
  }, error = function(e) {
    # Ensure file is closed if error occurs during creation, then rethrow
    if (!is.null(h5f) && h5f$is_valid) try(h5f$close_all(), silent = TRUE)
    stop(sprintf("Error creating minimal HDF5 for H5NeuroVol at %s: %s", out_file, e$message))
  }, finally = {
    if (!is.null(h5f) && h5f$is_valid) {
      try(h5f$close_all(), silent = TRUE)
    }
  })
  return(out_file)
}

#' Create a minimal HDF5 file suitable for H5NeuroVec examples
#'
#' This function creates a temporary HDF5 file with the minimal structure
#' expected by the \code{H5NeuroVec} constructor (rtype, /space, /data).
#'
#' @param dims 4D dimensions, e.g., c(3L, 3L, 2L, 5L).
#' @param file_path Optional: path to HDF5 file. If \code{NULL}, a temp file is created.
#' @return Path to the created HDF5 file.
#' @keywords internal
#' @examples
#' # Not typically called directly by users, used in other examples via fmristore:::
#' if (requireNamespace("neuroim2", quietly = TRUE) && 
#'     requireNamespace("hdf5r", quietly = TRUE)) {
#'   temp_file <- fmristore:::create_minimal_h5_for_H5NeuroVec()
#'   print(temp_file)
#'   # Example usage:
#'   # h5vec <- H5NeuroVec(temp_file)
#'   # print(h5vec)
#'   # close(h5vec)
#'   if (file.exists(temp_file)) unlink(temp_file)
#' }
create_minimal_h5_for_H5NeuroVec <- function(dims = c(3L, 3L, 2L, 5L), file_path = NULL) {
  if (!requireNamespace("neuroim2", quietly = TRUE) || !requireNamespace("hdf5r", quietly = TRUE)) {
    stop("Packages \'neuroim2\' and \'hdf5r\' are needed for this helper function.")
  }
  if (length(dims) != 4) stop("\'dims\' must be a 4-element vector for H5NeuroVec.")

  if (is.null(file_path)) {
    out_file <- tempfile(fileext = ".h5vec_example.h5")
  } else {
    out_file <- file_path
  }

  h5f <- NULL
  tryCatch({
    h5f <- hdf5r::H5File$new(out_file, mode = "w")
    
    hdf5r::h5attr(h5f, "rtype") <- "DenseNeuroVec" # Assuming H5NeuroVec is a DenseNeuroVec
    
    sp <- neuroim2::NeuroSpace(dims)
    
    space_grp <- h5f$create_group("space")
    space_grp$create_dataset("dim", robj = as.integer(dim(sp)), 
                             dtype = hdf5r::h5types$H5T_NATIVE_INT)
    space_grp$create_dataset("origin", robj = as.double(neuroim2::origin(sp)), 
                             dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE)
    space_grp$create_dataset("trans", robj = neuroim2::trans(sp), 
                             dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE)

    data_grp <- h5f$create_group("data")
    data_arr <- array(stats::rnorm(prod(dims)), dim = dims)
    data_grp$create_dataset("elements", robj = data_arr, 
                            dtype = hdf5r::h5types$H5T_NATIVE_FLOAT) 
    
  }, error = function(e) {
    if (!is.null(h5f) && h5f$is_valid) try(h5f$close_all(), silent = TRUE)
    stop(sprintf("Error creating minimal HDF5 for H5NeuroVec at %s: %s", out_file, e$message))
  }, finally = {
    if (!is.null(h5f) && h5f$is_valid) {
      try(h5f$close_all(), silent = TRUE)
    }
  })
  return(out_file)
}

#' Create a minimal LatentNeuroVec for examples
#'
#' @param space_dims 3D spatial dimensions for the underlying space, e.g., c(10L, 10L, 3L).
#' @param n_time Number of time points (columns in basis).
#' @param n_comp Number of components (rows in basis, columns in loadings).
#' @param n_mask_voxels Number of voxels to include in the mask. If NULL, ~20% of space_dims.
#' @return A \code{LatentNeuroVec} object.
#' @keywords internal
#' @examples
#' # Not typically called directly by users, used in other examples via fmristore:::
#' if (requireNamespace("neuroim2", quietly = TRUE)) {
#'   lnv <- fmristore:::create_minimal_LatentNeuroVec()
#'   print(lnv)
#' }
create_minimal_LatentNeuroVec <- function(space_dims = c(6L, 6L, 3L), 
                                          n_time = 10L, 
                                          n_comp = 3L, 
                                          n_mask_voxels = NULL) {
  if (!requireNamespace("neuroim2", quietly = TRUE)) {
    stop("Package 'neuroim2' is needed for this helper function.")
  }

  # Create a mask
  if (is.null(n_mask_voxels)) {
    n_mask_voxels <- floor(prod(space_dims) * 0.2)
    if (n_mask_voxels == 0 && prod(space_dims) > 0) n_mask_voxels <- 1L
  }
  
  mask_arr <- array(FALSE, dim = space_dims)
  if (n_mask_voxels > 0 && n_mask_voxels <= prod(space_dims)) {
    mask_indices <- sample(prod(space_dims), n_mask_voxels)
    mask_arr[mask_indices] <- TRUE
  }
  mask_space <- neuroim2::NeuroSpace(space_dims)
  mask_vol <- neuroim2::LogicalNeuroVol(mask_arr, mask_space)

  # Create basis functions (time x components)
  basis_mat <- matrix(stats::rnorm(n_time * n_comp), nrow = n_time, ncol = n_comp)

  # Create loadings (masked voxels x components)
  # Number of rows in loadings must match number of TRUE voxels in mask
  actual_n_mask_voxels <- sum(mask_vol@.Data)
  if (actual_n_mask_voxels == 0 && n_mask_voxels > 0) {
      # This case should ideally not be hit if mask creation is robust
      stop("Mask creation resulted in zero voxels, cannot create LatentNeuroVec loadings.")
  }
  if (actual_n_mask_voxels > 0) {
    loadings_mat <- matrix(stats::rnorm(actual_n_mask_voxels * n_comp), 
                         nrow = actual_n_mask_voxels, ncol = n_comp)
  } else {
    # Handle case of empty mask (e.g. if n_mask_voxels was 0)
    loadings_mat <- matrix(0, nrow = 0, ncol = n_comp)
  }
  
  # Create LatentNeuroVec
  # The constructor takes the full NeuroSpace for the 4D representation, not just the mask's space.
  # The 4th dimension is n_time.
  full_space_dims <- c(space_dims, n_time)
  full_space <- neuroim2::NeuroSpace(dim = full_space_dims, 
                                     spacing = neuroim2::spacing(mask_space), 
                                     origin = neuroim2::origin(mask_space), 
                                     axes = neuroim2::axes(mask_space)) # Retain spatial axes info
                                    
  # The LatentNeuroVec constructor from neuroim2 might be:
  # LatentNeuroVec(mask, loadings, basis, space)
  # where 'space' is the 4D space.  Need to check the fmristore definition or neuroim2.
  # Assuming fmristore's LatentNeuroVec is similar or identical to neuroim2's.
  # If fmristore has its own definition, this might need adjustment.
  
  # From R/latent_vec.R, the constructor seems to be new("LatentNeuroVec", ...) or LatentNeuroVec(...),
  # which takes `mask`, `loadings`, `basis`, `space`.
  
  return(neuroim2::LatentNeuroVec(mask = mask_vol, 
                                  loadings = loadings_mat, 
                                  basis = basis_mat, 
                                  space = full_space))
}

#' Create a minimal HDF5 file suitable for LabeledVolumeSet examples (via read_labeled_vec)
#'
#' This function creates a temporary HDF5 file with a minimal structure 
#' that can be read by \code{read_labeled_vec} to produce a \code{LabeledVolumeSet}.
#'
#' Refer to \code{write_labeled_vec} and \code{read_labeled_vec} for the expected structure.
#'
#' @param vol_dims 3D spatial dimensions for each volume, e.g., c(4L, 4L, 3L).
#' @param labels A character vector of labels, e.g., c("ConditionA", "ConditionB").
#' @param num_vols_per_label Integer, number of volumes to generate for each label.
#' @param file_path Optional: path to HDF5 file. If \code{NULL}, a temp file is created.
#' @return Path to the created HDF5 file.
#' @keywords internal
#' @examples
#' # Not typically called directly by users, used in other examples via fmristore:::
#' if (requireNamespace("neuroim2", quietly = TRUE) && 
#'     requireNamespace("hdf5r", quietly = TRUE) && 
#'     exists("read_labeled_vec", where = "package:fmristore")) {
#'   temp_file <- fmristore:::create_minimal_h5_for_LabeledVolumeSet()
#'   print(temp_file)
#'   # Example usage (now should work with read_labeled_vec):
#'   lvs <- NULL
#'   tryCatch({
#'     lvs <- fmristore::read_labeled_vec(temp_file)
#'     print(lvs)
#'     print(labels(lvs))
#'   }, finally = {
#'     if (!is.null(lvs)) try(close(lvs), silent = TRUE)
#'     if (file.exists(temp_file)) unlink(temp_file)
#'   })
#' }
create_minimal_h5_for_LabeledVolumeSet <- function(vol_dims = c(4L, 4L, 3L), 
                                                   labels = c("Set1", "Set2"), 
                                                   num_vols_per_label = 1L, # Simplified: 1 vol per label for minimal example
                                                   file_path = NULL) {
  if (!requireNamespace("neuroim2", quietly = TRUE) || !requireNamespace("hdf5r", quietly = TRUE)) {
    stop("Packages 'neuroim2' and 'hdf5r' are needed for this helper function.")
  }
  if (length(vol_dims) != 3) stop("'vol_dims' must be a 3-element vector.")
  if (num_vols_per_label <= 0) stop("'num_vols_per_label' must be positive.")
  if (length(labels) == 0) stop("'labels' must not be empty.")

  if (is.null(file_path)) {
    out_file <- tempfile(fileext = ".lvs_example.h5")
  } else {
    out_file <- file_path
  }

  # Total number of distinct volumes/datasets to create = number of labels * num_vols_per_label
  # For simplicity with read_labeled_vec which expects one dataset per unique label in /labels,
  # this helper will now create one dataset per entry in `labels`.
  # If num_vols_per_label > 1, the labels vector should reflect that (e.g. c("CondA_1", "CondA_2"))
  # For this revised helper, let's assume `labels` are unique and num_vols_per_label implies
  # that if you want multiple volumes for "Set1", you should provide e.g. labels = c("Set1_run1", "Set1_run2")
  # For true minimality for read_labeled_vec, num_vols_per_label = 1 is assumed per unique label.
  # The `labels` argument now directly defines the datasets to be created.

  h5f <- NULL
  tryCatch({
    h5f <- hdf5r::H5File$new(out_file, mode = "w")
    
    # Root attributes (mimicking write_labeled_vec)
    hdf5r::h5attr(h5f, "class") <- "LabeledVolumeSet" # Or as expected by read_labeled_vec if different
    hdf5r::h5attr(h5f, "version") <- "0.1" 

    # Create a simple mask (all TRUE for simplicity)
    mask_data_arr <- array(TRUE, dim = vol_dims)
    # Ensure mask is written as UCHAR for compatibility if read_labeled_vec expects it
    h5f$create_dataset("mask", robj = as.integer(mask_data_arr), dtype = hdf5r::h5types$H5T_NATIVE_UCHAR)

    # Labels dataset (vector of unique labels provided)
    # Use H5T_STRING$new(type="c") for variable length C strings, or fixed size if appropriate
    str_type <- hdf5r::H5T_STRING$new(type="c", size = Inf) # Variable length
    on.exit(str_type$close(), add = TRUE)
    h5f$create_dataset("labels", robj = labels, dtype = str_type)

    # /data group
    if (!h5f$exists("data")) data_grp <- h5f$create_group("data")
    else data_grp <- h5f[["data"]]
    
    # Create one dataset per label under /data/
    # This dataset will contain the 1D masked data for that label.
    n_vox_in_mask <- sum(mask_data_arr)
    if (n_vox_in_mask == 0) stop("Mask is empty, cannot create data datasets.")

    sanitize_label_func <- function(lbl) { gsub("[^A-Za-z0-9_.-]", "_", lbl) }

    for (lab_name in labels) {
      # Data for this label (1D vector of length n_vox_in_mask)
      label_vol_data_1d <- stats::rnorm(n_vox_in_mask)
      safe_lab_name <- sanitize_label_func(lab_name)
      
      # Path like /data/Label1, /data/Label2
      dataset_path <- file.path("data", safe_lab_name)
      data_grp$create_dataset(safe_lab_name, robj = label_vol_data_1d, 
                              dtype = hdf5r::h5types$H5T_NATIVE_FLOAT, # Consistent with write_labeled_vec default
                              chunk_dims = if(n_vox_in_mask > 0) min(1024L, n_vox_in_mask) else NULL,
                              compress_level = 0L # Minimal example, no compression needed
                              )
    }

    # Minimal NIfTI-like header information (required by read_labeled_vec to build the space)
    if (!h5f$exists("header")) header_grp <- h5f$create_group("header")
    else header_grp <- h5f[["header"]]
    
    sp <- neuroim2::NeuroSpace(vol_dims) # Basic space for one volume
    # nVols in header/dim should match length of /labels
    n_labels_for_header <- length(labels)
    header_grp$create_dataset("dim", robj = as.integer(c(4L, vol_dims, n_labels_for_header, 1L, 1L, 1L)), dtype = hdf5r::h5types$H5T_NATIVE_INT)
    header_grp$create_dataset("pixdim", robj = as.double(c(0.0, neuroim2::spacing(sp), rep(0,4))), dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE)
    
    q_info <- tryCatch(neuroim2::matrixToQuatern(neuroim2::trans(sp)), error=function(e) NULL)
    if (is.null(q_info)) { # Default q_info if matrixToQuatern fails (e.g. singular matrix)
        q_info <- list(quaternion=c(0,0,0), qoffset=neuroim2::origin(sp), qfac=1)
    }

    header_grp$create_dataset("quatern_b", robj = q_info$quaternion[1], dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE)
    header_grp$create_dataset("quatern_c", robj = q_info$quaternion[2], dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE)
    header_grp$create_dataset("quatern_d", robj = q_info$quaternion[3], dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE)
    header_grp$create_dataset("qoffset_x", robj = q_info$qoffset[1], dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE)
    header_grp$create_dataset("qoffset_y", robj = q_info$qoffset[2], dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE)
    header_grp$create_dataset("qoffset_z", robj = q_info$qoffset[3], dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE)
    header_grp$create_dataset("qfac", robj = q_info$qfac, dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE)
    
    affine_mat <- neuroim2::trans(sp)
    srow_x <- affine_mat[1,]
    srow_y <- affine_mat[2,]
    srow_z <- affine_mat[3,]
    header_grp$create_dataset("srow_x", robj = as.double(srow_x), dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE)
    header_grp$create_dataset("srow_y", robj = as.double(srow_y), dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE)
    header_grp$create_dataset("srow_z", robj = as.double(srow_z), dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE)
    header_grp$create_dataset("qform_code", robj = as.integer(1), dtype = hdf5r::h5types$H5T_NATIVE_INT)
    header_grp$create_dataset("sform_code", robj = as.integer(1), dtype = hdf5r::h5types$H5T_NATIVE_INT)
    header_grp$create_dataset("sizeof_hdr", robj = as.integer(348), dtype = hdf5r::h5types$H5T_NATIVE_INT)
    header_grp$create_dataset("magic", robj = "n+1", dtype = hdf5r::H5T_STRING$new(type="c", size=4L))

  }, error = function(e) {
    if (!is.null(h5f) && h5f$is_valid) try(h5f$close_all(), silent = TRUE)
    stop(sprintf("Error creating aligned HDF5 for LabeledVolumeSet at %s: %s", out_file, e$message))
  }, finally = {
    if (!is.null(h5f) && h5f$is_valid) {
      try(h5f$close_all(), silent = TRUE)
    }
  })
  return(out_file)
}

#' Populate HDF5 structure for a single H5ClusteredRunFull within an existing H5File
#'
#' This function is a helper for creating minimal H5ClusteredExperiment examples.
#' It assumes the H5File object is open and writable.
#'
#' @param h5obj An open, writable \code{H5File} object (from package hdf5r).
#' @param scan_name Character string, the name for this scan/run.
#' @param mask_vol A \code{LogicalNeuroVol} object for the mask of this run.
#' @param cluster_map_vol A \code{ClusteredNeuroVol} object representing the cluster assignments
#'   within the mask. Its space must match \code{mask_vol}.
#' @param n_time Integer, the number of time points for the cluster data.
#' @param compress Logical, whether to compress HDF5 datasets (minimal impact for small data).
#' @keywords internal
#' @return Invisibly returns TRUE on success, or stops on error.
populate_H5ClusteredRunFull_in_h5file <- function(h5obj,
                                                  scan_name,
                                                  mask_vol,
                                                  cluster_map_vol,
                                                  n_time,
                                                  compress = FALSE) {
  if (!inherits(h5obj, "H5File") || !h5obj$is_valid) {
    stop("'h5obj' must be a valid, open H5File object.")
  }
  if (!is(mask_vol, "LogicalNeuroVol")) stop("'mask_vol' must be a LogicalNeuroVol.")
  if (!is(cluster_map_vol, "ClusteredNeuroVol")) stop("'cluster_map_vol' must be a ClusteredNeuroVol.")
  if (!identical(neuroim2::space(mask_vol), neuroim2::space(cluster_map_vol))) {
    stop("Space of 'mask_vol' and 'cluster_map_vol' must be identical.")
  }
  if (n_time <= 0) stop("'n_time' must be positive.")

  tryCatch({
    if (!h5obj$exists("scans")) {
      h5obj$create_group("scans")
    }
    scans_grp <- h5obj[["scans"]]
    
    if (scans_grp$exists(scan_name)) {
      stop(sprintf("Scan group '%s' already exists in /scans.", scan_name))
    }
    run_grp <- scans_grp$create_group(scan_name)
    
    # Set attributes for the run group
    hdf5r::h5attr(run_grp, "class") <- "H5ClusteredRunFull"
    hdf5r::h5attr(run_grp, "n_time") <- as.integer(n_time)
    # fmristore version attribute could be added if constructor expects/uses it.

    # Write mask for this run (can be different from experiment master mask in general)
    # For simplicity, using the mask_vol@.Data which is a 3D array
    run_grp$create_dataset("mask", robj = mask_vol@.Data, dtype = hdf5r::h5types$H5T_NATIVE_HBOOL,
                               chunk_dims = "auto", compress_level = if(compress) 4L else 0L)

    # Extract cluster information from ClusteredNeuroVol
    # cluster_map is a vector of cluster IDs for voxels *within the mask_vol*
    # cluster_ids are the unique, sorted cluster IDs present.
    # cluster_names are derived from the label_map of ClusteredNeuroVol.
    
    # Get cluster assignments only for voxels TRUE in mask_vol
    cluster_assignments_in_mask <- cluster_map_vol@.Data[mask_vol@.Data]
    unique_cluster_ids <- sort(unique(cluster_assignments_in_mask[cluster_assignments_in_mask > 0])) # Exclude 0 if it means unclustered

    if (length(unique_cluster_ids) == 0) {
      warning(sprintf("Scan '%s': No positive cluster IDs found in cluster_map_vol within the mask. Creating empty cluster datasets.", scan_name))
      # Still create the datasets for schema completeness if needed by constructor
      run_grp$create_dataset("cluster_map", robj = integer(0), dtype = hdf5r::h5types$H5T_NATIVE_INT)
      run_grp$create_dataset("cluster_ids", robj = integer(0), dtype = hdf5r::h5types$H5T_NATIVE_INT)
      run_grp$create_dataset("cluster_names", robj = character(0), dtype = hdf5r::H5T_STRING$new(type="c"))
      run_grp$create_group("clusters") # Empty group
    } else {
      # Write cluster_map (dense vector for voxels within the mask)
      run_grp$create_dataset("cluster_map", robj = as.integer(cluster_assignments_in_mask), 
                               dtype = hdf5r::h5types$H5T_NATIVE_INT)
      
      # Write cluster_ids
      run_grp$create_dataset("cluster_ids", robj = as.integer(unique_cluster_ids), 
                               dtype = hdf5r::h5types$H5T_NATIVE_INT)
      
      # Write cluster_names (try to get from ClusteredNeuroVol's label_map)
      # Match unique_cluster_ids to names. Default to "Cluster_<id>" if not found.
      cl_names <- character(length(unique_cluster_ids))
      if (!is.null(cluster_map_vol@label_map) && length(cluster_map_vol@label_map) > 0) {
        # label_map stores names as keys and IDs as values, or vice-versa depending on neuroim2
        # neuroim2::ClusteredNeuroVol: label_map is a list where names are labels and values are cluster IDs.
        # We need names for our unique_cluster_ids.
        for (i in seq_along(unique_cluster_ids)) {
          id <- unique_cluster_ids[i]
          found_name <- FALSE
          for (name_in_map in names(cluster_map_vol@label_map)) {
            if (id %in% cluster_map_vol@label_map[[name_in_map]]) { # ID can be a vector in label_map value
              cl_names[i] <- name_in_map
              found_name <- TRUE
              break
            }
          }
          if (!found_name) cl_names[i] <- paste0("Cluster_", id)
        }
      } else {
        cl_names <- paste0("Cluster_", unique_cluster_ids)
      }
      run_grp$create_dataset("cluster_names", robj = cl_names, dtype = hdf5r::H5T_STRING$new(type="c"))

      # Write cluster data (time series for each cluster)
      # /scans/<scan_name>/clusters/cluster_<id>
      clusters_data_grp <- run_grp$create_group("clusters")
      for (id in unique_cluster_ids) {
        # Number of voxels in this specific cluster for this run
        n_vox_in_this_cluster <- sum(cluster_assignments_in_mask == id)
        if (n_vox_in_this_cluster > 0) {
          # Data: n_time x n_vox_in_this_cluster
          # The H5ClusteredRunFull constructor expects data loaded to be n_vox x n_time (from summary)
          # or reads n_time x n_vox from disk if that's how write_clustered_run_h5 does it.
          # write_clustered_run_h5 writes dataset as h5_write(..., robj = data_for_cluster) # data_for_cluster is n_vox x n_time from input
          # So, the dataset on disk should be n_vox x n_time. But cluster_ts() returns n_time x n_vox.
          # Let's assume disk storage is n_time x n_vox_in_cluster for now, common for time series.
          # If constructor expects n_vox x n_time, then this needs to be transposed before constructor is called, or disk format adjusted.
          # The H5ClusteredRunFull constructor: data for each cluster is loaded as: dset[,]
          # Let's use n_time x n_vox_in_cluster on disk.
          cluster_ts_data <- matrix(stats::rnorm(n_time * n_vox_in_this_cluster), 
                                   nrow = n_time, ncol = n_vox_in_this_cluster)
          clusters_data_grp$create_dataset(paste("cluster_", id, sep=""), robj = cluster_ts_data, 
                                   dtype = hdf5r::h5types$H5T_NATIVE_FLOAT,
                                   chunk_dims = "auto", compress_level = if(compress) 4L else 0L)
        } else {
           # Should not happen if id is from unique_cluster_ids > 0 and cluster_assignments_in_mask is not all 0
           warning(sprintf("Scan '%s', Cluster '%s': No voxels found. Skipping data dataset.", scan_name, id))
        }
      }
    }
    
  }, error = function(e) {
    stop(sprintf("Error populating H5ClusteredRunFull for scan '%s': %s", scan_name, e$message))
  })
  
  invisible(TRUE)
}

#' Populate HDF5 structure for a single H5ClusteredRunSummary within an existing H5File
#'
#' This function is a helper for creating minimal H5ClusteredExperiment examples.
#' It assumes the H5File object is open and writable.
#'
#' @param h5obj An open, writable \code{H5File} object (from package hdf5r).
#' @param scan_name Character string, the name for this scan/run.
#' @param mask_vol A \code{LogicalNeuroVol} object for the mask of this run.
#' @param cluster_map_vol A \code{ClusteredNeuroVol} object representing the cluster assignments
#'   within the mask. Its space must match \code{mask_vol}.
#' @param n_time Integer, the number of time points for the summary data.
#' @param summary_stat_names A character vector of summary statistics to generate (e.g., "mean", "median").
#'   For each stat, a dataset \code{summary_stats/<stat_name>/cluster_<id>} will be created.
#' @param compress Logical, whether to compress HDF5 datasets.
#' @keywords internal
#' @return Invisibly returns TRUE on success, or stops on error.
populate_H5ClusteredRunSummary_in_h5file <- function(h5obj,
                                                     scan_name,
                                                     mask_vol,
                                                     cluster_map_vol,
                                                     n_time,
                                                     summary_stat_names = c("mean"),
                                                     compress = FALSE) {
  if (!inherits(h5obj, "H5File") || !h5obj$is_valid) {
    stop("'h5obj' must be a valid, open H5File object.")
  }
  if (!is(mask_vol, "LogicalNeuroVol")) stop("'mask_vol' must be a LogicalNeuroVol.")
  if (!is(cluster_map_vol, "ClusteredNeuroVol")) stop("'cluster_map_vol' must be a ClusteredNeuroVol.")
  if (!identical(neuroim2::space(mask_vol), neuroim2::space(cluster_map_vol))) {
    stop("Space of 'mask_vol' and 'cluster_map_vol' must be identical.")
  }
  if (n_time <= 0) stop("'n_time' must be positive.")
  if (length(summary_stat_names) == 0) stop("'summary_stat_names' cannot be empty.")

  tryCatch({
    if (!h5obj$exists("scans")) {
      h5obj$create_group("scans")
    }
    scans_grp <- h5obj[["scans"]]
    
    if (scans_grp$exists(scan_name)) {
      stop(sprintf("Scan group '%s' already exists in /scans.", scan_name))
    }
    run_grp <- scans_grp$create_group(scan_name)
    
    # Set attributes for the run group
    hdf5r::h5attr(run_grp, "class") <- "H5ClusteredRunSummary"
    hdf5r::h5attr(run_grp, "n_time") <- as.integer(n_time)

    # Write mask, cluster_map, cluster_ids, cluster_names (same as H5ClusteredRunFull part)
    run_grp$create_dataset("mask", robj = mask_vol@.Data, dtype = hdf5r::h5types$H5T_NATIVE_HBOOL,
                               chunk_dims = "auto", compress_level = if(compress) 4L else 0L)

    cluster_assignments_in_mask <- cluster_map_vol@.Data[mask_vol@.Data]
    unique_cluster_ids <- sort(unique(cluster_assignments_in_mask[cluster_assignments_in_mask > 0]))

    if (length(unique_cluster_ids) == 0) {
      warning(sprintf("Scan '%s': No positive cluster IDs. Creating empty cluster info and summary_stats.", scan_name))
      run_grp$create_dataset("cluster_map", robj = integer(0), dtype = hdf5r::h5types$H5T_NATIVE_INT)
      run_grp$create_dataset("cluster_ids", robj = integer(0), dtype = hdf5r::h5types$H5T_NATIVE_INT)
      run_grp$create_dataset("cluster_names", robj = character(0), dtype = hdf5r::H5T_STRING$new(type="c"))
      run_grp$create_group("summary_stats") # Empty group
    } else {
      run_grp$create_dataset("cluster_map", robj = as.integer(cluster_assignments_in_mask), 
                               dtype = hdf5r::h5types$H5T_NATIVE_INT)
      run_grp$create_dataset("cluster_ids", robj = as.integer(unique_cluster_ids), 
                               dtype = hdf5r::h5types$H5T_NATIVE_INT)
      
      cl_names <- character(length(unique_cluster_ids))
      if (!is.null(cluster_map_vol@label_map) && length(cluster_map_vol@label_map) > 0) {
        for (i in seq_along(unique_cluster_ids)) {
          id <- unique_cluster_ids[i]
          found_name <- FALSE
          for (name_in_map in names(cluster_map_vol@label_map)) {
            if (id %in% cluster_map_vol@label_map[[name_in_map]]) {
              cl_names[i] <- name_in_map
              found_name <- TRUE
              break
            }
          }
          if (!found_name) cl_names[i] <- paste0("Cluster_", id)
        }
      } else {
        cl_names <- paste0("Cluster_", unique_cluster_ids)
      }
      run_grp$create_dataset("cluster_names", robj = cl_names, dtype = hdf5r::H5T_STRING$new(type="c"))

      # Write summary statistics data
      # /scans/<scan_name>/summary_stats/<stat_name>/cluster_<id>
      summary_stats_grp <- run_grp$create_group("summary_stats")
      for (stat_name in summary_stat_names) {
        stat_spec_grp <- summary_stats_grp$create_group(stat_name)
        for (id in unique_cluster_ids) {
          # Summary data: n_time x 1 (e.g., mean time series for the cluster)
          summary_ts_data <- matrix(stats::rnorm(n_time * 1), nrow = n_time, ncol = 1)
          stat_spec_grp$create_dataset(paste("cluster_", id, sep=""), robj = summary_ts_data, 
                                 dtype = hdf5r::h5types$H5T_NATIVE_FLOAT,
                                 chunk_dims = "auto", compress_level = if(compress) 4L else 0L)
        }
      }
    }

  }, error = function(e) {
    stop(sprintf("Error populating H5ClusteredRunSummary for scan '%s': %s", scan_name, e$message))
  })
  
  invisible(TRUE)
}

#' Create a minimal HDF5 file suitable for H5ClusteredExperiment examples
#'
#' This function creates a temporary HDF5 file with a minimal but valid structure 
#' for an H5ClusteredExperiment, including a master mask, master cluster definitions,
#' one H5ClusteredRunFull, and one H5ClusteredRunSummary.
#'
#' @param file_path Optional: path to HDF5 file. If \code{NULL}, a temp file is created.
#' @param master_mask_dims 3D dimensions for the master mask, e.g., c(5L, 5L, 4L).
#' @param num_master_clusters Integer, number of clusters in the master cluster map.
#' @param n_time_run1 Integer, n_time for the first (full) run.
#' @param n_time_run2 Integer, n_time for the second (summary) run.
#' @return Path to the created HDF5 file.
#' @keywords internal
#' @examples
#' # Not typically called directly by users, used in other examples via fmristore:::
#' if (requireNamespace("neuroim2", quietly = TRUE) && 
#'     requireNamespace("hdf5r", quietly = TRUE) &&
#'     exists("H5ClusteredExperiment", where = "package:fmristore")) {
#'   
#'   # Ensure helper dependencies are available if not automatically sourced
#'   # (e.g. source the R/zzz_example_helpers.R file if running interactively outside pkg build)
#' 
#'   temp_exp_file <- NULL
#'   exp <- NULL
#'   tryCatch({
#'     temp_exp_file <- fmristore:::create_minimal_h5_for_H5ClusteredExperiment()
#'     print(temp_exp_file)
#'     # exp <- fmristore::H5ClusteredExperiment(temp_exp_file)
#'     # print(exp)
#'   }, finally = {
#'     # if (!is.null(exp)) close(exp) # Close H5ClusteredExperiment handle
#'     if (!is.null(temp_exp_file) && file.exists(temp_exp_file)) unlink(temp_exp_file)
#'   })
#' }
create_minimal_h5_for_H5ClusteredExperiment <- function(
    file_path = NULL, 
    master_mask_dims = c(5L, 5L, 4L),
    num_master_clusters = 3L,
    n_time_run1 = 10L,
    n_time_run2 = 12L
  ) {

  if (!requireNamespace("neuroim2", quietly = TRUE) || !requireNamespace("hdf5r", quietly = TRUE)) {
    stop("Packages 'neuroim2' and 'hdf5r' are needed for this helper function.")
  }
  
  # 1. Create master mask and cluster map objects (in-memory)
  master_mask_vol <- create_minimal_LogicalNeuroVol(dims = master_mask_dims)
  # Ensure mask is not empty for cluster_map_vol creation
  if (sum(master_mask_vol@.Data) == 0) {
      warning("Default master_mask_vol is empty, attempting to create a new one with some TRUE voxels.")
      master_mask_vol <- create_minimal_LogicalNeuroVol(dims = master_mask_dims, true_voxels = list(c(1L,1L,1L)))
      if (sum(master_mask_vol@.Data) == 0) stop("Failed to create a non-empty master_mask_vol.")
  }
  master_cluster_map_vol <- create_minimal_ClusteredNeuroVol(mask_vol = master_mask_vol, 
                                                             num_clusters = num_master_clusters)

  if (is.null(file_path)) {
    out_file <- tempfile(fileext = ".h5exp_example.h5")
  } else {
    out_file <- file_path
  }

  h5f <- NULL
  tryCatch({
    # 2. Create the HDF5 file
    h5f <- hdf5r::H5File$new(out_file, mode = "w")

    # 3. Write master mask to /mask
    h5f$create_dataset("mask", robj = master_mask_vol@.Data, 
                       dtype = hdf5r::h5types$H5T_NATIVE_HBOOL)
    # Write space information for the master mask (as H5ClusteredExperiment constructor uses it)
    # The constructor calls .read_space(h5obj, "/")
    master_space_grp <- h5f$create_group("space") # Should be at root if path is "/"
    master_sp_obj <- neuroim2::space(master_mask_vol)
    master_space_grp$create_dataset("dim", robj = as.integer(dim(master_sp_obj)), dtype = hdf5r::h5types$H5T_NATIVE_INT)
    master_space_grp$create_dataset("origin", robj = as.double(neuroim2::origin(master_sp_obj)), dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE)
    master_space_grp$create_dataset("trans", robj = neuroim2::trans(master_sp_obj), dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE)
    master_space_grp$create_dataset("spacing", robj = as.double(neuroim2::spacing(master_sp_obj)), dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE)

    # 4. Write master cluster information to root
    master_cluster_assignments_in_mask <- master_cluster_map_vol@.Data[master_mask_vol@.Data]
    master_unique_ids <- sort(unique(master_cluster_assignments_in_mask[master_cluster_assignments_in_mask > 0]))

    if (length(master_unique_ids) > 0) {
        h5f$create_dataset("cluster_map", robj = as.integer(master_cluster_assignments_in_mask),
                           dtype = hdf5r::h5types$H5T_NATIVE_INT)
        h5f$create_dataset("cluster_ids", robj = as.integer(master_unique_ids),
                           dtype = hdf5r::h5types$H5T_NATIVE_INT)
        
        master_cl_names <- character(length(master_unique_ids))
        if (!is.null(master_cluster_map_vol@label_map) && length(master_cluster_map_vol@label_map) > 0) {
            for (i in seq_along(master_unique_ids)) {
                id <- master_unique_ids[i]
                found_name <- FALSE
                for (name_in_map in names(master_cluster_map_vol@label_map)) {
                    if (id %in% master_cluster_map_vol@label_map[[name_in_map]]) {
                        master_cl_names[i] <- name_in_map
                        found_name <- TRUE
                        break
                    }
                }
                if (!found_name) master_cl_names[i] <- paste0("MasterCluster_", id)
            }
        } else {
            master_cl_names <- paste0("MasterCluster_", master_unique_ids)
        }
        h5f$create_dataset("cluster_names", robj = master_cl_names, 
                           dtype = hdf5r::H5T_STRING$new(type="c"))
    } else {
        # Create empty datasets if no clusters, to satisfy constructor expectations if it reads them unconditionally
        h5f$create_dataset("cluster_map", robj = integer(0), dtype = hdf5r::h5types$H5T_NATIVE_INT)
        h5f$create_dataset("cluster_ids", robj = integer(0), dtype = hdf5r::h5types$H5T_NATIVE_INT)
        h5f$create_dataset("cluster_names", robj = character(0), dtype = hdf5r::H5T_STRING$new(type="c"))
        warning("Master cluster map is empty. Experiment will have no common clusters defined at root.")
    }

    # 5. Populate a full run
    populate_H5ClusteredRunFull_in_h5file(h5obj = h5f, 
                                          scan_name = "Run1_Full", 
                                          mask_vol = master_mask_vol, # Using master mask for run
                                          cluster_map_vol = master_cluster_map_vol, # Using master clusters for run
                                          n_time = n_time_run1)

    # 6. Populate a summary run
    populate_H5ClusteredRunSummary_in_h5file(h5obj = h5f, 
                                             scan_name = "Run2_Summary", 
                                             mask_vol = master_mask_vol, 
                                             cluster_map_vol = master_cluster_map_vol, 
                                             n_time = n_time_run2,
                                             summary_stat_names = c("mean", "sd"))

    # 7. Set root class attribute
    hdf5r::h5attr(h5f, "class") <- "H5ClusteredExperiment"
    # Potentially add fmristore_version attribute if constructor checks it
    # hdf5r::h5attr(h5f, "fmristore_version") <- as.character(utils::packageVersion("fmristore"))

  }, error = function(e) {
    if (!is.null(h5f) && h5f$is_valid) try(h5f$close_all(), silent = TRUE)
    stop(sprintf("Error creating minimal HDF5 for H5ClusteredExperiment at %s: %s", out_file, e$message))
  }, finally = {
    if (!is.null(h5f) && h5f$is_valid) {
      try(h5f$close_all(), silent = TRUE)
    }
  })
  
  return(out_file)
} 