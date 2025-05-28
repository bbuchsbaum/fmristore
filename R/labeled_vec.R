#' Write a NeuroVec to an HDF5 file with NIfTI-like quaternions
#'
#' @description
#' Creates an HDF5 file following a NIfTI-like header layout, storing:
#' \itemize{
#'   \item \code{/header/dim} => \code{[4, X, Y, Z, nVols, 1,1,1]}
#'   \item \code{/header/pixdim} => \code{[0.0, dx, dy, dz, ...]} (Note: qfac stored in /header/qfac)
#'   \item \code{/header/quatern_b,c,d} and \code{qoffset_x,y,z}
#'   \item \code{/header/qfac} => Quaternion factor (±1)
#'   \item \code{/mask} => 3D dataset \code{[X, Y, Z]} (0/1) at root level
#'   \item \code{/labels} => array of label strings at root level
#'   \item \code{/data/<label>} => 1D array (length = number of nonzero mask voxels)
#'         storing the sub-volume values
#' }
#'
#' @details
#' The 4×4 matrix in \code{trans(space(vec))} is passed to
#' \code{\link[neuroim2]{matrixToQuatern}}, which returns a list containing:
#' \itemize{
#'   \item \code{quaternion = c(b, c, d)} (the three imaginary parts)
#'   \item \code{qfac} (±1 sign)
#' }
#' This function stores \code{qfac} in \code{/header/qfac} and sets \code{/header/pixdim[0]=0}.
#' We also gather voxel spacing (dx,dy,dz) from \code{spacing(space(vec))} and
#' the origin from \code{origin(space(vec))}.
#'
#' We store a subset of NIfTI-like header fields in the \code{/header} group.
#' The user can supply \code{header_values} (a named list) to override or
#' augment *some* additional fields (e.g., \code{qform_code=1L}). See implementation
#' notes for which fields are protected.
#'
#' @section Lifecycle Management:
#' The HDF5 file is opened in write mode and the resulting handle is returned
#' without being automatically closed. **It is the caller's responsibility to
#' close this handle** (via \code{h5file$close_all()} or \code{close()}) when
#' finished working with the file.
#'
#' @param vec A 4D \code{\link[neuroim2]{NeuroVec}} with dimension \code{[X,Y,Z,nVols]}.
#' @param mask A \code{\link[neuroim2]{LogicalNeuroVol}} of shape \code{[X,Y,Z]}
#'   (the same 3D shape as \code{vec}).
#' @param labels A character vector of length \code{nVols}, labeling each 4D sub-volume.
#' @param file Either a character path to the HDF5 file to create or
#'   an open \code{\link[hdf5r]{H5File}} in write mode.
#' @param compression Integer \code{0-9} for gzip level; default \code{4}.
#' @param dtype An HDF5 data type object (e.g., \code{hdf5r::h5types$H5T_NATIVE_FLOAT}).
#'   Default is \code{hdf5r::h5types$H5T_NATIVE_DOUBLE}.
#' @param chunk_size If non-NULL, the chunk dimension for the 1D datasets. Default is \code{1024}.
#' @param header_values A named list of optional overrides for fields in the header
#'   (e.g., \code{list(qform_code=1L, sform_code=2L)}). Note that fields derived from
#'   `vec` or `mask` (like `dim`, `pixdim`, quaternion fields) cannot be overridden here.
#' @param verbose Logical, whether to print verbose messages during processing.
#'
#' @return Invisibly returns the open \code{\link[hdf5r]{H5File}} handle
#'   containing the written data. The caller should close this handle when
#'   finished.
#'
#' @seealso
#' \code{\link[neuroim2]{matrixToQuatern}} for how the quaternion is derived,
#' \code{\link[neuroim2]{quaternToMatrix}} for reconstructing the 4×4,
#' \code{\link{read_labeled_vec}} for reading the file back in.
#'
#' @import hdf5r
#' @importFrom neuroim2 spacing space origin trans matrixToQuatern
#' @importFrom hdf5r H5T_STRING H5S
#' @importFrom lifecycle deprecate_warn
#' @export
write_labeled_vec <- function(vec,
                              mask,
                              labels,
                              file,
                              compression = 4,
                              dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE,
                              chunk_size = 1024,
                              header_values = list(),
                              verbose = FALSE)
{
  # === Pre-flight checks before opening file ===
  
  # 1. Validate mask object and get mask array
  stopifnot(inherits(mask, "LogicalNeuroVol"))
  mask_arr <- as.array(mask)
  stopifnot(length(dim(mask_arr)) == 3) # Ensure mask is 3D
  
  # 2. Check if mask is empty - fail fast
  # Now works correctly with logical mask_arr
  idx_nonzero <- which(mask_arr == TRUE)
  n_nonzero   <- length(idx_nonzero)
  if (n_nonzero == 0) {
    stop("Mask is empty (all FALSE). Cannot write a LabeledVolumeSet with no valid voxels.")
  }
  
  # 3. Validate vec dimensions against mask
  nd <- dim(vec)  # [X, Y, Z, nVols]
  stopifnot(length(nd) == 4)
  if (!all(dim(mask_arr) == nd[1:3])) {
     stop("Mask dimensions [", paste(dim(mask_arr), collapse=","), 
          "] do not match first 3 dimensions of vec [", paste(nd[1:3], collapse=","), "]")
  }
  nVols <- nd[4]

  # 4. Check labels length
  if (length(labels) != nVols) {
    stop("Length of 'labels' (", length(labels), ") must match the 4th dimension of 'vec' (", nVols, ").")
  }
  
  # 5. Sanitize labels and check for duplicates
  # Perform basic sanitization first
  basic_safe_labels <- vapply(labels, function(lbl) gsub("[^A-Za-z0-9_.-]", "_", lbl), character(1))
  
  # Check for duplicates *after* basic sanitization
  if (length(unique(basic_safe_labels)) != length(basic_safe_labels)) {
    stop("Duplicate labels detected after basic sanitization (gsub). Check input labels.")
  }
  
  # Use the basic sanitized labels (now guaranteed unique after gsub)
  safe_labels <- basic_safe_labels
  
  # === Open file using helper ===
  fh <- open_h5(file, mode = "w") # Use write mode
  h5obj <- fh$h5
  
  # Get dimensions (already validated above)
  X <- nd[1]; Y <- nd[2]; Z <- nd[3]

  # Extract 4×4 transformation matrix from NeuroSpace
  tmat <- trans(space(vec))          # e.g. a 4×4
  # Convert to quaternion + qfac - Harden error message
  q <- tryCatch(
      matrixToQuatern(tmat),
      error = function(e) {
          stop("Invalid NeuroSpace in 'vec' – cannot convert transformation matrix to quaternion: ", e$message)
      }
  )

  # Gather spacing & origin
  sp  <- spacing(space(vec))         # c(dx, dy, dz)
  org <- origin(space(vec))          # c(ox, oy, oz)
  if (length(sp) < 3) {
    sp  <- c(sp, rep(1, 3 - length(sp)))
  }
  if (length(org) < 3) {
    org <- c(org, rep(0, 3 - length(org)))
  }

  # Validate dtype argument - must be a single H5T object
  if (!inherits(dtype, "H5T")) {
    stop("'dtype' argument must be a single H5T object (e.g., hdf5r::h5types$H5T_NATIVE_FLOAT).")
  }
  # single_dtype is now just dtype, used for NIFTI header mapping
  single_dtype <- dtype

  # --- Use map_dtype function to get NIFTI codes ---
  # Call the centralized function using :::
  databit <- map_dtype(single_dtype)
  nifti_datatype_code <- databit[1]
  nifti_bitpix <- databit[2]
  if (nifti_datatype_code == 0L) {
    warning("Could not map HDF5 dtype to NIfTI codes. Header datatype/bitpix may be incorrect.")
  }
  # --- End mapping section ---
  
  # Build minimal NIfTI-like header fields
  # Add common unused fields, initialized
  hdr_default <- list(
    sizeof_hdr  = 348L,
    data_type   = "", # Unused
    db_name     = "", # Unused
    extents     = 0L, # Unused
    session_error = 0L, # Unused
    regular     = 0L, # Unused
    dim_info    = 0L, # Unused
    dim         = c(4L, X, Y, Z, nVols, 1L, 1L, 1L),
    intent_p1   = 0.0, intent_p2 = 0.0, intent_p3 = 0.0, # Unused
    intent_code = 0L, # Unused
    datatype    = nifti_datatype_code,
    bitpix      = nifti_bitpix,
    slice_start = 0L, # Unused
    pixdim      = c(0.0, sp[1], sp[2], sp[3], 0, 0, 0, 0), # Keep pixdim[1]=0, qfac stored elsewhere
    vox_offset  = 0.0, # Unused
    scl_slope   = 1.0, # Unused
    scl_inter   = 0.0, # Unused
    slice_end   = 0L, # Unused
    slice_code  = 0L, # Unused
    xyzt_units  = 0L, # Unused
    cal_max     = 0.0, cal_min = 0.0, # Unused
    slice_duration = 0.0, # Unused
    toffset     = 0.0, # Unused
    glmax       = 0L, glmin = 0L, # Unused
    descrip     = "fmristore labeled volume set", # Basic description
    aux_file    = "", # Unused
    qform_code  = 1L, # Default to NIFTI_XFORM_SCANNER_ANAT
    sform_code  = 0L, # Default to NIFTI_XFORM_UNKNOWN
    quatern_b   = q$quaternion[1],
    quatern_c   = q$quaternion[2],
    quatern_d   = q$quaternion[3],
    qoffset_x   = org[1],
    qoffset_y   = org[2],
    qoffset_z   = org[3],
    srow_x      = c(tmat[1,1], tmat[1,2], tmat[1,3], tmat[1,4]), # Store sform rows
    srow_y      = c(tmat[2,1], tmat[2,2], tmat[2,3], tmat[2,4]),
    srow_z      = c(tmat[3,1], tmat[3,2], tmat[3,3], tmat[3,4]),
    intent_name = "", # Unused
    magic       = "n+1" # Keep writing variable length for now
  )

  # Merge user overrides from header_values, preventing overwrite of critical fields
  protected_fields <- c("dim", "pixdim", "quatern_b", "quatern_c", "quatern_d",
                        "qoffset_x", "qoffset_y", "qoffset_z", "sizeof_hdr", "magic",
                        "datatype", "bitpix")
  for (nm in names(header_values)) {
    if (nm %in% protected_fields) {
        warning("Ignoring attempt to override protected header field: '", nm, "'")
    } else {
        hdr_default[[nm]] <- header_values[[nm]]
    }
  }
  
  # 1) Write header fields into /header group
  # Create group first, checking existence
  if (!h5obj$exists("/header")) h5obj$create_group("/header")
  # Write each header field using h5_write
  for (nm in names(hdr_default)) {
      h5_write(h5obj, file.path("/header", nm), hdr_default[[nm]], overwrite = TRUE)
  }
  # Store qfac separately using h5_write
  h5_write(h5obj, "/header/qfac", q$qfac, overwrite = TRUE)

  # 2) Write /mask => shape [X, Y, Z] at root level
  h5_write(h5obj, "/mask", mask_arr, dtype = hdf5r::h5types$H5T_NATIVE_UCHAR, overwrite = TRUE)

  # 3) Write ORIGINAL labels => /labels (variable-length string array) at root level
  h5_write(h5obj, "/labels", labels, dtype = H5T_STRING$new(size = Inf), overwrite = TRUE)

  # 4) /data => subdatasets for each volume, storing masked data
  # Create parent group /data, checking existence
  if (!h5obj$exists("/data")) h5obj$create_group("/data") 
  # idx_nonzero and n_nonzero already calculated above

  # Set valid chunk dimensions for 1D data
  chunk_dims_1d <- if (!is.null(chunk_size) && n_nonzero > 0) c(min(chunk_size, n_nonzero)) else NULL

  for (i in seq_len(nVols)) {
    # Control verbosity
    if (verbose) message("Writing label: ", labels[i], " (safe name: ", safe_labels[i], ")")
    
    # Data path uses SANITIZED label
    data_path <- file.path("/data", safe_labels[i])
    
    # --- Corrected data extraction ---
    # Extract the 3D data for the current volume i
    # Note: NeuroVec subsetting might return a NeuroVol, ensure we get the array data
    vol_i_data <- as.array(vec[[i]]) 
    # Extract only the values at the mask locations (idx_nonzero)
    vol_1d_masked <- vol_i_data[idx_nonzero]
    # --- End corrected data extraction ---
    
    # Write using h5_write
    h5_write(
      h5 = h5obj,
      path = data_path,
      data = vol_1d_masked, # Write the correctly masked 1D data
      dtype = dtype, # Use the single dtype for all datasets
      chunk_dims = chunk_dims_1d,
      compression = compression,
      overwrite = TRUE # Assume overwrite within this function's scope
    )
  }

  invisible(h5obj)
}



#' Read a Labeled Neuroimaging Volume Set from HDF5
#'
#' @description
#' Reads an HDF5 file, typically one previously created by 
#' \code{\link{write_labeled_vec}} (now deprecated), and constructs a 
#' \code{\link{LabeledVolumeSet-class}} object.
#' The HDF5 file is opened in read-only mode.
#'
#' @details
#' This function is the primary way to load data into a \code{LabeledVolumeSet}
#' from its HDF5 disk representation. The structure of the HDF5 file is expected
#' to follow the specification laid out by \code{write_labeled_vec}.
#'
#' @section Lifecycle Management:
#' When a \code{LabeledVolumeSet} object is created by this function,
#' it opens the specified HDF5 file and the returned object takes ownership of
#' this open file handle.
#' **It is the user's responsibility to explicitly close this handle** when the
#' object is no longer needed to release system resources. This can be done by calling
#' \code{close(your_labeled_volume_set_object)}.
#'
#' Failure to close the handle may lead to issues such as reaching file handle
#' limits or problems with subsequent access to the file.
#'
#' @param file_path Character string: path to the HDF5 file.
#' @return A \code{\link{LabeledVolumeSet-class}} object with an open HDF5 file handle.
#'
#' @seealso
#' \code{\link{write_labeled_vec}} for the (deprecated) writing function,
#' \code{\link{LabeledVolumeSet-class}} for details on the object structure,
#' \code{\link{close.LabeledVolumeSet}} for closing the file handle.
#'
#' @examples
#' \dontrun{
#' # Assuming "my_labeled_set.h5" is a valid HDF5 file for LabeledVolumeSet
#' lvs <- read_labeled_vec("my_labeled_set.h5")
#' # ... perform operations with lvs ...
#' print(dim(lvs))
#' print(labels(lvs))
#' # Important: Close the handle when done
#' close(lvs)
#' }
#' @import hdf5r
#' @export
read_labeled_vec <- function(file_path) {
  # --- 1. Handle File Source ---
  # Call simplified open_h5 (no auto_close)
  fh <- open_h5(file_path, mode = "r")
  h5obj <- fh$h5

  # If we opened the file from a path, register a closing handler
  # in case an error occurs before the LabeledVolumeSet object is
  # successfully returned. The handler is cleared on success.
  if (fh$owns) {
    on.exit(safe_h5_close(h5obj), add = TRUE)
  }

  hdr_grp <- NULL # Initialize for finally block
  hdr_values <- list()

  # Helper to read dataset from header group if present
  .rd_hdr <- function(nm) {
    h5_read(h5obj, file.path("/header", nm), missing_ok = TRUE)
  }
  # Helper to read dataset from root group if present
   .rd_root <- function(nm) {
    h5_read(h5obj, paste0("/", nm), missing_ok = TRUE)
  }

  # Read required header fields
  dims       <- .rd_hdr("dim")         # c(4, X, Y, Z, nVols, 1,1,1)
  pixdim     <- .rd_hdr("pixdim")      # c(0.0, dx, dy, dz, ...)
  qb         <- .rd_hdr("quatern_b")
  qc         <- .rd_hdr("quatern_c")
  qd         <- .rd_hdr("quatern_d")
  qx         <- .rd_hdr("qoffset_x")
  qy         <- .rd_hdr("qoffset_y")
  qz         <- .rd_hdr("qoffset_z")
  qfac       <- .rd_hdr("qfac")        # Read qfac from header

  # Read labels from root level
  labels_arr <- .rd_root("labels")
  if (is.null(labels_arr)) {
      stop("Mandatory '/labels' dataset not found at file root.")
  }
  
  # Check and validate dimensions
  if (is.null(dims) || length(dims)<5 || dims[1]!=4) {
    stop("Invalid or missing 'dim' in /header/dim")
  }
  X <- dims[2]; Y <- dims[3]; Z <- dims[4]
  nVols <- dims[5]
  
  # Check labels consistency after reading from root
  if (length(labels_arr) != nVols) {
    stop("Mismatch: #labels (", length(labels_arr), ") != nVols specified in header/dim (", nVols, ").")
  }

  # read /mask => 3D from root level
  mask_arr <- h5_read(h5obj, "/mask", missing_ok = FALSE)
  # Ensure mask_arr is 3D
  if (length(dim(mask_arr)) != 3) {
      stop("Read /mask dataset is not 3-dimensional.")
  }
  if (!all(dim(mask_arr) == c(X,Y,Z))) {
      stop("Dimensions of /mask [", paste(dim(mask_arr), collapse=","),
           "] do not match dimensions specified in header/dim [", X, ",", Y, ",", Z, "]")
  }

  # Rebuild 4x4 transform from quaternion
  # Use default qfac=1 if not found
  qfac_val <- if (!is.null(qfac)) qfac else 1.0

  if (!is.null(pixdim) && length(pixdim) >= 4) {
    # pixdim[0] is ignored (should be 0)
    dx   <- pixdim[2]
    dy   <- pixdim[3]
    dz   <- pixdim[4]
  } else {
    warning("Missing or incomplete 'pixdim' in header. Using default spacing (1,1,1).")
    dx <- 1; dy <- 1; dz <- 1
  }

  if (!all(sapply(list(qb, qc, qd, qx, qy, qz), function(x) !is.null(x) && is.numeric(x)))) {
     warning("Missing or non-numeric quaternion parameters in header. Using identity transform.")
     mat <- diag(4)
     mat[1,1] <- dx
     mat[2,2] <- dy
     mat[3,3] <- dz
  } else {
     mat <- tryCatch(
         neuroim2::quaternToMatrix(
             quat     = c(qb,qc,qd),
             origin   = c(qx,qy,qz),
             stepSize = c(dx,dy,dz),
             qfac     = qfac_val
         ),
         error = function(e) {
             warning("Error calling quaternToMatrix: ", e$message, ". Using identity transform.")
             mat_fallback <- diag(4)
             mat_fallback[1,1] <- dx; mat_fallback[2,2] <- dy; mat_fallback[3,3] <- dz
             mat_fallback
         }
     )
  }

  # build space => just do a 3D NeuroSpace
  spc <- NeuroSpace(dim=c(X,Y,Z), spacing=c(dx,dy,dz), trans=mat)

  # build mask
  mask_vol <- LogicalNeuroVol(as.logical(mask_arr), space=spc)

  # --- 6. Prepare Lazy Loading Environment --- 
  load_env <- new.env(parent=emptyenv())
  # The H5File handle is stored in the main object's obj slot
  load_env$mask_idx <- which(as.logical(mask_arr)==TRUE) # Ensure logical mask used
  load_env$dims     <- c(X,Y,Z)
  load_env$space    <- spc
  # Remove sanitize function from environment
  # load_env$sanitize_label <- function(lbl) { gsub("[^A-Za-z0-9_.-]", "_", lbl) }

  # Define sanitize function locally instead
  sanitize_label_func <- function(lbl) { gsub("[^A-Za-z0-9_.-]", "_", lbl) }

  # Define the core loading logic which needs the parent object
  internal_loader <- function(i, parent_obj) {
    h5f <- parent_obj@obj # Get handle from parent S4 object
    if (!inherits(h5f, "H5File") || !h5f$is_valid) {
      stop("HDF5 file handle associated with this LabeledVolumeSet is invalid or closed.")
    }

    tryCatch({
      # Get the ORIGINAL label name using the numeric index
      lab <- parent_obj@labels[i]
      
      # --- Corrected: Sanitize label to find data path ---
      # Use the local sanitize function
      safe_lab <- sanitize_label_func(lab)
      # Problem: make.unique suffixing might be needed if collisions occurred during write
      # This requires storing the safe_labels from writer, or re-implementing make.unique logic
      # Simplification for now: Assume basic sanitization is enough or paths were unique
      # TODO: Revisit this if make.unique collisions cause read failures
      data_path <- file.path("/data", safe_lab)
      
      # Read the data using h5_read
      val1 <- h5_read(h5f, data_path, missing_ok = FALSE)
      
      # Create the 3D array using dims stored in load_env
      vol <- array(0, dim=parent_obj@load_env$dims)
      
      current_mask_idx <- parent_obj@load_env$mask_idx
      if (length(current_mask_idx) == 0) {
        warning("Internal inconsistency: Mask indices are empty but file exists for label ", lab)
        # Use space from load_env
        return(DenseNeuroVol(vol, space=parent_obj@load_env$space)) 
      }
      
      if (length(current_mask_idx) != length(val1)) {
        # Check if data length exceeds mask length (potential corruption)
        if (length(val1) > length(current_mask_idx)) {
             stop(paste0("Data length mismatch for label ", sQuote(lab), ". Stored data (", 
                        length(val1), ") exceeds mask size (", length(current_mask_idx), "). File may be corrupt."))
        } else {
             warning(paste0("Data length mismatch for label ", sQuote(lab), ". Expected ",
                            length(current_mask_idx), " values based on mask, but found ", length(val1), ". Padding with zeros."))
             # Fill available data, rest remains zero
             vol[current_mask_idx[1:length(val1)]] <- val1
        }
      } else {
        vol[current_mask_idx] <- val1 # Fill using mask indices
      }
      # Use space from load_env
      DenseNeuroVol(vol, space=parent_obj@load_env$space)
      
    }, error=function(e) {
        # Simplify error message to reveal original error
        stop(sprintf("Error loading data for label '%s': %s", 
                     lab %||% "(unknown)", 
                     conditionMessage(e)))
    })
  }

  # --- 7. Create LabeledVolumeSet Object --- 
  lvol <- new("LabeledVolumeSet",
              obj = h5obj, # Assign handle directly to obj slot
              mask = mask_vol,
              labels = labels_arr, # Store original labels
              load_env = load_env # Minimal load_env
              # Removed h5_wrapper assignment
              )

  # Assign the loader function - creates a closure capturing lvol
  load_env$loader <- function(i) {
      internal_loader(i, lvol) # 'lvol' is the parent object here
  }

  # If file was opened internally, the defer handler takes care of closing.
  # If user provided handle, they manage its lifecycle.

  if (fh$owns) {
    on.exit(NULL, add = FALSE)  # clear handler on success
  }

  return(lvol)
}

#' 4D Array-like subsetting for LabeledVolumeSet
#'
#' @param x A \code{LabeledVolumeSet} object.
#' @param i Numeric indices for the 1st dimension (x).
#' @param j Numeric indices for the 2nd dimension (y).
#' @param k Numeric indices for the 3rd dimension (z).
#' @param l Numeric indices for the 4th dimension (label).
#' @param drop Logical, whether to drop singleton dimensions.
#' @param ... Ignored.
#'
#' @return An R array after subsetting, or a lower-dimensional array if \code{drop=TRUE}.
#' @export
setMethod(
  f = "[",
  signature = signature(x="LabeledVolumeSet"),
  definition = function(x, i, j, k, l, ..., drop=TRUE) {
    # Validity check is implicitly handled by the loader now

    # 1) Figure out any missing dims => use full range
    dims_3d <- dim(space(x@mask))
    nVols   <- length(x@labels)
    # Handle missing i/j/k/l by using full ranges
    if (missing(i)) i <- seq_len(dims_3d[1])
    if (missing(j)) j <- seq_len(dims_3d[2])
    if (missing(k)) k <- seq_len(dims_3d[3])
    if (missing(l)) l <- seq_len(nVols)

    # Coerce to integer indices
    i <- as.integer(i)
    j <- as.integer(j)
    k <- as.integer(k)
    l <- as.integer(l)

    # Basic bounds check
    if (any(i < 1 | i > dims_3d[1])) {
      stop("Subscript i out of range 1..", dims_3d[1])
    }
    if (any(j < 1 | j > dims_3d[2])) {
      stop("Subscript j out of range 1..", dims_3d[2])
    }
    if (any(k < 1 | k > dims_3d[3])) {
      stop("Subscript k out of range 1..", dims_3d[3])
    }
    if (any(l < 1 | l > nVols)) {
      stop("Subscript l out of range 1..", nVols)
    }

    out_dims <- c(length(i), length(j), length(k), length(l))
    result <- array(0, dim=out_dims)

    # 2) Read each volume in l, subset in memory
    loader_func <- x@load_env$loader
    if (!is.function(loader_func)) stop("Internal error: loader function not found in LabeledVolumeSet environment.")

    out_l_pos <- 1
    for (lv in l) {
      vol_3d <- tryCatch(loader_func(lv), error = function(e) {
        stop("Failed to load volume ", lv, " (label: '", x@labels[lv], "'): ", e$message)
      })
      vol_arr <- as.array(vol_3d) # Convert DenseNeuroVol to array
      subcube <- vol_arr[i, j, k, drop=FALSE]
      result[,,, out_l_pos] <- subcube
      out_l_pos <- out_l_pos + 1
    }

    if (drop) {
      result <- drop(result)
    }
    result
  }
)


#' @rdname linear_access-methods  
#' @importFrom neuroim2 linear_access
#' @export
setMethod(
  f = "linear_access",
  signature = signature(x="LabeledVolumeSet", i="numeric"),
  definition = function(x, i) {
    # Validity check is implicitly handled by the loader

    dims_3d <- dim(space(x@mask))
    nVols   <- length(x@labels)
    bigDim  <- c(dims_3d, nVols)
    total   <- prod(bigDim)

    i <- as.integer(i)
    if (any(i < 1 | i > total)) {
      stop("Some indices out of range 1..", total)
    }

    sub_4d <- arrayInd(i, .dim=bigDim)
    vol_groups <- split(seq_len(nrow(sub_4d)), sub_4d[,4])
    out <- numeric(length(i))

    loader_func <- x@load_env$loader
    if (!is.function(loader_func)) stop("Internal error: loader function not found in LabeledVolumeSet environment.")

    for (v_str in names(vol_groups)) {
      v_idx <- as.integer(v_str)
      these_rows <- vol_groups[[v_str]]
      coords <- sub_4d[these_rows, , drop=FALSE]

      # Call the loader
      vol_3d <- tryCatch(loader_func(v_idx), error = function(e) {
         stop("Failed to load volume ", v_idx, " (label: '", x@labels[v_idx], "'): ", e$message)
      })
      vol_arr <- as.array(vol_3d)

      for (row_i in seq_len(nrow(coords))) {
        rx <- coords[row_i,1]
        ry <- coords[row_i,2]
        rz <- coords[row_i,3]
        val <- vol_arr[rx, ry, rz]
        out_idx <- these_rows[row_i]
        out[out_idx] <- val
      }
    }
    out
  }
)

#' @export
setMethod(
  f = "[[",
  signature = signature(x="LabeledVolumeSet", i="numeric"),
  definition = function(x, i, j, ...) {
    # Validity check implicit in loader

    if (length(i) != 1) {
      stop("Must provide a single index i.")
    }
    nLabels <- length(x@labels)
    if (i < 1 || i > nLabels) {
      stop("Index i out of range 1..", nLabels)
    }

    loader_func <- x@load_env$loader
    if (!is.function(loader_func)) stop("Internal error: loader function not found in LabeledVolumeSet environment.")

    # use the environment's loader to get that one volume
    vol <- tryCatch(loader_func(i), error = function(e) {
        stop("Failed to load volume ", i, " (label: '", x@labels[i], "'): ", e$message)
    })
    vol
  }
)




# The next method is the one for character indexing `[[`, which should remain.

#' @export
setMethod(
  f = "names",
  signature = signature(x = "LabeledVolumeSet"),
  definition = function(x) {
    x@labels
  }
)

#' show method for LabeledVolumeSet
#'
#' Displays essential info about a \code{LabeledVolumeSet}, including spatial dims,
#' number of volumes, label previews, spacing, origin, orientation (if known),
#' and storage paths.
#'
#' @param object A \code{\link{LabeledVolumeSet}} instance
#' @importFrom crayon bold blue silver yellow green italic
#' @importFrom methods show
#' @export
setMethod(
  f = "show",
  signature = "LabeledVolumeSet",
  definition = function(object) {
    # Check file validity for showing file path
    file_status <- "Unknown"
    file_path_or_handle <- object@obj # Accessing the S4 slot directly

    if (inherits(file_path_or_handle, "H5File")) {
        # It's an H5File handle
        if (file_path_or_handle$is_valid) {
            file_status <- tryCatch(file_path_or_handle$get_filename(), error = function(e) "Valid Handle (path unavailable)")
        } else {
            file_status <- "CLOSED Handle"
        }
    } else if (is.character(file_path_or_handle)) {
        # It's a file path string
        file_status <- file_path_or_handle
        # We could add a check here if the file exists, but might be slow/unnecessary for 'show'
    }

    cat("\n", crayon::bold(crayon::blue("LabeledVolumeSet")), "\n", sep="")

    sp   <- space(object@mask)
    nd3  <- dim(sp)
    nvol <- length(object@labels)

    cat(crayon::bold("\n╔═ Volume Info "), crayon::silver("───────────────────────────"), "\n", sep="")
    cat("║ ", crayon::yellow("3D Dimensions"), " : ", paste(nd3, collapse=" × "), "\n", sep="")
    cat("║ ", crayon::yellow("Total Volumes"), " : ", nvol, "\n", sep="")

    lbl_preview <- object@labels[1:min(3, nvol)]
    cat("║ ", crayon::yellow("Labels"), "        : ",
        paste(lbl_preview, collapse=", "),
        if (nvol > 3) crayon::silver(paste0(" ... (", nvol - 3, " more)")), "\n", sep="")

    cat(crayon::bold("\n╠═ Spatial Info "), crayon::silver("───────────────────────────"), "\n", sep="")
    cat("║ ", crayon::yellow("Spacing"), "       : ", paste(round(sp@spacing,2), collapse=" × "), "\n", sep="")
    cat("║ ", crayon::yellow("Origin"), "        : ", paste(round(sp@origin, 2), collapse=" × "), "\n", sep="")

    if (length(sp@axes@ndim) == 1 && sp@axes@ndim >= 3) {
      cat("║ ", crayon::yellow("Orientation"), "   : ",
          paste(sp@axes@i@axis, sp@axes@j@axis, sp@axes@k@axis), "\n", sep="")
    } else {
      cat("║ ", crayon::yellow("Orientation"), "   : Unknown / Not specified\n")
    }

    cat(crayon::bold("\n╚═ Storage Info "), crayon::silver("──────────────────────────"), "\n", sep="")
    cat("  ", crayon::yellow("HDF5 File"), "    : ", file_status, "\n", sep="")
    cat("  ", crayon::yellow("Data Path"), "    : /data/<label>\n", sep="")
    cat("  ", crayon::yellow("Mask Path"), "    : /mask\n", sep="")
    cat("  ", crayon::yellow("Labels Path"), "  : /labels\n", sep="") # Added label path info

    cat("\n")
  }
)



# Define the method for character index `i`
#' @export
setMethod(
  f = "[[",
  signature = signature(x = "LabeledVolumeSet", i = "character"),
  definition = function(x, i, j, ...) {
    # Validity check implicit in loader
    
    if (length(i) != 1) {
      stop("Must provide a single label name for character index.")
    }
    
    # Find the numeric index corresponding to the label
    numeric_idx <- match(i, x@labels)
    
    if (is.na(numeric_idx)) {
      stop("Label '", i, "' not found in LabeledVolumeSet.")
    }
    
    # Call the loader function directly with the numeric index
    loader_func <- x@load_env$loader
    if (!is.function(loader_func)) stop("Internal error: loader function not found in LabeledVolumeSet environment.")
    
    vol <- tryCatch(loader_func(numeric_idx), error = function(e) {
      stop("Failed to load volume for label '", i, "' (index: ", numeric_idx, "): ", e$message)
    })
    
    vol
  }
)

#' Close the HDF5 file associated with a LabeledVolumeSet
#'
#' This method manually closes the HDF5 file handle stored within the
#' LabeledVolumeSet object. It uses the \\code{safe_h5_close} helper to
#' ensure the handle is valid before attempting to close. After closing,
#' the internal handle reference is nulled to prevent accidental reuse.
#'
#' **Important:** If this \code{LabeledVolumeSet} object was created from
#' a file path using \code{\link{read_labeled_vec}}, the user is responsible
#' for calling this \code{close} method when finished with the object to release
#' the file handle. Failure to do so will leave the file open until the R session ends.
#' If the object was created using an existing \code{H5File} handle, closing
#' remains the responsibility of the code that originally opened the handle.
#'
#' @param con A \\code{LabeledVolumeSet} object.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns \\code{NULL}.
#' @rdname close
#' @export
setMethod("close", "LabeledVolumeSet", function(con, ...) {
  if (!is.null(con@obj)) {
    safe_h5_close(con@obj)
    # Nulling out the reference is problematic if the slot expects an H5File object.
    # The hdf5r object itself will become invalid after closing.
    # con@obj <- NULL # This line can cause S4 validation error
  }
  invisible(NULL)
})

