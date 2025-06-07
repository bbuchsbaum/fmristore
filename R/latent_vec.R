#' @importFrom neuroim2 matricized_access concat axes indices space origin spacing trans
#' @importFrom Matrix tcrossprod
NULL

#' Latent Space Representation of Neuroimaging Data
#'
#' @name LatentNeuroVec-class
#' @description
#' The \code{LatentNeuroVec} class provides a memory-efficient representation of neuroimaging
#' data in a latent space (e.g., from PCA or ICA). It stores data as a product of basis vectors
#' and loadings, allowing for efficient storage and computation of high-dimensional
#' neuroimaging data.
#'
#' @section Implementation Details:
#' The class implements a matrix factorization approach where the data is represented as:
#' \deqn{X = B \times L^T + c}
#' where:
#' \itemize{
#'   \item B is the basis matrix (\eqn{n \times k})
#'   \item L is the loadings matrix (\eqn{p \times k})
#'   \item c is an optional offset vector
#'   \item n is the number of time points
#'   \item p is the number of voxels
#'   \item k is the number of components
#' }
#'
#' @seealso
#' \code{\link[neuroim2]{NeuroVec-class}} for the base 4D brain image class.
#' \code{\link[neuroim2]{AbstractSparseNeuroVec-class}} for the sparse representation framework.
NULL

#' Internal Constructor for LatentNeuroVecSource
#'
#' @description
#' Creates a new LatentNeuroVecSource object for handling file-backed storage.
#'
#' @param file_name Character string specifying the path to the HDF5 file.
#' @return A new \code{LatentNeuroVecSource} object.
#' @keywords internal
#' @noRd
LatentNeuroVecSource <- function(file_name) {
  if (!is.character(file_name) || length(file_name) != 1) {
    stop("'file_name' must be a single character string")
  }
  if (!file.exists(file_name)) {
    stop("File '", file_name, "' does not exist")
  }
  new("LatentNeuroVecSource", file_name=file_name)
}

#' Construct a LatentNeuroVec Object
#'
#' @title Create a Latent Space Representation of Neuroimaging Data
#' @description
#' Constructs a \code{\link{LatentNeuroVec-class}} object, which provides a memory-efficient
#' representation of neuroimaging data using matrix factorization. This is particularly useful
#' for dimensionality reduction techniques (e.g., PCA or ICA).
#'
#' @param basis A numeric or \code{Matrix} object (\eqn{n \times k}) containing the temporal basis.
#' @param loadings A numeric or \code{Matrix} object (\eqn{p \times k}) containing spatial loadings.
#' @param space A \code{\link[neuroim2]{NeuroSpace-class}} defining the spatial/temporal dims.
#' @param mask A \code{\link[neuroim2]{LogicalNeuroVol-class}} defining the brain mask.
#' @param offset Optional numeric vector of length \eqn{p} (voxel-wise offsets).
#' @param label Optional character label for the object.
#'
#' @return A new \code{\link{LatentNeuroVec-class}} instance.
#'
#' @examples
#' \dontrun{
#' library(Matrix)
#' library(neuroim2)
#'
#' # Example data
#' n_timepoints <- 100
#' n_components <- 10
#' n_voxels <- 1000
#'
#' # Create basis & loadings
#' basis <- Matrix(rnorm(n_timepoints * n_components),
#'                 nrow = n_timepoints,
#'                 ncol = n_components)
#' loadings <- Matrix(rnorm(n_voxels * n_components),
#'                    nrow = n_voxels,
#'                    ncol = n_components,
#'                    sparse = TRUE)
#'
#' # Create space (10x10x10 volume, 100 timepoints)
#' spc <- NeuroSpace(c(10,10,10,n_timepoints))
#'
#' # Create mask
#' mask_array <- array(TRUE, dim=c(10,10,10))
#' mask_vol   <- LogicalNeuroVol(mask_array, NeuroSpace(c(10,10,10)))
#'
#' # Construct LatentNeuroVec
#' lvec <- LatentNeuroVec(basis = basis,
#'                        loadings = loadings,
#'                        space = spc,
#'                        mask = mask_vol)
#' }
#'
#' @export
LatentNeuroVec <- function(basis, loadings, space, mask, offset = NULL, label = "") {
  # Validate 'space'
  if (!inherits(space, "NeuroSpace")) {
    stop("'space' must be a NeuroSpace object")
  }
  # Validate 'basis' / 'loadings'
  if (!is.matrix(basis) && !inherits(basis, "Matrix")) {
    stop("'basis' must be a matrix or Matrix object")
  }
  if (!is.matrix(loadings) && !inherits(loadings, "Matrix")) {
    stop("'loadings' must be a matrix or Matrix object")
  }

  # Ensure we have a LogicalNeuroVol for mask
  # Initialize space_for_map using the main 3D space first
  space_3d <- drop_dim(space)
  space_for_map <- space_3d 

  
  if (!inherits(mask, "LogicalNeuroVol")) {
    # Construct 3D space from the main 4D space
    # space_3d already created above
    mask <- LogicalNeuroVol(as.logical(mask), space_3d)
  } else {
    # If mask is already a LogicalNeuroVol, ensure its space matches the 3D part of the main space
    # space_3d already created above
    mask_space <- neuroim2::space(mask)
    # Use check_same_dims for the dimension check
    check_same_dims(mask_space, space_3d, 
                    msg = "Space dimensions of provided mask do not match dimensions derived from the main 4D space.")
    # Keep the check for space equality separate for now
    if (!isTRUE(all.equal(mask_space, space_3d))) { 
        stop("Space object of provided mask does not match the space derived from the main 4D space. Cannot create IndexLookupVol.")
    }
  }

  cardinality <- sum(mask)

  # Handle offset: treat NULL or empty numeric offset as no offset
  if (is.null(offset) || length(offset) == 0) {
    # No offset provided or empty offset, use numeric(0)
    offset <- numeric(0)
  } else {
    # Check that provided offset matches loadings rows
    if (length(offset) != nrow(loadings)) {
      stop("'offset' length must match number of rows in 'loadings'")
    }
  }

  # Dimension checks - check columns first
  if (ncol(loadings) != ncol(basis)) {
    stop("'basis' and 'loadings' must have the same number of columns")
  }
  if (nrow(basis) != dim(space)[4]) {
    stop("'basis' must have ", dim(space)[4], " rows (the 4th dimension of space)")
  }
  if (nrow(loadings) != cardinality) {
    stop(paste0("'loadings' must have ", cardinality, " rows (i.e. #non-zero in mask)"))
  }

  # Ensure all numeric inputs are finite
  if (!all(is.finite(basis))) {
    stop("'basis' must contain only finite values")
  }
  if (!all(is.finite(loadings))) {
    stop("'loadings' must contain only finite values")
  }
  if (length(offset) > 0 && !all(is.finite(offset))) {
    stop("'offset' must contain only finite values")
  }

  # Convert basis/loadings to Matrix objects, choosing dense/sparse based on density
  if (is.matrix(basis) && !is(basis, "Matrix")) {
      density_basis <- sum(basis != 0) / length(basis)
      if (density_basis > 0.5) {
          message("Input 'basis' is dense (", round(density_basis * 100), "% non-zero); storing as dense dgeMatrix.")
          basis <- Matrix::Matrix(basis, sparse = FALSE)
      } else {
          basis <- Matrix::Matrix(basis)
      }
  } # else: Already a Matrix object
  
  if (is.matrix(loadings) && !is(loadings, "Matrix")) {
      density_loadings <- sum(loadings != 0) / length(loadings)
      if (density_loadings > 0.5) {
           message("Input 'loadings' is dense (", round(density_loadings * 100), "% non-zero); storing as dense dgeMatrix.")
           loadings <- Matrix::Matrix(loadings, sparse = FALSE)
      } else {
          loadings <- Matrix::Matrix(loadings)
      }
  } # else: Already a Matrix object

  # Check component count
  k <- ncol(basis)
  if (k < 1) {
    stop("Number of components (columns in basis and loadings) must be >= 1")
  }

  # Create the object
  new("LatentNeuroVec",
      basis    = basis,
      loadings = loadings,
      space    = space,
      mask     = mask,
      map      = IndexLookupVol(space_for_map, as.integer(which(mask@.Data))), # Use determined space
      offset   = offset,
      label    = label)
}

#' Write LatentNeuroVec to HDF5 File
#'
#' @description
#' Writes a \code{LatentNeuroVec} to an HDF5 file with optional compression.
#'
#' @param x A \code{LatentNeuroVec} to write.
#' @param file_name \code{character} file path to the output HDF5.
#' @param nbit \code{logical}; if TRUE, uses N-bit compression (default: FALSE).
#' @param compression \code{integer} in [1..9] specifying compression level (default: 9).
#'
#' @return Invisible \code{NULL}, called for side effects (writes to disk).
#'
#' @details
#' This function saves:
#' \itemize{
#'   \item \code{basis} matrix
#'   \item \code{loadings} matrix
#'   \item \code{offset} vector
#'   \item Spatial metadata
#'   \item Mask information
#' }
#' all inside an HDF5 file for future loading.
#'
#' @examples
#' \dontrun{
#' if (requireNamespace("neuroim2", quietly = TRUE) &&
#'     requireNamespace("hdf5r", quietly = TRUE) &&
#'     requireNamespace("Matrix", quietly = TRUE) && 
#'     exists("write_vec", where = "package:neuroim2") &&
#'     !is.null(fmristore:::create_minimal_LatentNeuroVec)) {
#'   
#'   lnv <- NULL
#'   temp_h5_file <- NULL
#'   
#'   tryCatch({
#'     # Create a minimal LatentNeuroVec
#'     lnv <- fmristore:::create_minimal_LatentNeuroVec(
#'       space_dims = c(3L, 3L, 2L),
#'       n_time = 5L,
#'       n_comp = 2L
#'     )
#'     
#'     # Create a temporary file for output
#'     temp_h5_file <- tempfile(fileext = ".h5")
#'     
#'     # Write the LatentNeuroVec to HDF5
#'     write_vec(lnv, temp_h5_file, compression = 6)
#'     
#'     # Verify file was created
#'     if (file.exists(temp_h5_file)) {
#'       cat("Successfully wrote LatentNeuroVec to:", temp_h5_file, "\n")
#'     }
#'     
#'   }, error = function(e) {
#'     message("write_vec example failed: ", e$message)
#'   }, finally = {
#'     # Clean up temporary file
#'     if (!is.null(temp_h5_file) && file.exists(temp_h5_file)) {
#'       unlink(temp_h5_file)
#'     }
#'   })
#' }
#' }
#'
#' @seealso
#' \code{\link[neuroim2]{read_vec}} for reading back the file.
#'
#' @importFrom neuroim2 write_vec
#' @rdname write_vec-methods
#' @export
setMethod(
  f = "write_vec",
  signature = signature(x="LatentNeuroVec", file_name="character", format="missing", data_type="missing"),
  definition = function(x, file_name, compression=9, nbit=FALSE) {
    if (!is.character(file_name) || length(file_name) != 1) {
      stop("'file_name' must be a single character string")
    }
    if (!is.numeric(compression) || compression < 1 || compression > 9) {
      stop("'compression' must be an integer between 1 and 9")
    }

    # Call internal writer. It handles opening/closing the file via on.exit.
    h5obj <- to_h5_latentvec(
      vec        = x,
      file_name  = file_name,
      compression= compression,
      nbit = nbit
    )
    # obj$close() # REMOVED: to_h5_latentvec manages closure via on.exit
    
    invisible(NULL)
  }
)

#' Internal Matrix-Based Access for LatentNeuroVec
#'
#' @description
#' Internal method providing efficient matrix-based access to elements.
#'
#' @param x A \code{LatentNeuroVec}.
#' @param i Index specification.
#' @return Computed values.
#'
#' @keywords internal
#' @noRd
setMethod(
  f = "matricized_access",
  signature = signature(x="LatentNeuroVec", i="matrix"),
  definition = function(x, i) {
    if (!is.numeric(i) || ncol(i) != 2L)
      stop("`i` must be a numeric matrix with 2 columns (time, spatial-index)")

    ## -- 1. split and sanity-check the two index columns -------------------
    t_idx <- as.integer(i[, 1L])                       # time rows in @basis
    s_idx <- as.integer(i[, 2L])                       # spatial indices 1..X·Y·Z

    nt  <- nrow(x@basis)
    nxy <- prod(dim(x)[1:3])

    if (any(t_idx < 1L | t_idx > nt))
      stop("time index out of bounds")
    if (any(s_idx < 1L | s_idx > nxy))
      stop("spatial index out of bounds")

    ## -- 2. map spatial -> mask rows  (0 means 'outside the mask') ---------
    v_idx <- lookup(x@map, s_idx)                      # 0 / 1..nVox

    inside <- v_idx > 0L
    out    <- numeric(length(t_idx))                   # zeros by default

    if (any(inside)) {

      ## -- 3. gather the relevant rows  --------------------------------------
      b1 <- x@basis   [ t_idx[inside] , , drop = FALSE ]   # n_inside × k
      b2 <- x@loadings[ v_idx[inside], , drop = FALSE ]    # n_inside × k

      ## -- 4. pair-wise dot product + offset  --------------------------------
      # Handle potential NA from sparse matrix multiplication
      # Ensure we have matrices for rowSums - handle Matrix objects properly
      if (is.vector(b1) || (inherits(b1, "Matrix") && prod(dim(b1)) == length(b1))) {
        b1 <- as.matrix(b1)
        if (is.vector(b1)) b1 <- matrix(b1, nrow = 1)
      }
      if (is.vector(b2) || (inherits(b2, "Matrix") && prod(dim(b2)) == length(b2))) {
        b2 <- as.matrix(b2)
        if (is.vector(b2)) b2 <- matrix(b2, nrow = 1)
      }
      # Ensure both are regular matrices for element-wise multiplication
      if (inherits(b1, "Matrix")) b1 <- as.matrix(b1)
      if (inherits(b2, "Matrix")) b2 <- as.matrix(b2)
      dot_products <- rowSums(b1 * b2)
      # Replace NA with 0 (happens when all elements in a row are 0)
      dot_products[is.na(dot_products)] <- 0
      
      if (length(x@offset) > 0) {
        out[inside] <- dot_products + x@offset[v_idx[inside]]
      } else {
        out[inside] <- dot_products
      }
    }

    out
  }
)

#' @keywords internal
#' @noRd
#' @importFrom neuroim2 matricized_access
#' @importFrom Matrix tcrossprod
setMethod(
  f = "matricized_access",
  signature = signature(x="LatentNeuroVec", i="integer"),
  definition = function(x, i) {
    # Import functions needed
    tcrossprod <- Matrix::tcrossprod
    
    # Ensure component dimensions match
    stopifnot("[matricized_access,LatentNeuroVec,integer] Number of components mismatch." = 
              ncol(x@basis) == ncol(x@loadings))
              
    if (any(i < 1) || any(i > nrow(x@loadings))) {
      stop("Index out of bounds for 'loadings'")
    }
    b1 <- x@basis
    b2 <- x@loadings[as.integer(i),, drop=FALSE] # Use as.integer explicitly
    
    # Convert to regular matrices to avoid dispatch issues
    if (inherits(b1, "Matrix")) b1 <- as.matrix(b1)
    if (inherits(b2, "Matrix")) b2 <- as.matrix(b2)
    
    # Use regular tcrossprod
    out <- b1 %*% t(b2)
    # Add offsets
    if (length(x@offset) > 0) {
      as.matrix(sweep(out, 2, x@offset[as.integer(i)], "+"))
    } else {
      as.matrix(out)
    }
  }
)

#' @keywords internal
#' @noRd
setMethod(
  f = "matricized_access",
  signature = signature(x="LatentNeuroVec", i="numeric"),
  definition = function(x, i) {
    # Coerce to integer and call the integer method using callNextMethod
    callNextMethod(x = x, i = as.integer(i)) 
  }
)

#' @export 
#' @rdname linear_access-methods
setMethod(
  f = "linear_access",
  signature = signature(x="LatentNeuroVec", i="numeric"),
  definition = function(x, i) {
    # Coerce to integer and call the integer method using callNextMethod
    callNextMethod(x = x, i = as.integer(i)) 
  }
)

#' @export 
#' @rdname linear_access-methods
setMethod(
  f = "linear_access",
  signature = signature(x="LatentNeuroVec", i="integer"),
  definition = function(x, i) {
    
    dims_full <- dim(x) # Get 4D dimensions
    nels_4d <- prod(dims_full)
    nels_3d <- prod(dims_full[1:3])
    n_time <- dims_full[4]
    
    # 1. Validate 4D linear indices
    if (!is.numeric(i) || any(is.na(i))) {
      stop("[linear_access,LatentNeuroVec] Index `i`` must be numeric without NA values")
    }
    if (any(i < 1) || any(i > nels_4d)) {
      stop(paste0("[linear_access,LatentNeuroVec] Index out of bounds for 4D volume [1..", nels_4d, "]"))
    }
    
    # 2. Convert 4D linear indices to 3D spatial index + time index
    # Time index (1-based)
    time_idx <- ceiling(i / nels_3d) 
    # 3D spatial linear index (1-based)
    spatial_idx_3d <- i %% nels_3d
    spatial_idx_3d[spatial_idx_3d == 0] <- nels_3d

    # 3. Map 3D spatial indices to mask indices (rows in loadings/offset)
    rowmap <- lookup(x@map, spatial_idx_3d) # Length = length(i)
    
    # 4. Identify unique needed mask indices and time indices
    unique_valid_mask_idx <- unique(rowmap[rowmap > 0])
    unique_time_idx <- unique(time_idx)
    
    # 5. If no requested indices fall within the mask, return zeros
    if (length(unique_valid_mask_idx) == 0) {
      return(numeric(length(i))) # Return vector of zeros
    }

    # --- Optimization Start ---
    # Pre-calculate transposed loadings (k x nVox_mask) - used multiple times
    t_loadings <- t(x@loadings) 
    basis_subset <- x@basis[unique_time_idx, , drop = FALSE] # [nUniqueTime, k]
    
    # Shortcut for single time point query
    if (length(unique_time_idx) == 1L) {
      # Compute values only for valid mask indices for this single time point
      # basis_subset [1,k] %*% t_loadings[k, unique_valid_mask_idx] => [1, nUniqueValidMask]
      computed_vals <- drop(basis_subset %*% t_loadings[, unique_valid_mask_idx, drop = FALSE])
      # Add offset if it exists
      if (length(x@offset) > 0) {
          computed_vals <- computed_vals + x@offset[unique_valid_mask_idx]
      }
      
      # Map these computed values back to the original query indices
      ovals <- numeric(length(i)) # Initialize with zeros
      # Find which original indices correspond to the computed valid mask indices
      original_indices_for_computed <- match(rowmap, unique_valid_mask_idx)
      # Select only those that were actually computed (non-NA in the match)
      valid_output_positions <- !is.na(original_indices_for_computed)
      # Place computed values into the correct spots in the output vector
      ovals[valid_output_positions] <- computed_vals[original_indices_for_computed[valid_output_positions]]
      return(ovals)
    }
    # --- Optimization End ---

    # 6. Calculate the required data block for multiple time points
    # Result is [length(unique_time_idx), length(unique_valid_mask_idx)]
    # Use pre-calculated t_loadings and basis_subset
    data_block <- basis_subset %*% t_loadings[, unique_valid_mask_idx, drop = FALSE]
    
    # 7. Add offsets if they exist
    if (length(x@offset) > 0) {
        data_block <- sweep(data_block, 2, x@offset[unique_valid_mask_idx], "+")
    }
    
    # 8. Create the output vector
    ovals <- numeric(length(i))
    
    # 9. Map results back to the original requested indices
    #    Need to find the correct row (time) and column (mask index) in data_block
    #    for each element of the original index 'i'.
    
    # Create lookup maps for row/column indices in data_block
    time_map <- match(time_idx, unique_time_idx)
    mask_map <- match(rowmap, unique_valid_mask_idx) # NA for out-of-mask
    
    # Identify indices that were *in* the mask
    in_mask_selector <- which(rowmap > 0)
    
    # Extract the corresponding row and column indices *within data_block*
    row_indices_in_block <- time_map[in_mask_selector]
    col_indices_in_block <- mask_map[in_mask_selector]
    
    # Convert 2D indices (row, col) in data_block to linear indices
    linear_indices_in_block <- row_indices_in_block + (col_indices_in_block - 1) * nrow(data_block)
    
    # Assign values from data_block to the output vector
    ovals[in_mask_selector] <- data_block[linear_indices_in_block]
    
    return(ovals)
  }
)

#' Extract a Single Volume from LatentNeuroVec
#'
#' @description
#' Extracts a single volume from a \code{LatentNeuroVec} as a \code{SparseNeuroVol}.
#'
#' @param x A \code{\link{LatentNeuroVec-class}} object.
#' @param i A numeric index specifying which volume to extract (must be a single value).
#'
#' @return A \code{\link[neuroim2]{SparseNeuroVol-class}} containing:
#' \itemize{
#'   \item The computed volume data
#'   \item A new 3D \code{NeuroSpace} object
#'   \item The original spatial indices
#' }
#'
#' @details
#' 1. Validates index \code{i}
#' 2. Computes the volume data with matrix operations
#' 3. Adds the offset
#' 4. Builds a new 3D \code{NeuroSpace}
#' 5. Returns a \code{SparseNeuroVol}
#'
#' @examples
#' \dontrun{
#' # Extract volumes
#' vol1 <- lvec[[1]]
#' vol_mid <- lvec[[dim(lvec)[4] / 2]]
#' vol_last <- lvec[[dim(lvec)[4]]]
#' }
#'
#' @seealso
#' \code{\link[neuroim2]{SparseNeuroVol-class}} for the return type,
#' \code{\link[neuroim2]{NeuroSpace-class}} for spatial metadata.
#'
#' @importFrom neuroim2 SparseNeuroVol
#' @rdname extract-methods
#' @export
# single volume at time i
setMethod("[[", signature(x="LatentNeuroVec", i="numeric"),
          function(x, i) {
            if (length(i) != 1) stop("Index must be a single number")
            if (i < 1 || i > dim(x)[4]) stop("Index out of range")

            # basis[i,,drop=FALSE] => shape (1 x k)
            # loadings => shape (p x k)
            # We want (1 x p).
            # Calculate using tcrossprod(B, L) for B %*% t(L)
            # tcrossprod(x@basis[i,,drop=FALSE], x@loadings) => (1 x k) %*% (k x p) => (1 x p)
            # Ensure both are Matrix objects for tcrossprod
            b1 <- x@basis[i,,drop=FALSE]
            b2 <- x@loadings
            # Matrix objects already - no need to convert
            dat <- as.numeric(b1 %*% t(b2))
            dat <- dat + x@offset  # length p

            # Now place them in a SparseNeuroVol with the known mask
            newdim <- dim(x)[1:3]
            bspace <- NeuroSpace(newdim,
                                 spacing=neuroim2::spacing(x),
                                 origin=neuroim2::origin(x),
                                 axes=neuroim2::axes(x@space),
                                 trans=neuroim2::trans(x))

            SparseNeuroVol(dat, bspace, indices=neuroim2::indices(x))
          }
)

#' @rdname extract-methods
#' @export
setMethod(
  f = "[",
  signature = signature(x="LatentNeuroVec", i="numeric", j="numeric", drop="ANY"),
  definition = function(x, i, j, k, l, ..., drop=TRUE) {

    # Determine full dimensions
    dims_full <- dim(x)
    if (missing(i)) i <- seq_len(dims_full[1])
    if (missing(j)) j <- seq_len(dims_full[2])
    if (missing(k)) k <- seq_len(dims_full[3])
    if (missing(l)) l <- seq_len(dims_full[4])

    # Basic bounds check
    if (any(i < 1 | i > dims_full[1]) || any(j < 1 | j > dims_full[2]) ||
        any(k < 1 | k > dims_full[3]) || any(l < 1 | l > dims_full[4])) {
      stop("Subscript out of range for LatentNeuroVec.")
    }

    i <- as.integer(i)
    j <- as.integer(j)
    k <- as.integer(k)
    l <- as.integer(l)

    # Calculate the output dimensions
    out_dim <- c(length(i), length(j), length(k), length(l))
    n_vox_req <- prod(out_dim[1:3]) # Number of voxels requested in the output slice
    n_time_req <- out_dim[4]

    # 1. Map requested (i, j, k) coordinates to 1D indices in the full 3D space
    #    Avoid `expand.grid(i, j, k)` which can allocate a very large
    #    intermediate data frame.  Instead compute linear indices in chunks
    #    with `k` as the outer loop and the i/j grid precomputed once.
    dims_space_3d <- dims_full[1:3]
    linear_idx_3d <- integer(n_vox_req)

    ij_grid <- outer(i, j, function(ii, jj) ii + (jj - 1L) * dims_space_3d[1])
    ij_step <- length(i) * length(j)

    for (kk in seq_along(k)) {
      start <- (kk - 1L) * ij_step + 1L
      end <- kk * ij_step
      linear_idx_3d[start:end] <-
        ij_grid + (k[kk] - 1L) * dims_space_3d[1] * dims_space_3d[2]
    }
    
    # 2. Map these 3D linear indices to rows in loadings/offset using the mask map
    #    `rowmap` will have 0 for out-of-mask voxels.
    rowmap <- lookup(x@map, linear_idx_3d) # Length = n_vox_req
    valid_mask_indices <- rowmap[rowmap > 0] # Indices within loadings/offset needed
    map_req_to_valid <- which(rowmap > 0)   # Mapping from requested voxel pos to valid pos

    # Initialize result array with zeros
    result <- array(0, dim = out_dim)

    # 3. If no requested voxels are within the mask, return the zero array
    if (length(valid_mask_indices) == 0) {
        if (drop) return(drop(result))
        return(result)
    }

    # 4. Calculate the requested time points for the *valid* voxels
    #    basis_subset (length(l) x k) %*% loadings[valid_mask_indices, k]^T 
    #    => result_valid (length(l) x length(valid_mask_indices))
    basis_subset <- x@basis[l, , drop = FALSE]
    # Use tcrossprod, subsetting loadings directly
    # Ensure both are Matrix objects for tcrossprod
    b1 <- basis_subset
    b2 <- x@loadings[valid_mask_indices, , drop = FALSE]
    # Matrix objects already - no need to convert
    result_valid <- b1 %*% t(b2)

    # 5. Add offset (only for the valid mask indices)
    #    Use sweep along the columns (MARGIN=2) of result_valid
    result_valid <- sweep(result_valid, 2, x@offset[valid_mask_indices], "+")

    # 6. Place the calculated values into the final result array using linear indexing
    #    `map_req_to_valid` gives the linear index (1..n_vox_req) within each 3D slice
    #    for each valid column in `result_valid`.
    #    The `result` array is arranged [x, y, z, t].
    #    Linear index formula: voxel_linear_index + (time_index - 1) * n_voxels_per_slice
    n_vox_per_slice <- prod(out_dim[1:3])
    for (col_idx in seq_along(valid_mask_indices)) {
      # Get the linear index for this voxel within a *single* 3D slice of the output
      voxel_linear_idx_in_slice <- map_req_to_valid[col_idx]
      # Calculate the linear indices for *all time points* for this specific voxel
      # result_valid[, col_idx] holds the time series for this voxel
      linear_indices_in_result <- voxel_linear_idx_in_slice + (seq_len(n_time_req) - 1L) * n_vox_per_slice
      # Assign the time series to the calculated linear indices in the result array
      result[linear_indices_in_result] <- result_valid[, col_idx]
    }

    # --- Remove old loop ---
    # for(t_idx in seq_len(n_time_req)) {
    #     # Create a temporary 3D slice filled with zeros
    #     slice_3d <- array(0, dim = out_dim[1:3])
    #     # Place the valid voxel values for this time point into the slice
    #     # The t_idx row of result_valid corresponds to the current time point
    #     # The columns correspond to valid_mask_indices, mapped by map_req_to_valid
    #     slice_3d[map_req_to_valid] <- result_valid[t_idx, ]
    #     # Assign the populated slice to the result array
    #     result[,,, t_idx] <- slice_3d
    # }
    
    # 7. Drop dimensions if requested
    if (drop) result <- drop(result)
    result
  }
)

#' Concatenate LatentNeuroVec Objects
#'
#' @description
#' Concatenates two or more \code{LatentNeuroVec} objects along the temporal dimension.
#'
#' @param x First \code{LatentNeuroVec}.
#' @param y Second \code{LatentNeuroVec}.
#' @param ... Additional \code{LatentNeuroVec} objects to concatenate.
#'
#' @return A \code{\link[neuroim2]{NeuroVecSeq-class}} containing all input vectors.
#'
#' @details
#' 1. Checks for compatibility
#' 2. Creates a sequence container
#' 3. Preserves individual object properties
#'
#' @examples
#' \dontrun{
#' # Concatenate multiple LatentNeuroVec objects
#' combined <- concat(lvec1, lvec2, lvec3)
#' dim(combined)
#' }
#'
#' @seealso
#' \code{\link[neuroim2]{NeuroVecSeq-class}}
#'
#' @importFrom neuroim2 NeuroVecSeq
#' @rdname concat-methods
#' @export
setMethod(
  f = "concat",
  signature = signature(x="LatentNeuroVec", y="LatentNeuroVec"),
  definition = function(x, y, ...) {
    # Get additional objects if any
    additional <- list(...)
    all_objects <- c(list(x, y), additional)
    all_lvecs <- all(sapply(all_objects, is, "LatentNeuroVec"))
    
    if (!all_lvecs) {
      # If not all objects are LatentNeuroVec, fall back to NeuroVecSeq
      return(do.call(NeuroVecSeq, list(x, y, ...)))
    }
    
    # Check if all objects have compatible spatial components
    # 1. Same dimensions (except time)
    # 2. Same mask
    # 3. Same column count in loadings (k value)
    
    # Check dimensions
    x_space <- space(x)
    x_dims_3d <- dim(x_space)[1:3]
    compatible_dims <- TRUE
    
    for (obj in all_objects[-1]) {  # Skip x
      obj_space <- space(obj)
      obj_dims_3d <- dim(obj_space)[1:3]
      # Use validate_same_dims and check if the result is NULL (success)
      validation_result <- validate_same_dims(x_dims_3d, obj_dims_3d, dims_to_compare = 1:3)
      if (!is.null(validation_result)) { # If result is NOT NULL, there was an error
        compatible_dims <- FALSE
        # Optionally, we could log the validation_result here
        # message("Dimension mismatch detected: ", validation_result)
        break
      }
    }
    
    if (!compatible_dims) {
      # Fall back to NeuroVecSeq if dimensions don't match
      return(do.call(NeuroVecSeq, list(x, y, ...)))
    }
    
    # Check masks
    x_mask_array <- as.array(mask(x))
    compatible_masks <- TRUE
    
    for (obj in all_objects[-1]) {  # Skip x
      obj_mask_array <- as.array(mask(obj))
      if (!identical(x_mask_array, obj_mask_array)) {
        compatible_masks <- FALSE
        break
      }
    }
    
    if (!compatible_masks) {
      # Fall back to NeuroVecSeq if masks don't match
      return(do.call(NeuroVecSeq, list(x, y, ...)))
    }
    
    # Check k values (component count)
    x_k <- ncol(x@loadings)
    compatible_k <- TRUE
    
    for (obj in all_objects[-1]) {  # Skip x
      obj_k <- ncol(obj@loadings)
      if (x_k != obj_k) {
        compatible_k <- FALSE
        break
      }
    }
    
    if (!compatible_k) {
      # Fall back to NeuroVecSeq if component counts don't match
      return(do.call(NeuroVecSeq, list(x, y, ...)))
    }
    
    # Check loadings - they must be identical for this optimization
    compatible_loadings <- TRUE
    x_loadings <- x@loadings
    
    for (obj in all_objects[-1]) {  # Skip x
      if (!identical(as.matrix(x_loadings), as.matrix(obj@loadings))) {
        compatible_loadings <- FALSE
        break
      }
    }
    
    if (!compatible_loadings) {
      # Fall back to NeuroVecSeq if loadings don't match
      return(do.call(NeuroVecSeq, list(x, y, ...)))
    }
    
    # If we get here, all spatial components are compatible
    # We can create a new LatentNeuroVec with concatenated basis matrices
    
    # Collect all basis matrices
    all_basis <- lapply(all_objects, function(obj) obj@basis)
    
    # Collect all time dimensions
    time_dims <- sapply(all_objects, function(obj) dim(obj@space)[4])
    total_time <- sum(time_dims)
    
    # Create new space with updated time dimension
    new_space_dims <- dim(x_space)
    new_space_dims[4] <- total_time
    new_space <- NeuroSpace(
      dim = new_space_dims,
      spacing = spacing(x_space),
      origin = origin(x_space),
      axes = axes(x_space),
      trans = trans(x_space)
    )
    
    # Concatenate basis matrices
    # We need to rbind them while preserving the Matrix class if possible
    if (all(sapply(all_basis, is, "sparseMatrix"))) {
      # Sparse case - use rbind from the Matrix package
      new_basis <- do.call(rbind, all_basis)
    } else {
      # Mixed or dense case - convert to matrices first
      all_basis_matrix <- lapply(all_basis, as.matrix)
      new_basis_matrix <- do.call(rbind, all_basis_matrix)
      # Convert back to Matrix object
      new_basis <- Matrix::Matrix(new_basis_matrix, sparse = (Matrix::nnzero(new_basis_matrix) / length(new_basis_matrix) < 0.5))
    }
    
    # Create a new LatentNeuroVec with concatenated basis
    new_label <- paste0(
      x@label,
      ifelse(nchar(x@label) > 0, "_plus_", ""),
      length(all_objects) - 1,
      "_more"
    )
    
    # Create new LatentNeuroVec
    LatentNeuroVec(
      basis = new_basis,
      loadings = x@loadings,  # Use loadings from the first object (all are identical)
      space = new_space,
      mask = x@mask,          # Use mask from the first object (all are identical)
      offset = x@offset,      # Use offset from the first object (should check if identical?)
      label = new_label
    )
  }
)

#' Convert LatentNeuroVec to HDF5 Format (Spec Compliant - Refactored)
#' @keywords internal
#' @noRd
to_h5_latentvec <- function(vec, file_name=NULL, data_type="FLOAT",
                            compression=6, nbit=FALSE) { 

  assert_that(inherits(vec, "LatentNeuroVec"))


  if (!is.null(file_name) && !endsWith(file_name, ".lv.h5")) {
    file_name <- paste0(file_name, ".lv.h5")
  } else if (is.null(file_name)) {
    file_name <- tempfile(fileext = ".lv.h5")
  }

  fh <- open_h5(file_name, mode = "w")
  h5obj <- fh$h5
  defer(if (fh$owns) h5obj$close_all(), envir = parent.frame())

  tryCatch({
      message("[to_h5_latentvec] Writing LatentNeuroVec to: ", file_name)

      hdf5r::h5attr(h5obj, "latent_spec_version") <- "1.0"
      hdf5r::h5attr(h5obj, "rtype") <- class(vec)
      hdf5r::h5attr(h5obj, "voxel_order") <- "column-major"

      sp <- neuroim2::space(vec)
      dims_vec <- dim(sp)
      if (length(dims_vec) != 4) stop("LatentNeuroVec space must be 4D.")
      X <- dims_vec[1]; Y <- dims_vec[2]; Z <- dims_vec[3]; T_vec <- dims_vec[4]

      tmat <- trans(sp)
      q <- tryCatch(neuroim2::matrixToQuatern(tmat),
                    error = function(e) {
                       warning("Failed to convert transform matrix to quaternions: ", e$message, ". Using default.")
                       list(quaternion=c(0,0,0), qoffset=c(0,0,0), qfac=1)
                    })
      sp_spacing <- spacing(sp)
      if (length(sp_spacing) < 3) sp_spacing <- c(sp_spacing, rep(1, 3 - length(sp_spacing)))
      TR <- attr(sp, "TR") %||% 0.0

      h5dtype_internal <- switch(toupper(data_type),
                        "FLOAT"   = hdf5r::h5types$H5T_NATIVE_FLOAT,
                        "DOUBLE"  = hdf5r::h5types$H5T_NATIVE_DOUBLE,
                        hdf5r::h5types$H5T_NATIVE_FLOAT)
      if (is.null(h5dtype_internal)) stop(paste0("Invalid data_type specified: ", data_type))

      dtype_text <- h5dtype_internal$to_text()
      nifti_dt_map <- list("H5T_NATIVE_FLOAT"=16L, "H5T_IEEE_F32LE"=16L,
                           "H5T_NATIVE_DOUBLE"=64L, "H5T_IEEE_F64LE"=64L,
                           "H5T_NATIVE_INT"=8L, "H5T_STD_I32LE"=8L,
                           "H5T_NATIVE_SHORT"=4L)
      nifti_bp_map <- list("H5T_NATIVE_FLOAT"=32L, "H5T_IEEE_F32LE"=32L,
                           "H5T_NATIVE_DOUBLE"=64L, "H5T_IEEE_F64LE"=64L,
                           "H5T_NATIVE_INT"=32L, "H5T_STD_I32LE"=32L,
                           "H5T_NATIVE_SHORT"=16L)
      nifti_datatype_code <- nifti_dt_map[[dtype_text]] %||% 0L
      nifti_bitpix <- nifti_bp_map[[dtype_text]] %||% 0L
      if (nifti_datatype_code == 0L) {
          warning("Could not map data_type ", data_type, " (", dtype_text, ") to NIfTI codes.")
      }

      hdr_fields <- list(
        sizeof_hdr = 348L, data_type = "", db_name = "", extents = 0L,
        session_error = 0L, regular = "", dim_info = 0L,
        dim = c(4L, X, Y, Z, T_vec, 1L, 1L, 1L),
        intent_p1 = 0.0, intent_p2 = 0.0, intent_p3 = 0.0, intent_code = 0L,
        datatype = nifti_datatype_code, bitpix = nifti_bitpix,
        slice_start = 0L,
        pixdim = c(q$qfac %||% 1.0, sp_spacing[1], sp_spacing[2], sp_spacing[3], TR, 0.0, 0.0, 0.0),
        vox_offset = 0.0, scl_slope = 1.0, scl_inter = 0.0,
        slice_end = as.integer(Z - 1), slice_code = 0L, xyzt_units = 10L,
        cal_max = 0.0, cal_min = 0.0, slice_duration = 0.0, toffset = 0.0,
        glmax = 0L, glmin = 0L,
        descrip = paste("LatentNeuroVec data:", vec@label %||% "(no label)"),
        aux_file = "", qform_code = 1L, sform_code = 0L,
        quatern_b = if(!is.null(q$quaternion)) q$quaternion[1] else 0.0,
        quatern_c = if(!is.null(q$quaternion)) q$quaternion[2] else 0.0,
        quatern_d = if(!is.null(q$quaternion)) q$quaternion[3] else 0.0,
        qoffset_x = if(!is.null(q$qoffset)) q$qoffset[1] else 0.0,
        qoffset_y = if(!is.null(q$qoffset)) q$qoffset[2] else 0.0,
        qoffset_z = if(!is.null(q$qoffset)) q$qoffset[3] else 0.0,
        srow_x = c(tmat[1,1], tmat[1,2], tmat[1,3], tmat[1,4]),
        srow_y = c(tmat[2,1], tmat[2,2], tmat[2,3], tmat[2,4]),
        srow_z = c(tmat[3,1], tmat[3,2], tmat[3,3], tmat[3,4]),
        intent_name = "", magic = "n+1"
      )

      .write_header(h5obj, hdr_fields, q$qfac %||% 1.0)

      nvox <- .write_mask(h5obj, vec@mask, compression)

      basis_info <- .write_basis(h5obj, vec@loadings, h5dtype_internal,
                                 compression, vec@offset, nvox)

      scan_name <- vec@label %||% tools::file_path_sans_ext(basename(file_name))
      .write_scans(h5obj, as.matrix(vec@basis), scan_name, TR,
                   basis_info$k, T_vec, h5dtype_internal, compression)

      message("[to_h5_latentvec] HDF5 write SUCCESSFUL.")
      return(h5obj)

  }, error = function(e) {
      warning(paste0("[to_h5_latentvec] ERROR during HDF5 write: ", e$message))
      return(NULL)
  })
}

#' @keywords internal
#' @noRd
.write_header <- function(h5, hdr_fields, qfac) {
  message("[to_h5_latentvec] Writing Header Group...")
  for (nm in names(hdr_fields)) {
    h5_write(h5, file.path("/header", nm), hdr_fields[[nm]],
             dtype = guess_h5_type(hdr_fields[[nm]]), overwrite = TRUE)
  }
  h5_write(h5, "/header/qfac", qfac,
           dtype = h5types$H5T_NATIVE_DOUBLE, overwrite = TRUE)
  message("[to_h5_latentvec] Header Group DONE.")
}

#' @keywords internal
#' @noRd
.write_mask <- function(h5, mask_vol, compression) {
  message("[to_h5_latentvec] Writing Mask Dataset...")
  if (!inherits(mask_vol, "LogicalNeuroVol"))
    stop("vec@mask is not a LogicalNeuroVol")
  mask_arr <- array(as.integer(mask_vol@.Data), dim = dim(mask_vol))
  dims <- dim(mask_vol)
  mask_chunk_dim <- c(min(32, dims[1]), min(32, dims[2]), min(32, dims[3]))
  h5_write(h5, "/mask", mask_arr,
           dtype = hdf5r::h5types$H5T_NATIVE_UCHAR,
           chunk_dims = mask_chunk_dim, compression = compression,
           overwrite = TRUE)
  nvox <- sum(mask_vol)

  message("[to_h5_latentvec] Writing Voxel Coords...")
  mask_indices <- which(mask_arr == 1L)
  coords_to_write <- if (length(mask_indices) > 0) {
    coords <- arrayInd(mask_indices, .dim = dim(mask_arr))
    matrix(as.integer(coords - 1), nrow = length(mask_indices), ncol = 3)
  } else {
    matrix(integer(), nrow = 0, ncol = 3)
  }
  if (nrow(coords_to_write) != nvox)
    stop("Internal dimension mismatch: voxel_coords rows != non-zero count in mask")
  coord_chunk <- if (nrow(coords_to_write) > 0) c(min(1024, nrow(coords_to_write)), 3) else NULL
  h5_write(h5, "/voxel_coords", coords_to_write,
           dtype = hdf5r::h5types$H5T_NATIVE_INT32,
           chunk_dims = coord_chunk, compression = compression,
           overwrite = TRUE)
  return(nvox)
}

#' @keywords internal
#' @noRd
.write_basis <- function(h5, loadings, h5dtype, compression, offset, nvox_mask) {
  message("[to_h5_latentvec] Writing Basis Group (spatial components)...")
  h5$create_group("/basis")
  t_loadings <- Matrix::t(loadings)
  k <- nrow(t_loadings)
  nvox <- ncol(t_loadings)
  if (nvox != nvox_mask)
    stop(paste0("Internal dimension mismatch: spatial basis columns (", nvox,
                ") != nVox in mask (", nvox_mask, ")"))
  density <- Matrix::nnzero(t_loadings) / length(t_loadings)
  write_sparse <- density < 0.30
  message(paste0("  Spatial basis density: ", round(density*100, 2),
                 "%. Writing as ", if(write_sparse) "SPARSE" else "DENSE", "."))
  if (!write_sparse) {
    dense_chunk <- c(k, min(1024, nvox))
    h5_write(h5, "/basis/basis_matrix", as.matrix(t_loadings),
             dtype = h5dtype,
             chunk_dims = dense_chunk, compression = compression,
             overwrite = TRUE)
  } else {
    sparse_grp <- h5$create_group("/basis/basis_matrix_sparse")
    hdf5r::h5attr(sparse_grp, "storage") <- "csc"
    hdf5r::h5attr(sparse_grp, "shape") <- dim(t_loadings)
    if (!inherits(t_loadings, "dgCMatrix"))
      t_loadings <- methods::as(t_loadings, "CsparseMatrix")
    nnz <- length(t_loadings@x)
    target_min_chunk <- 128 * 1024
    element_size <- h5dtype$get_size()
    chunk_len <- max(1L, floor(target_min_chunk / element_size))
    if (nnz > 0 && chunk_len > nnz) chunk_len <- nnz else if(nnz == 0) chunk_len <- NULL
    h5_write(h5, "/basis/basis_matrix_sparse/data", t_loadings@x, h5dtype,
             chunk_dims = chunk_len, compression = compression,
             overwrite = TRUE)
    h5_write(h5, "/basis/basis_matrix_sparse/indices", as.integer(t_loadings@i),
             hdf5r::h5types$H5T_NATIVE_INT32,
             chunk_dims = chunk_len, compression = compression,
             overwrite = TRUE)
    h5_write(h5, "/basis/basis_matrix_sparse/indptr", as.integer(t_loadings@p),
             hdf5r::h5types$H5T_NATIVE_INT32,
             chunk_dims = NULL, compression = 0,
             overwrite = TRUE)
  }
  message("[to_h5_latentvec] Basis Group DONE.")
  message("[to_h5_latentvec] Writing Offset...")
  if (length(offset) > 0) {
    if (length(offset) != nvox)
      stop(paste0("Offset length (", length(offset), ") does not match spatial basis nVox (", nvox, ")"))
    h5_write(h5, "/offset", offset, dtype = h5dtype,
             chunk_dims = NULL, compression = 0, overwrite = TRUE)
  } else {
    message("  Offset is empty, skipping write.")
  }
  invisible(list(k = k, nVox = nvox))
}

#' @keywords internal
#' @noRd
.write_scans <- function(h5, embedding_matrix, scan_name, TR, k, T_vec,
                         h5dtype, compression) {
  message("[to_h5_latentvec] Writing Scans Group...")
  h5$create_group("/scans")
  scan_name <- gsub("[^a-zA-Z0-9_.-]", "_", scan_name)
  if (!nzchar(scan_name)) scan_name <- "scan_1"

  meta_path <- file.path("/scans", scan_name, "metadata")
  h5$create_group(file.path("/scans", scan_name))
  h5$create_group(meta_path)
  run_length <- nrow(embedding_matrix)
  h5_write(h5, file.path(meta_path, "run_length"), as.integer(run_length), overwrite = TRUE)
  if (TR > 0) h5_write(h5, file.path(meta_path, "TR"), as.double(TR), overwrite = TRUE)
  h5_write(h5, file.path(meta_path, "subject_id"), "", dtype=hdf5r::H5T_STRING$new(size = Inf), overwrite = TRUE)
  h5_write(h5, file.path(meta_path, "task"), "", dtype=hdf5r::H5T_STRING$new(size = Inf), overwrite = TRUE)
  h5_write(h5, file.path(meta_path, "session"), "", dtype=hdf5r::H5T_STRING$new(size = Inf), overwrite = TRUE)

  message("  Writing embedding matrix...")
  embed_dims <- dim(embedding_matrix)
  if (embed_dims[1] != T_vec)
    warning(paste0("Embedding time points (", embed_dims[1], ") does not match header time dim (", T_vec, "). Using embedding dim."))
  if (embed_dims[2] != k)
    stop(paste0("Embedding components (", embed_dims[2], ") mismatch spatial basis components (", k, ")"))
  embed_chunk <- c(min(128, embed_dims[1]), k)
  h5_write(h5, file.path("/scans", scan_name, "embedding"), embedding_matrix,
           dtype = h5dtype, chunk_dims = embed_chunk,
           compression = compression, overwrite = TRUE)
  message("[to_h5_latentvec] Scans Group DONE.")
}
#' Internal helper to write each field from a list into an HDF5 group
#' @noRd
.write_nifti_header_fields <- function(hdr_grp, fields_list) {
  for (nm in names(fields_list)) {
    val <- fields_list[[nm]]
    if (is.null(val)) next
    .write_field(hdr_grp, nm, val)
  }
}

#' Internal helper for writing a single field
#' @noRd
.write_field <- function(hdr_grp, nm, val) {
  if (is.character(val)) {
    stype <- hdf5r::H5T_STRING$new(size = Inf)
    on.exit(stype$close(), add=TRUE)
    ds    <- hdr_grp$create_dataset(nm, dims=1, dtype=stype)
    on.exit(if (!is.null(ds) && ds$is_valid) ds$close(), add = TRUE)
    ds[]  <- if (length(val) > 1) paste(val, collapse=" ") else val
  } else if (is.integer(val)) {
    ds <- hdr_grp$create_dataset(nm, dims=length(val),
                                 dtype=hdf5r::h5types$H5T_NATIVE_INT)
    on.exit(if (!is.null(ds) && ds$is_valid) ds$close(), add = TRUE)
    ds[] <- as.integer(val)
  } else if (is.numeric(val)) {
    ds <- hdr_grp$create_dataset(nm, dims=length(val),
                                 dtype=hdf5r::h5types$H5T_NATIVE_DOUBLE)
    on.exit(if (!is.null(ds) && ds$is_valid) ds$close(), add = TRUE)
    ds[] <- as.double(val)
  } else {
    warning("Writing field '", nm, "' as string due to unrecognized type: ", class(val))
    stype <- hdf5r::H5T_STRING$new(size = Inf)
    on.exit(stype$close(), add=TRUE)
    ds    <- hdr_grp$create_dataset(nm, dims=1, dtype=stype)
    on.exit(if (!is.null(ds) && ds$is_valid) ds$close(), add = TRUE)
    ds[]  <- as.character(val)
  }
}

#' Load data from a LatentNeuroVecSource object (Spec Compliant)
#'
#' @description
#' Constructs a `LatentNeuroVec` from a spec-compliant HDF5 file.
#' Reads header, mask, basis, offset, and the embedding for a selected scan.
#'
#' @param x A `LatentNeuroVecSource` specifying the HDF5 file path.
#' @param scan_name Optional `character` string specifying which scan to load under `/scans/`. 
#'   If `NULL` (default), the first scan found will be loaded.
#'
#' @return A `LatentNeuroVec` object for the specified scan.
#'
#' @details
#' Reads data according to `BasisEmbeddingSpec.yaml`:
#' 1. Reads `/header` to build `NeuroSpace`.
#' 2. Reads `/mask` to build `LogicalNeuroVol`.
#' 3. Reads `/basis/basis_matrix` -> transposes -> `object@loadings`.
#' 4. Reads `/offset` (optional) -> `object@offset`.
#' 5. Reads `/scans/<scan_name>/embedding` -> `object@basis`.
#' 6. Performs dimension consistency checks.
#'
#' @note Requires `hdf5r`.
#' @seealso `LatentNeuroVecSource`, `LatentNeuroVec`, `validate_latent_file`
#'
#' @importFrom hdf5r H5File h5attr
#' @importFrom neuroim2 NeuroSpace NeuroVol LogicalNeuroVol IndexLookupVol drop_dim quaternToMatrix lookup
#' @noRd
setMethod(
  f = "load_data",
  signature = c("LatentNeuroVecSource"),
  definition = function(x, scan_name = NULL) {
    # --- 1. Handle File Source ---
    fh <- open_h5(x@file_name, mode = "r") # Use helper
    h5obj <- fh$h5
    # Defer closing only if we opened the file via path
    defer(if (fh$owns) h5obj$close_all(), envir = parent.frame())
        
    # --- TRY CATCH REMOVED FOR PARSING DIAGNOSIS --- 
    
    # --- 2. Check Version Attribute --- 
    version_attr <- NULL
    if (h5obj$attr_exists("latent_spec_version")) { 
        version_attr <- hdf5r::h5attr(h5obj, "latent_spec_version")
    }
    if (is.null(version_attr)) {
        warning("[load_data,LatentNeuroVecSource] HDF5 file '", x@file_name, 
                "' is missing the 'latent_spec_version' attribute. Assuming version 1.0 format.")
    } else if (version_attr != "1.0") {
        warning("[load_data,LatentNeuroVecSource] HDF5 file '", x@file_name, 
                "' has unexpected 'latent_spec_version'='", version_attr, 
                "'. Attempting to load assuming version 1.0 format.")
    }
        
    # --- 3. Read Header and Reconstruct Space --- 
    .rd_hdr <- function(nm, required=TRUE) {
      h5_read(h5obj, file.path("/header", nm), missing_ok = !required)
    }
    dims_hdr   <- .rd_hdr("dim", required=TRUE)
    pixdim_hdr <- .rd_hdr("pixdim", required=FALSE) # Allow missing, default later
    qb         <- .rd_hdr("quatern_b", required=FALSE)
    qc         <- .rd_hdr("quatern_c", required=FALSE)
    qd         <- .rd_hdr("quatern_d", required=FALSE)
    qx         <- .rd_hdr("qoffset_x", required=FALSE)
    qy         <- .rd_hdr("qoffset_y", required=FALSE)
    qz         <- .rd_hdr("qoffset_z", required=FALSE)
    qfac       <- .rd_hdr("qfac", required=FALSE) # Default to 1.0 if missing
    
    if (length(dims_hdr) < 5 || dims_hdr[1] != 4) {
        stop("Invalid '/header/dim' dimensions found.")
    }
    dims_3d <- dims_hdr[2:4]
    nTime_hdr <- dims_hdr[5]
    
    # Rebuild transform matrix
    # qfac_val <- qfac %||% 1.0 # Use utils::`%||%` or ifelse
    qfac_val <- if (is.null(qfac)) 1.0 else qfac
    spacing_3d <- if (!is.null(pixdim_hdr) && length(pixdim_hdr) >= 4) pixdim_hdr[2:4] else c(1,1,1)
    # origin_3d <- c(qx %||% 0, qy %||% 0, qz %||% 0)
    origin_3d <- c(if(is.null(qx)) 0 else qx, 
                   if(is.null(qy)) 0 else qy, 
                   if(is.null(qz)) 0 else qz)
    
    if (!all(sapply(list(qb, qc, qd), function(q) !is.null(q) && is.numeric(q)))) {
         warning("Missing or non-numeric quaternion b,c,d parameters. Using default orientation.")
         trans_mat <- diag(4)
         trans_mat[1,1] <- spacing_3d[1]; trans_mat[2,2] <- spacing_3d[2]; trans_mat[3,3] <- spacing_3d[3]
         trans_mat[1:3, 4] <- origin_3d
    } else {
         trans_mat <- tryCatch(
             neuroim2::quaternToMatrix(
                 quat     = c(qb, qc, qd),
                 origin   = origin_3d,
                 stepSize = spacing_3d,
                 qfac     = qfac_val
             ), error = function(e) {
                 warning("Error calling quaternToMatrix: ", e$message, ". Using default orientation.")
                 mat_fallback <- diag(4)
                 mat_fallback[1,1] <- spacing_3d[1]; mat_fallback[2,2] <- spacing_3d[2]; mat_fallback[3,3] <- spacing_3d[3]
                 mat_fallback[1:3, 4] <- origin_3d
                 mat_fallback
             })
    }
    
    # Create 4D NeuroSpace (needed for LatentNeuroVec)
    full_space <- NeuroSpace(dim = dims_hdr[2:5], # Use 4D dims + time
                             spacing = spacing_3d, # Use 3D spacing
                             origin = origin_3d,
                             trans = trans_mat)
    space_3d <- drop_dim(full_space) # 3D space for mask

    # --- 4. Read Mask --- 
    mask_arr <- h5_read(h5obj, "/mask", missing_ok = FALSE)
    # Check dimensions of mask array against header dims - check_same_dims stops on error
    check_same_dims(dim(mask_arr), dims_3d, dims_to_compare = 1:3, 
                      msg = paste("load_data: Dimensions of /mask mismatch header dims"))
    
    mask_vol <- LogicalNeuroVol(as.logical(mask_arr), space = space_3d)
    nVox_mask <- sum(mask_vol)

    # --- 5. Read Voxel Coords --- 
    voxel_coords <- h5_read(h5obj, "/voxel_coords", missing_ok = TRUE)
    if (!is.null(voxel_coords)) {
        if (length(dim(voxel_coords)) != 2 || dim(voxel_coords)[2] != 3) {
            warning("/voxel_coords dataset does not have dimensions [nVox, 3]. Discarding.")
            voxel_coords <- NULL
        } else if (nrow(voxel_coords) != nVox_mask) {
            warning(paste0("/voxel_coords rows (", nrow(voxel_coords),
                             ") mismatch non-zero count in mask (", nVox_mask, "). Discarding."))
            voxel_coords <- NULL
        }
    }
    if (is.null(voxel_coords) && nVox_mask > 0) {
        # Derive if missing or failed to read correctly
        warning("Deriving voxel coordinates from mask.")
        voxel_coords <- t(arrayInd(which(mask_vol@.Data), .dim = dim(mask_vol))) 
    } else if (is.null(voxel_coords) && nVox_mask == 0) {
        voxel_coords <- matrix(integer(), nrow=0, ncol=3) # Empty matrix if no voxels
    }
    # `voxel_coords` now holds the [nVox, 3] matrix, either read or derived.
    # Not strictly needed for LatentNeuroVec construction if mask & map are correct, but good practice to load.

    
    # --- 6. Read Basis Matrix -> object@loadings --- 
    basis_grp_exists <- h5obj$exists("/basis")
    if (!basis_grp_exists) stop("Mandatory group '/basis' not found.")
    internal_loadings <- NULL
    has_dense <- h5obj$exists("/basis/basis_matrix")
    has_sparse_grp <- h5obj$exists("/basis/basis_matrix_sparse")
    
    if (has_dense && has_sparse_grp) {
        stop("Found both dense '/basis/basis_matrix' and sparse '/basis/basis_matrix_sparse' group. File is invalid.")
    } else if (has_dense) {
        basis_matrix <- h5_read(h5obj, "/basis/basis_matrix", missing_ok = FALSE)
        if (length(dim(basis_matrix)) != 2) stop("Dense basis matrix is not 2-dimensional.")
        k_basis <- dim(basis_matrix)[1]
        nVox_basis <- dim(basis_matrix)[2]
        if (nVox_basis != nVox_mask) {
            stop("Dense basis matrix nVox (", nVox_basis, ") mismatch non-zero count in mask (", nVox_mask, ")")
        }
        # NOTE: Mapping HDF5 /basis/basis_matrix (spatial, k x p) -> R @loadings (spatial, p x k)
        internal_loadings <- Matrix::Matrix(t(basis_matrix))
        message("[load_data] Dense spatial basis loaded.")

    } else if (has_sparse_grp) {
        sparse_grp <- h5obj[["/basis/basis_matrix_sparse"]]
        on.exit(if (!is.null(sparse_grp) && sparse_grp$is_valid) sparse_grp$close(), add=TRUE) # Keep this specific on.exit
        # ... (read attributes) ...
        data_data <- h5_read(h5obj, "/basis/basis_matrix_sparse/data", missing_ok = FALSE)
        indices_data <- h5_read(h5obj, "/basis/basis_matrix_sparse/indices", missing_ok = FALSE)
        indptr_data <- h5_read(h5obj, "/basis/basis_matrix_sparse/indptr", missing_ok = FALSE)
        # ... (sparse matrix reconstruction assigns to internal_loadings) ...
        message("[load_data] Sparse triplet datasets read.")
        
        # Reconstruct k x p matrix
        reconstructed_k_p <- NULL
        if (storage_fmt == "csc") {
            expected_indptr_len <- shape_attr[2] + 1
            if (length(indptr_data) != expected_indptr_len) stop("CSC indptr length mismatch.")
            reconstructed_k_p <- Matrix::sparseMatrix(i = indices_data, p = indptr_data, x = data_data, 
                                                      dims = shape_attr, index1 = FALSE) # Assume 0-based indices from HDF5
        } else if (storage_fmt == "csr") {
            expected_indptr_len <- shape_attr[1] + 1
             if (length(indptr_data) != expected_indptr_len) stop("CSR indptr length mismatch.")
            # Reconstruct as if it were CSC (i=row_indices, p=row_pointers), then transpose
            # Spec: CSR stores column indices in 'indices' and row pointers in 'indptr'
            # Matrix::sparseMatrix needs row indices (i) and column pointers (p) for CSC
            # We have column indices (j) and row pointers (p_row). To build k x p using sparseMatrix:
            # This requires iterating through row pointers to determine row index for each element. 
            # SIMPLER: Build the p x k matrix directly using CSR components
            #          (Need j = col_indices, p = row_pointers)
            message("  Reconstructing CSR matrix (p x k) from triplet...")
            reconstructed_p_k <- Matrix::sparseMatrix(j = indices_data, 
                                                      p = indptr_data, 
                                                      x = data_data, 
                                                      dims = rev(shape_attr), # Use p x k dimensions
                                                      index1 = FALSE) # HDF5 indices are 0-based
            # Transpose is not needed here, as we built p x k directly
            # internal_loadings <- Matrix::t(reconstructed_p_k) # No! reconstructed_p_k is already p x k
            internal_loadings <- reconstructed_p_k # Use p x k directly
        } else {
            stop(paste0("Unsupported sparse storage format: ", storage_fmt))
        }
        
        message("[load_data] Sparse spatial basis loaded.") # No transpose needed if built correctly
        
    } else {
         stop("No basis matrix found in '/basis'.")
    }

    # --- 7. Read Offset (Optional) --- 
    offset_vec <- h5_read(h5obj, "/offset", missing_ok = TRUE)
    if (!is.null(offset_vec) && length(offset_vec) != nVox_basis) {
        warning("Offset length (", length(offset_vec), ") mismatch basis nVox (", nVox_basis, "). Using zero offset.")
        offset_vec <- NULL
    }
    if (is.null(offset_vec)) {
        offset_vec <- numeric(nVox_basis) # Default to zero offset
    }

    # --- 8. Select Scan and Read Embedding -> object@basis --- 
    scans_grp_exists <- h5obj$exists("/scans")
    if (!scans_grp_exists) stop("Mandatory group '/scans' not found.")
    scans_grp <- h5obj[["/scans"]]
    available_scans <- scans_grp$ls()$name
    if (length(available_scans) == 0) stop("No scans found under '/scans'.")
    
    target_scan_name <- scan_name
    if (is.null(target_scan_name)) {
        target_scan_name <- available_scans[1]
        if (length(available_scans) > 1) {
            warning("'scan_name' not specified, loading first scan found: ", target_scan_name)
        }
    } else if (!(target_scan_name %in% available_scans)) {
        stop("Requested scan_name '", target_scan_name, "' not found in file. Available: ", 
             paste(available_scans, collapse=", "))
    }
    
    # No need to open scan_grp, h5_read handles full path
    embedding_path <- file.path("/scans", target_scan_name, "embedding")
    # NOTE: Mapping HDF5 /scans/.../embedding (temporal, n x k) -> R @basis (temporal, n x k)
    internal_basis <- h5_read(h5obj, embedding_path, missing_ok = FALSE)
    if (length(dim(internal_basis)) != 2) stop("Embedding matrix is not 2-dimensional.")
    nTime_embed <- dim(internal_basis)[1]
    k_embed <- dim(internal_basis)[2]
    
    # Final consistency checks
    if (k_embed != k_basis) {
        stop("Component mismatch: embedding k (", k_embed, ") != basis k (", k_basis, ")")
    }
    if (nTime_embed != nTime_hdr) {
        stop("Time mismatch: embedding time (", nTime_embed, ") != header time (", nTime_hdr, ")")
    }

    # --- 9. Construct LatentNeuroVec --- 
    LatentNeuroVec(basis    = Matrix::Matrix(internal_basis), 
                   loadings = internal_loadings, 
                   space    = full_space, 
                   mask     = mask_vol, 
                   offset   = offset_vec,
                   label    = target_scan_name)
    # --- END OF REMOVED TRY CATCH BLOCK --- 
} # End definition function
) # End setMethod

#' Validate HDF5 File Against Latent Spec
#'
#' @description
#' Performs basic checks on an HDF5 file to verify if it conforms to the
#' essential structure and dimension consistency defined by the
#' `BasisEmbeddingSpec.yaml` specification.
#'
#' @param file_path Character string specifying the path to the HDF5 file.
#'
#' @return Logical `TRUE` if basic checks pass, `FALSE` otherwise. Issues
#'   warnings for inconsistencies found. Throws an error if the file cannot
#'   be opened or fundamental groups/datasets are missing.
#'
#' @details
#' Checks performed:
#' \itemize{
#'   \item File existence and HDF5 readability.
#'   \item Presence of mandatory groups: `/header`, `/basis`, `/scans`.
#'   \item Presence of mandatory datasets: `/header/dim`, `/mask`, `/basis/basis_matrix`,
#'         at least one scan group under `/scans`, and its `embedding` dataset.
#'   \item Dimension consistency:
#'     \itemize{
#'       \item `/header/dim[2:4]` vs `/mask` dimensions.
#'       \item `/basis/basis_matrix` columns vs number of non-zero elements in `/mask`.
#'       \item `/basis/basis_matrix` rows (k) vs embedding columns (k).
#'       \item `/header/dim[5]` (time) vs embedding rows (time).
#'       \item Optional: `/offset` length vs `/basis/basis_matrix` columns.
#'     }
#' }
#' Does NOT validate header field *values* extensively beyond dimensions, nor does it
#' check data types rigorously.
#'
#' @examples
#' \dontrun{
#' # Create a temporary latent neuroimaging HDF5 file for validation
#' temp_file <- tempfile(fileext = ".h5")
#' 
#' # Create minimal latent data and write to HDF5
#' lnv <- fmristore:::create_minimal_LatentNeuroVec()
#' 
#' # Write to HDF5 (using internal structure expected by validate_latent_file)
#' # Note: This is a simplified example - real files would use write_vec methods
#' h5f <- hdf5r::H5File$new(temp_file, mode = "w")
#' 
#' # Write required groups and datasets
#' header_grp <- h5f$create_group("header")
#' header_grp$create_dataset("dim", robj = c(4L, dim(lnv@mask), dim(lnv@basis)[1]))
#' 
#' basis_grp <- h5f$create_group("basis")
#' basis_grp$create_dataset("basis_matrix", robj = lnv@basis)
#' 
#' scans_grp <- h5f$create_group("scans")
#' h5f$create_dataset("mask", robj = as.array(lnv@mask))
#' 
#' h5f$close_all()
#' 
#' # Now validate the file
#' result <- validate_latent_file(temp_file)
#' 
#' if (result$is_valid) {
#'   message("File is valid!")
#' } else {
#'   message("File validation failed: ", result$error_message)
#' }
#' 
#' # Clean up
#' unlink(temp_file)
#' }
#' 
#' @importFrom hdf5r H5File h5attr
#' @export
validate_latent_file <- function(file_path) {
  if (!is.character(file_path) || length(file_path) != 1 || !nzchar(file_path)) {
    stop("[validate_latent_file] 'file_path' must be a single character string.")
  }
  if (!file.exists(file_path)) {
    stop("[validate_latent_file] File not found: '", file_path, "'")
  }

  is_valid <- TRUE
  error_message <- NULL
  fh <- NULL # Initialize fh for potential use in finally block

  tryCatch({
    # --- 1. Handle File Source ---
    fh <- open_h5(file_path, mode = "r") # Use helper
    h5obj <- fh$h5
    # Defer closing only if we opened the file via path
    # CRITICAL: Ensure this defer is correctly placed to run after tryCatch completes
    defer(if (!is.null(fh) && fh$owns) try(h5obj$close_all(), silent = TRUE), envir = parent.frame())

    # --- Check for mandatory groups/datasets using h5obj$exists ---
    required_groups <- c("header", "basis", "scans")
    required_datasets <- c("mask") # Root datasets
    required_header_datasets <- c("dim")
    # required_basis_datasets <- c("basis_matrix") # Or basis_matrix_sparse

    for (grp_name in required_groups) {
      if (!h5obj$exists(grp_name)) {
        stop(paste0("Mandatory group '/", grp_name, "' not found."))
      }
    }
    for (ds_name in required_datasets) {
      if (!h5obj$exists(ds_name)) {
        stop(paste0("Mandatory dataset '/", ds_name, "' not found."))
      }
    }
    # Check header datasets using full path
    for (ds_name in required_header_datasets) {
       hdr_ds_path <- file.path("/header", ds_name)
       if (!h5obj$exists(hdr_ds_path)) {
           stop(paste0("Mandatory dataset '", hdr_ds_path, "' not found."))
       }
    }
    # Check basis dataset/group using full path
    basis_dense_path <- "/basis/basis_matrix"
    basis_sparse_path <- "/basis/basis_matrix_sparse"
    has_dense_basis <- h5obj$exists(basis_dense_path)
    has_sparse_basis <- h5obj$exists(basis_sparse_path)
    if (!has_dense_basis && !has_sparse_basis) {
        stop("Mandatory dataset '/basis/basis_matrix' or group '/basis/basis_matrix_sparse' not found.")
    }
    if (has_dense_basis && has_sparse_basis) {
         stop("Found both dense '/basis/basis_matrix' and sparse '/basis/basis_matrix_sparse'. File is invalid.")
    }
    # Check scans and embedding
    if (!h5obj$exists("/scans")) stop("Mandatory group '/scans' not found.")
    # Need to list scans without opening scans_grp handle long-term
    scan_names <- tryCatch(h5obj[["scans"]]$ls()$name, finally = try(h5obj[["scans"]]$close(), silent=TRUE))
    if (length(scan_names) == 0) {
        stop("No scan subgroups found under '/scans'.")
    }
    first_scan_name <- scan_names[1]
    embed_path <- file.path("/scans", first_scan_name, "embedding")
    if (!h5obj$exists(embed_path)) {
        stop(paste0("Mandatory dataset '", embed_path, "' not found."))
    }

    # --- Read dimensions for consistency checks using h5_read or direct dim access ---
    header_dim <- h5_read(h5obj, "/header/dim", missing_ok = FALSE)
    mask_data <- h5_read(h5obj, "/mask", missing_ok = FALSE)
    mask_dim <- dim(mask_data)

    basis_dim <- NULL
    k_basis <- NA_integer_
    nVox_basis <- NA_integer_
    if (has_dense_basis) {
        # Read dims directly from dataset handle for efficiency
        basis_dset <- NULL
        tryCatch({
           basis_dset <- h5obj[[basis_dense_path]]
           basis_dim <- basis_dset$dims
        }, finally= if(!is.null(basis_dset) && basis_dset$is_valid) try(basis_dset$close(), silent=TRUE))
    } else { # has_sparse_basis
         # Read shape attribute
         sparse_grp <- NULL
         tryCatch({
            sparse_grp <- h5obj[[basis_sparse_path]]
            if (!"shape" %in% names(hdf5r::h5attributes(sparse_grp))) stop("Sparse basis group missing 'shape' attribute.")
            basis_dim <- hdf5r::h5attr(sparse_grp, "shape") # Should be [k, nVox_basis]
         }, finally= if(!is.null(sparse_grp) && sparse_grp$is_valid) try(sparse_grp$close(), silent=TRUE))
    }
    if (is.null(basis_dim) || length(basis_dim) != 2) stop("Failed to read valid 2D dimensions for basis.")
    k_basis <- basis_dim[1]
    nVox_basis <- basis_dim[2]

    # Read embedding dims directly
    embed_dset <- NULL
    embed_dim <- NULL
     tryCatch({
       embed_dset <- h5obj[[embed_path]]
       embed_dim <- embed_dset$dims # Should be [T_emb, k]
     }, finally= if(!is.null(embed_dset) && embed_dset$is_valid) try(embed_dset$close(), silent=TRUE))
    if (is.null(embed_dim) || length(embed_dim) != 2) stop("Failed to read valid 2D dimensions for embedding.")
    T_emb <- embed_dim[1]
    k_embed <- embed_dim[2]


    # --- Perform Consistency Checks ---
    valid_checks <- list()

    # Check 1: Header dim[0] should be 4
    if (length(header_dim) < 5 || header_dim[1] != 4) {
        valid_checks$hdr_dim0 <- paste0("'/header/dim' should start with 4, but starts with ", header_dim[1])
    }
    # Check 2: Header spatial dims vs Mask dims
    dim_check_msg <- validate_same_dims(header_dim[2:4], mask_dim, dims_to_compare = 1:3)
    if (!is.null(dim_check_msg)) {
        valid_checks$hdr_mask_dims <- paste0("'/header/dim[2:4]' (", paste(header_dim[2:4], collapse=","),
                                             ") mismatch '/mask' dims (", paste(mask_dim, collapse=","), ")")
    }
    # Check 3: Basis components (k) vs Embedding components (k)
    if (k_basis != k_embed) {
        valid_checks$k_mismatch <- paste0("Component dimension mismatch: basis k (", k_basis,
                                           ") != embedding k (", k_embed, ")")
    }
    # Check 4: Header time (T_hdr) vs Embedding time (T_emb)
    T_hdr <- header_dim[5]
    if (T_hdr != T_emb) {
        valid_checks$time_mismatch <- paste0("Time dimension mismatch: header dim[5] (", T_hdr,
                                             ") != embedding rows (", T_emb, ")")
    }
    # Check 5: Basis nVox vs Mask non-zero count
    nVox_mask <- sum(mask_data > 0) # Assuming mask is 0/1
    if (nVox_basis != nVox_mask) {
        valid_checks$basis_nvox_mask <- paste0("Basis nVox mismatch: basis nVox (", nVox_basis,
                                               ") != non-zero count in /mask (", nVox_mask, ")")
    }
    # Check 6: Optional Offset length
    offset_val <- h5_read(h5obj, "/offset", missing_ok=TRUE)
    if (!is.null(offset_val)) {
        offset_len <- length(offset_val)
        if (offset_len != nVox_basis) {
            valid_checks$offset_len <- paste0("Offset length mismatch: offset length (", offset_len,
                                              ") != basis nVox (", nVox_basis, ")")
        }
    } # else: offset doesn't exist, no check needed

    if (length(valid_checks) > 0) {
        is_valid <- FALSE
        warning("[validate_latent_file] Validation failed for '", file_path, "' with issues:\n",
                paste("  - ", unlist(valid_checks), collapse="\n"))
    }

  }, error = function(e) {
    is_valid <<- FALSE # Modify variable in parent env
    error_message <<- e$message # Capture error message without prefix
    # Ensure h5obj is closed by on.exit, which should still trigger
  })

  if (!is.null(error_message)) {
      # Re-throw the error to match test expectations
      stop(error_message)
  }

  return(is_valid)
}


#' Extract method for LatentNeuroVec
#' @rdname extract-methods
#' @importFrom neuroim2 lookup
#' @export
setMethod(
  f = "[",
  signature = signature(x="LatentNeuroVec", i="ANY", j="ANY", drop="ANY"),
  definition = function(x, i, j, k, l, ..., drop=TRUE) {
    # Get dimensions
    dims <- dim(x)
    
    # Handle missing indices
    if (missing(i)) i <- seq_len(dims[1])
    if (missing(j)) j <- seq_len(dims[2])
    if (missing(k)) k <- seq_len(dims[3])
    if (missing(l)) l <- seq_len(dims[4])
    
    # Convert to integer
    i <- as.integer(i)
    j <- as.integer(j)
    k <- as.integer(k)
    l <- as.integer(l)
    
    # Calculate output dimensions
    out_dims <- c(length(i), length(j), length(k), length(l))
    result <- array(0, dim = out_dims)
    
    # Calculate the linear indices for the 3D spatial locations
    spatial_indices <- array(0, dim = out_dims[1:3])
    idx <- 1
    for (kk in seq_along(k)) {
      for (jj in seq_along(j)) {
        for (ii in seq_along(i)) {
          # Calculate linear index in original 3D space
          lin_idx <- i[ii] + (j[jj] - 1) * dims[1] + (k[kk] - 1) * dims[1] * dims[2]
          spatial_indices[ii, jj, kk] <- lin_idx
        }
      }
    }
    
    # Map spatial indices to mask indices
    spatial_vec <- as.vector(spatial_indices)
    mask_indices <- lookup(x@map, spatial_vec)
    
    # Get unique valid mask indices
    valid_mask_idx <- unique(mask_indices[mask_indices > 0])
    
    if (length(valid_mask_idx) > 0) {
      # Extract the subset of basis and loadings we need
      basis_sub <- x@basis[l, , drop = FALSE]  # [n_time_sub, k]
      loadings_sub <- x@loadings[valid_mask_idx, , drop = FALSE]  # [n_valid_vox, k]
      
      # Calculate values: basis_sub %*% t(loadings_sub)
      # Result is [n_time_sub, n_valid_vox]
      values <- basis_sub %*% t(loadings_sub)
      
      # Add offset
      if (length(x@offset) > 0) {
        values <- sweep(values, 2, x@offset[valid_mask_idx], "+")
      }
      
      # Now map these values back to the output array
      for (t_idx in seq_along(l)) {
        slice_3d <- array(0, dim = out_dims[1:3])
        # Map values back using mask_indices
        for (v_idx in seq_along(valid_mask_idx)) {
          # Find which spatial positions map to this mask index
          positions <- which(mask_indices == valid_mask_idx[v_idx])
          if (length(positions) > 0) {
            slice_3d[positions] <- values[t_idx, v_idx]
          }
        }
        result[,,,t_idx] <- slice_3d
      }
    }
    
    if (drop) {
      result <- drop(result)
    }
    
    result
  }
)

#' Double bracket extract method for LatentNeuroVec
#' @rdname extract-methods
#' @importFrom neuroim2 NeuroSpace SparseNeuroVol spacing origin trans
#' @export
setMethod(
  f = "[[",
  signature = signature(x="LatentNeuroVec", i="numeric"),
  definition = function(x, i, ...) {
    if (length(i) != 1) {
      stop("[[ can only extract a single time point")
    }
    
    # Get the full 4D array for just this time point
    result <- x[, , , i, drop = FALSE]
    
    # Create a SparseNeuroVol from the 3D result
    dims_3d <- dim(x)[1:3]
    space_3d <- NeuroSpace(dim = dims_3d, 
                           spacing = spacing(space(x))[1:3],
                           origin = origin(space(x))[1:3],
                           trans = trans(space(x)))
    
    # Find non-zero voxels
    vol_3d <- drop(result)
    non_zero_idx <- which(vol_3d != 0)
    
    if (length(non_zero_idx) > 0) {
      SparseNeuroVol(data = vol_3d[non_zero_idx],
                     space = space_3d,
                     indices = non_zero_idx)
    } else {
      # Return empty SparseNeuroVol
      SparseNeuroVol(data = numeric(0),
                     space = space_3d,
                     indices = integer(0))
    }
  }
)

#' @export
#' @rdname series-methods  
#' @importFrom neuroim2 series lookup
setMethod(
  f = "series",
  signature = c(x="LatentNeuroVec", i="integer"),
  definition = function(x, i, j, k, ..., drop = TRUE) {
    # Handle missing arguments
    if (missing(j)) j <- NULL
    if (missing(k)) k <- NULL
    
    # Extract any remaining arguments from dots
    dots <- list(...)
    if ("drop" %in% names(dots)) drop <- dots$drop
    
    # Check if j and k were actually provided
    has_j <- !is.null(j)
    has_k <- !is.null(k)
    
    nTime  <- dim(x)[4]
    nels3d <- prod(dim(x)[1:3])

    # CASE A: user gave only i -> interpret as multiple 3D voxel indices
    if (!has_j && !has_k) {
      if (any(i < 1 | i > nels3d)) {
        stop("Some voxel index in 'i' is out of range [1..(X*Y*Z)].")
      }
      n_vox_req <- length(i) # Number of requested voxels

      # Map requested 3D indices (i) to rows in loadings/offset (mask indices)
      # rowmap will have 0 for voxels outside the mask
      rowmap <- lookup(x@map, i)
      valid_mask_indices <- rowmap[rowmap > 0] # Indices within loadings/offset that are needed
      map_req_to_valid <- which(rowmap > 0)   # Mapping from requested index pos to valid index pos

      if (length(valid_mask_indices) == 0) {
          # All requested voxels are outside the mask
          if (n_vox_req == 1) {
            return(numeric(nTime)) # Return vector if only 1 voxel requested
          } else {
            return(matrix(0, nrow = nTime, ncol = n_vox_req)) # Return matrix otherwise
          }
      }

      # Use tcrossprod: basis (nTime x k) %*% loadings[valid_mask_indices, k]^T => (nTime x length(valid_mask_indices))
      # Ensure both are Matrix objects for tcrossprod
      b1 <- x@basis
      b2 <- x@loadings[valid_mask_indices, , drop = FALSE]
      # Already using matrix multiplication, no need to convert
      valid_vox_series <- b1 %*% t(b2)

      # Add offset (only for the valid mask indices)
      # Need to transpose offset subset to broadcast correctly with sweep
      valid_vox_series <- sweep(valid_vox_series, 2, x@offset[valid_mask_indices], "+")

      # Create the final output matrix, filled with 0s
      out_mat <- matrix(0, nrow = nTime, ncol = n_vox_req)
      # Place the calculated series into the correct columns
      out_mat[, map_req_to_valid] <- as.numeric(valid_vox_series)

      if (drop && n_vox_req == 1) {
        return(drop(out_mat))
      } else {
        return(out_mat)
      }

    } else {
      # CASE B: user gave i,j,k => each must be length 1 => single voxel
      if (!(length(i) == 1 && length(j) == 1 && length(k) == 1)) {
        stop("series(x, i,j,k): i,j,k must each be a single integer for one voxel.")
      }
      # convert (i,j,k) => 3D linear index
      idx_1d <- i + (j-1)*dim(x)[1] + (k-1)*dim(x)[1]*dim(x)[2]
      if (idx_1d<1 || idx_1d>nels3d) {
        stop("Voxel subscript (i,j,k) out of range for LatentNeuroVec.")
      }
      # map => row in loadings
      mr <- lookup(x@map, idx_1d)
      # if out-of-mask => entire time series is 0
      if (mr < 1) {
        return(numeric(nTime))
      }

      # Use tcrossprod: basis (nTime x k) %*% loadings[mr, k]^T => (nTime x 1)
      # Ensure both are Matrix objects for tcrossprod
      b1 <- x@basis
      b2 <- x@loadings[mr, , drop = FALSE]
      # Matrix objects already - no need to convert
      out_vec <- b1 %*% t(b2)
      # Add the single offset value
      out_vec <- out_vec + x@offset[mr]

      if (drop) as.vector(out_vec) else out_vec
    }
  }
)

#' @export
#' @rdname series-methods
setMethod(
  f = "series",
  signature = signature(x="LatentNeuroVec", i="numeric"),
  definition = function(x, i, j, k, ..., drop=TRUE) {
    # Cast to integer and call the integer method directly
    i <- as.integer(i)
    if (!missing(j)) {
      j <- as.integer(j)
      if (!missing(k)) {
        k <- as.integer(k)
        series(x, i, j, k, ..., drop=drop)
      } else {
        series(x, i, j, ..., drop=drop)
      }
    } else {
      series(x, i, ..., drop=drop)
    }
  }
)

#' @export
#' @rdname series-methods
setMethod(
  f = "series",
  signature = signature(x="LatentNeuroVec"),
  definition = function(x, i, j, k, ..., drop=TRUE) {
    # Default method - handle when i is missing or of unknown type
    if (missing(i)) {
      stop("series requires at least one index argument 'i'")
    }
    # Convert all to integer if numeric
    if (is.numeric(i)) i <- as.integer(i)
    
    # Call the appropriate method based on what arguments are provided
    if (!missing(j)) {
      j <- as.integer(j)
      if (!missing(k)) {
        k <- as.integer(k)
        series(x, i, j, k, ..., drop=drop)
      } else {
        series(x, i, j, ..., drop=drop)
      }
    } else {
      series(x, i, ..., drop=drop)
    }
  }
)


#' @importFrom crayon bold red green blue yellow silver
#' @importFrom utils object.size
#' @rdname show-methods
#' @export
setMethod(
  f = "show",
  signature = "LatentNeuroVec",
  definition = function(object) {
    # Header
    cat("\n", crayon::bold(crayon::blue("LatentNeuroVec Object")), "\n")
    cat(crayon::silver("======================\n"))

    # Dimensions
    dims <- dim(object)
    spatial_dims <- paste(dims[1:3], collapse=" x ")
    cat("\n", crayon::yellow("Dimensions:"), "\n")
    cat(" ", crayon::silver("*"), " Spatial: ", crayon::green(spatial_dims), "\n")
    cat(" ", crayon::silver("*"), " Temporal: ", crayon::green(dims[4]), "\n")

    # Components
    n_components <- ncol(object@basis)
    first_basis_coeffs <- format(object@basis[1:min(5, nrow(object@basis)), 1], digits=3)
    cat("\n", crayon::yellow("Components:"), "\n")
    cat(" ", crayon::silver("*"), " Number: ", crayon::green(n_components), "\n")
    cat(" ", crayon::silver("*"), " First component, first 5 coeffs: ",
        crayon::green(paste(first_basis_coeffs, collapse=", ")), 
        if (length(first_basis_coeffs) < 5) "" else "...", "\n") # Adjusted description

    # Memory Usage
    basis_size    <- format(object.size(object@basis), units="auto")
    loadings_size <- format(object.size(object@loadings), units="auto")
    total_size    <- format(object.size(object),          units="auto")

    cat("\n", crayon::yellow("Memory Usage:"), "\n")
    cat(" ", crayon::silver("*"), " Basis: ",    crayon::green(basis_size),    "\n")
    cat(" ", crayon::silver("*"), " Loadings: ", crayon::green(loadings_size), "\n")
    cat(" ", crayon::silver("*"), " Total: ",    crayon::green(total_size),    "\n")

    # Sparsity
    n_nonzero <- sum(object@mask)
    sparsity <- round(100 * n_nonzero / prod(dims[1:3]), 2)
    cat("\n", crayon::yellow("Sparsity:"), "\n")
    cat(" ", crayon::silver("*"), " Non-zero voxels: ", crayon::green(n_nonzero), "\n")
    cat(" ", crayon::silver("*"), " Coverage: ",       crayon::green(sparsity), "%\n")

    # Space Info
    sp <- space(object)
    spacing_str <- paste(round(spacing(sp), 2), collapse=" x ")
    origin_str  <- paste(round(origin(sp), 2), collapse=" x ")
    cat("\n", crayon::yellow("Space Information:"), "\n")
    cat(" ", crayon::silver("*"), " Spacing: ", crayon::green(spacing_str), "\n")
    cat(" ", crayon::silver("*"), " Origin:  ", crayon::green(origin_str),  "\n")

    # Footer
    cat("\n", crayon::bold("Data Access:"), "\n")
    cat("\n", crayon::yellow("Reconstructed Space Access:"), "\n")
    cat(" ", crayon::silver("*"), " Extract volume: ",
        crayon::blue("object[[i]]"),
        crayon::silver("  # 3D volume at timepoint i\n"))
    cat(" ", crayon::silver("*"), " Get value: ",
        crayon::blue("object[i]"),
        crayon::silver("  # Value at linear index i\n"))
    cat(" ", crayon::silver("*"), " Subset: ",
        crayon::blue("object[mask]"),
        crayon::silver("  # Values at mask positions\n"))

    cat("\n", crayon::yellow("Latent Space Access:"), "\n")
    cat(" ", crayon::silver("*"), " Basis vectors: ",
        crayon::blue("basis(object)"),
        crayon::silver("  # nxk temporal basis\n"))
    cat(" ", crayon::silver("*"), " Loadings: ",
        crayon::blue("loadings(object)"),
        crayon::silver("  # pxk spatial loadings\n"))
    cat(" ", crayon::silver("*"), " Components: ",
        crayon::blue("components(object)"),
        crayon::silver("  # List of k component volumes\n"))

    cat("\n", crayon::yellow("Conversions:"), "\n")
    cat(" ", crayon::silver("*"), " as.array(object): ",
        crayon::silver("4D reconstruction\n"))
    cat(" ", crayon::silver("*"), " as.matrix(object): ",
        crayon::silver("nxp matrix of reconstructed values\n"))

    cat("\n", crayon::silver("Note: All access methods reconstruct data (X = B x L^T + offset)"),
        "\n", crayon::silver("unless you're explicitly accessing latent space."), "\n\n")
  }
)

#' Validity Check for LatentNeuroVec Objects
#' 
#' @param object A LatentNeuroVec object
#' @return TRUE if the object is valid, otherwise a character vector of error messages.
#' @keywords internal
#' @noRd
.validate_LatentNeuroVec <- function(object) {
    
    errors <- character()
    
    # Check types
    if (!inherits(object@basis, "Matrix")) {
        errors <- c(errors, "Slot @basis must be a Matrix object.")
    }
    if (!inherits(object@loadings, "Matrix")) {
        errors <- c(errors, "Slot @loadings must be a Matrix object.")
    }
    if (!is.numeric(object@offset)) {
        errors <- c(errors, "Slot @offset must be numeric.")
    }
    if (!inherits(object@mask, "LogicalNeuroVol")) {
        errors <- c(errors, "Slot @mask must be a LogicalNeuroVol object.")
    }
    if (!inherits(object@map, "IndexLookupVol")) {
        errors <- c(errors, "Slot @map must be an IndexLookupVol object.")
    }
    if (!is.character(object@label) || length(object@label) != 1) {
        errors <- c(errors, "Slot @label must be a single character string.")
    }
    if (!inherits(object@space, "NeuroSpace")) {
         errors <- c(errors, "Slot @space must be a NeuroSpace object.")
         # Stop further checks if space is invalid, as dims depend on it
         return(errors)
    }
    
    # If types are okay, check dimensions and consistency
    if (length(errors) == 0) {
        s_dims <- dim(object@space)
        if (length(s_dims) != 4) {
            errors <- c(errors, "Slot @space must have 4 dimensions.")
        } else {
            if (ncol(object@basis) != ncol(object@loadings)) {
                errors <- c(errors, paste0("Component mismatch: ncol(@basis) = ", ncol(object@basis), 
                                           " != ncol(@loadings) = ", ncol(object@loadings)))
            }
            
            if (nrow(object@basis) != s_dims[4]) {
                 errors <- c(errors, paste0("Time mismatch: nrow(@basis) = ", nrow(object@basis), 
                                            " != dim(@space)[4] = ", s_dims[4]))
            }
            
            dim_check_result <- validate_same_dims(
                object@mask, 
                object@space, 
                dims_to_compare = 1:3, 
                msg = "[.validate_LatentNeuroVec] Mask/Space dim mismatch:"
            )
            if (!is.null(dim_check_result)) {
                errors <- c(errors, dim_check_result)
            }

           
            nVox_mask <- sum(object@mask)
            if (nrow(object@loadings) != nVox_mask) {
                errors <- c(errors, paste0("Loadings rows (", nrow(object@loadings), 
                                           ") mismatch non-zero count in mask (", nVox_mask, ")"))
            }
           
            if (length(object@offset) > 0 && length(object@offset) != nrow(object@loadings)) {
                errors <- c(errors, paste0("Offset length (", length(object@offset), 
                                           ") mismatch number of rows in loadings (", nrow(object@loadings), ")"))
            }
           
            if (length(object@map@indices) != nVox_mask) {
                 errors <- c(errors, paste0("Map indices length (", length(object@map@indices), 
                                            ") mismatch non-zero count in mask (", nVox_mask, ")"))
            }
        }
    }

    if (length(errors) == 0) TRUE else errors
}


#' @keywords internal
setValidity("LatentNeuroVec", .validate_LatentNeuroVec)

# --- Accessor Methods for LatentNeuroVec --- 

#' @export
#' @rdname basis-methods
setMethod("basis", "LatentNeuroVec", function(x) x@basis)

#' @export
#' @rdname loadings-methods
setMethod("loadings", "LatentNeuroVec", function(x) x@loadings)

#' @export
#' @rdname offset-methods
setMethod("offset", "LatentNeuroVec", function(x) x@offset)

#' @export
#' @rdname mask-methods
setMethod("mask", "LatentNeuroVec", function(x) x@mask)

#' @export
#' @rdname map-methods
setMethod("map", "LatentNeuroVec", function(x) x@map)



#' Create HDF5 Dataset Property List
#' 
#' Internal helper to create and configure an HDF5 dataset creation property list (H5P_DATASET_CREATE).
#' Handles setting fill value, chunking, and compression (deflate/gzip) or nbit filter.
#'
#' @param dims The dimensions of the dataset for which the plist is being created.
#' @param dtype An H5T object representing the dataset data type.
#' @param chunk_dims A numeric vector specifying chunk dimensions. If NULL, no chunking is set. 
#'                   Must match the length of `dims`. Values are capped by `dims`.
#' @param compression Integer [0..9] specifying gzip compression level. 0 means no compression. Default 6. Ignored if nbit=TRUE.
#' @param nbit Logical, whether to use n-bit filter. Default FALSE. Takes precedence over compression.
#' @return A *copy* of the configured H5P_DATASET_CREATE object. The original is closed via on.exit.
#' @keywords internal
#' @noRd
.create_dset_plist <- function(dims, dtype, chunk_dims = NULL, compression = 6, nbit = FALSE) {
    plist <- hdf5r::H5P_DATASET_CREATE$new()
    # IMPORTANT: Ensure plist is closed when this function exits, even on error
    # Do NOT add the plist itself to the calling function's on.exit stack if returning a copy
    on.exit(plist$close(), add = TRUE) 
    
    # Set fill value based on dtype
    fill_val <- tryCatch({
        # First check that dtype is an H5T object with methods
        if (!inherits(dtype, "H5T") || !is.function(dtype$get_class)) {
            # For non-H5T objects, use a default value and issue a message
            message("Using default fill value 0 for non-H5T dtype")
            return(0)
        }
        
        dtype_class_char <- dtype$get_class()$to_text()
        # Expanded handling of HDF5 types with additional type classes 
        switch(dtype_class_char,
               # Float types (32-bit)
               "H5T_FLOAT"      = 0.0,
               "H5T_NATIVE_FLOAT" = 0.0,
               "H5T_IEEE_F32LE" = 0.0,
               "H5T_IEEE_F32BE" = 0.0,
               
               # Double types (64-bit)
               "H5T_DOUBLE"     = 0.0,
               "H5T_NATIVE_DOUBLE" = 0.0,
               "H5T_IEEE_F64LE" = 0.0,
               "H5T_IEEE_F64BE" = 0.0,
               
               # Integer types
               "H5T_INTEGER"    = 0L,
               "H5T_NATIVE_INT" = 0L,
               "H5T_NATIVE_INT8" = 0L,
               "H5T_NATIVE_INT16" = 0L,
               "H5T_NATIVE_INT32" = 0L,
               "H5T_NATIVE_INT64" = 0L,
               "H5T_STD_I8LE"   = 0L,
               "H5T_STD_I16LE"  = 0L,
               "H5T_STD_I32LE"  = 0L,
               "H5T_STD_I64LE"  = 0L,
               
               # Unsigned integer types
               "H5T_NATIVE_UINT" = 0L,
               "H5T_NATIVE_UINT8" = 0L,
               "H5T_NATIVE_UINT16" = 0L,
               "H5T_NATIVE_UINT32" = 0L,
               "H5T_NATIVE_UINT64" = 0L,
               "H5T_STD_U8LE"   = 0L,
               "H5T_STD_U16LE"  = 0L,
               "H5T_STD_U32LE"  = 0L,
               "H5T_STD_U64LE"  = 0L,
               
               # Character/string types
               "H5T_STRING"     = "",
               "H5T_C_S1"       = "",
               
               # Catch-all for other types
               {
                 # More informative message with dtype info
                 # Safely check for the existence of get_size method
                 dtype_size <- if (is.function(dtype$get_size)) dtype$get_size() else "unknown"
                 message("Using default fill value 0 for HDF5 dtype: ", dtype_class_char, 
                         " (", dtype_size, " bytes)")
                 0L  # Default to integer zero for unknown types
               }
        )
        }, error = function(e) {
            # Improved error handling
            message(paste0("Could not determine HDF5 fill value, using default 0. Error: ", e$message))
            0  # Default to numeric zero on error
        })
    
    tryCatch(plist$set_fill_value(dtype = dtype, value = fill_val), 
             error = function(e) warning("Could not set fill value for HDF5 dataset: ", e$message))
    
    # Set chunking if provided and valid
    if (!is.null(chunk_dims)) {
          # Get rank, handling NULL dims (scalar) and vectors (rank 1)
          data_rank <- if (is.null(dims)) 0 else length(dims)
          chunk_rank <- length(chunk_dims)

          if (chunk_rank == data_rank) {
              # Ensure chunk dims don't exceed dataset dims and are >= 1
              dims_vec <- dims # Use original dims for pmin
              safe_chunks <- pmax(1L, as.integer(pmin(chunk_dims, dims_vec)))
              plist$set_chunk(safe_chunks)
          } else {
              warning("Chunk dimensions length (", chunk_rank, ") does not match dataset rank (", data_rank, "). Ignoring chunking.")
          }
     }

    # Set compression: nbit takes precedence
    if (nbit) {
        # hdf5r doesn't have a direct set_nbit(). User needs external plugin.
        # We warn but don't explicitly set a filter here.
        warning("N-bit filter requested (nbit=TRUE). This requires a registered HDF5 filter plugin. Ensure it is available in your HDF5 library.")
        # Placeholder: If hdf5r gains support or a known filter ID exists:
        # tryCatch(plist$set_filter(filter_id = H5FILTER_NBIT, dcpl = ...), error = ... ) 
        if (compression > 0) {
            warning("N-bit filter requested; gzip compression level (", compression, ") will be ignored.")
        }
    } else if (compression > 0 && compression <= 9) {
        # Only set deflate if nbit is FALSE and compression > 0
        plist$set_deflate(as.integer(compression))
    }
    
    # Return a copy for the caller to use (caller is responsible for closing the copy)
    return(plist$copy())
}

#' Write Dense HDF5 Dataset
#'
#' Internal helper function to write a dense R array/vector to an HDF5 dataset.
#'
#' @param h5_group An open H5Group or H5File object where the dataset will be created.
#' @param name The name for the new dataset.
#' @param data The R data (vector, matrix, array) to write.
#' @param dtype An H5T object for the dataset datatype.
#' @param chunk_dims Optional chunk dimensions (numeric vector). See `.create_dset_plist`.
#' @param compression Integer [0..9] for gzip level. Default 6. Ignored if nbit=TRUE.
#' @param nbit Logical, whether to use n-bit filter. Default FALSE.
#' @return Invisible NULL. Called for side effects.
#' @keywords internal
#' @noRd
.write_h5_dataset <- function(h5_group, name, data, dtype, chunk_dims = NULL, compression = 6, nbit = FALSE) {
    dims <- dim(data)
    if (is.null(dims)) dims <- length(data) # Handle 1D vectors
    
    if (any(dims == 0) && !is.null(dims)) { # Handle zero-extent dimensions (e.g., empty matrix) 
        message(paste0("  Skipping write for dataset '", h5_group$get_obj_name(), "/", name, "' due to zero-length dimension."))
        # Create empty dataset for spec compliance
        space <- hdf5r::H5S$new(dims = dims, maxdims = dims)
        on.exit(space$close(), add = TRUE)
        # No chunk/compress/nbit for empty datasets
        plist <- .create_dset_plist(dims, dtype, chunk_dims = NULL, compression = 0, nbit = FALSE) 
        on.exit(plist$close(), add = TRUE) # Close the *copied* plist
        dset <- h5_group$create_dataset(name = name, space = space, dtype = dtype, dataset_create_pl = plist)
        dset$close()
        return(invisible(NULL))
    }    

    space <- NULL; plist <- NULL; dset <- NULL # Initialize for on.exit
    on.exit({
        if (!is.null(dset) && dset$is_valid) dset$close()
        if (!is.null(plist) && plist$is_valid) plist$close() # Close the *copied* plist
        if (!is.null(space) && space$is_valid) space$close()
    }, add = TRUE)

    space <- hdf5r::H5S$new(dims = dims, maxdims = dims)
    plist <- .create_dset_plist(dims, dtype, chunk_dims, compression, nbit) # Pass nbit
    dset <- h5_group$create_dataset(name = name, space = space, dtype = dtype, dataset_create_pl = plist)
    
    # Write data using appropriate subsetting based on rank
    rank <- length(dims)
    if (rank == 0 || rank == 1) { 
        dset[] <- data 
    } else if (rank == 2) {
        dset[,] <- data
    } else if (rank == 3) {
        dset[,,] <- data
    } else if (rank == 4) {
        dset[,,,] <- data
    } else {
        stop(".write_h5_dataset currently only supports up to 4 dimensions.")
    }
    message(paste0("  Dataset '", h5_group$get_obj_name(), "/", name, "' written."))
    invisible(NULL)
}

#' Convert LatentNeuroVec to Array
#'
#' @description
#' Converts a \code{LatentNeuroVec} to a standard R 4D array.
#' 
#' @param x A \code{LatentNeuroVec} object to convert.
#' @param ... Not used.
#'
#' @return A 4D array representing the full reconstructed data.
#'
#' @details
#' This method reconstructs the full 4D array from the latent representation, which can be 
#' memory-intensive for large datasets. It calculates:
#' \deqn{array[,,, t] = Map(basis[t,] \%*\% t(loadings)) + offset}
#' for each time point, where values outside the mask are set to zero.
#'
#' @export
setMethod(
  f = "as.array",
  signature = signature(x = "LatentNeuroVec"),
  definition = function(x, ...) {
    # Get dimensions from the space
    space_dims <- dim(x)
    
    # Initialize the output array
    result_array <- array(0, dim = space_dims)
    
    # Mask information
    mask_array <- as.logical(as.array(x@mask))
    mask_indices <- which(mask_array)

    for (t in seq_len(space_dims[4])) {
      # Get the basis vector for this time point
      basis_t <- x@basis[t, , drop = FALSE]
      
      # Calculate time point using matrix multiplication: basis_t %*% t(loadings)
      # tcrossprod is more efficient for this operation
      # Ensure both are Matrix objects for tcrossprod
      b1 <- basis_t
      b2 <- x@loadings
      # Matrix objects already - no need to convert
      values_in_mask <- as.vector(b1 %*% t(b2)) + x@offset
      
      # Map these values back to 3D space using the mask 
      # Initialize a 3D array for this time point
      vol_3d <- array(0, dim = space_dims[1:3])
      
      # Place values at the correct indices
      vol_3d[mask_indices] <- values_in_mask
      
      # Assign to the corresponding time slice in the result
      result_array[,,,t] <- vol_3d
    }
    
    return(result_array)
  }
)
