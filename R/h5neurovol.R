#' @include all_class.R
#' @importFrom assertthat assert_that
#' @importFrom neuroim2 NeuroSpace space origin trans
#' @importFrom lifecycle deprecate_warn
#' @importFrom hdf5r h5attr H5File h5types
#' @importFrom crayon bold blue silver yellow green italic red
NULL


#' H5NeuroVol Constructor
#'
#' @description
#' Constructs an \code{\link{H5NeuroVol}} object representing a 3D brain volume
#' stored in an HDF5 file. The HDF5 file is opened in read-only mode.
#' 
#' @details
#' This constructor is typically used for reading existing HDF5 files that conform
#' to the H5NeuroVol specification.
#'
#' @section Lifecycle Management:
#' When an \code{H5NeuroVol} object is created by providing a \code{file_name},
#' it opens the specified HDF5 file and maintains an open handle to it.
#' **It is the user's responsibility to explicitly close this handle** when the
#' object is no longer needed to release system resources. This can be done by calling
#' \code{close(your_h5neurovol_object)}.
#'
#' Failure to close the handle may lead to issues such as reaching file handle
#' limits or problems with subsequent access to the file.
#'
#' @param file_name A \code{character} string giving the path to an existing 3D HDF5 neuroimaging file.
#' @return A new \code{\link{H5NeuroVol-class}} instance with an open HDF5 file handle.
#'
#' @seealso \code{\link{close.H5NeuroVol}} for closing the file handle, \code{\link[neuroim2]{NeuroVol-class}}
#'
#' @examples
#' \dontrun{
#' # Assuming "my_volume.h5" is a valid H5NeuroVol HDF5 file
#' h5vol <- H5NeuroVol("my_volume.h5")
#' # ... perform operations with h5vol ...
#' print(dim(h5vol))
#' # Important: Close the handle when done
#' close(h5vol)
#' }
#' @export
H5NeuroVol <- function(file_name) {
  assert_that(is.character(file_name))
  assert_that(file.exists(file_name))

  h5obj <- hdf5r::H5File$new(file_name)

  # Check the "rtype" attribute
  rtype <- try(hdf5r::h5attr(h5obj, which="rtype"), silent=TRUE)
  if (!is.character(rtype) || rtype != "DenseNeuroVol") {
    stop("Invalid HDF5 file for H5NeuroVol: ", file_name)
  }

  # Check dimension count
  if (length(h5obj[["space/dim"]][]) != 3) {
    stop(
      "Cannot create H5NeuroVol: file must have 3 dimensions; found: ",
      paste(h5obj[["space/dim"]][], collapse=" ")
    )
  }

  # Build NeuroSpace
  sp <- NeuroSpace(
    dim    = h5obj[["space/dim"]][],
    origin = h5obj[["space/origin"]][],
    trans  = h5obj[["space/trans"]][,]
  )

  new("H5NeuroVol", space=sp, h5obj=h5obj)
}

#' Convert a NeuroVol to HDF5 Format (as_h5 Method)
#'
#' @description
#' Saves a NeuroVol to an HDF5 file with minimal necessary metadata
#' to reconstruct an H5NeuroVol.
#'
#' @param object A NeuroVol object (3D)
#' @param file Path to the output file (if NULL, uses tempfile)
#' @param data_type Character: "FLOAT", "DOUBLE", "INT", etc.
#' @param chunk_dim Numeric vector specifying chunk sizes
#' @param compression Integer [1..9], default 6
#'
#' @return A new \code{H5NeuroVol} referencing the written file.
#'   The returned object contains an open read-mode HDF5 handle.
#'   **Important:** The user is responsible for closing this handle using
#'   \code{close()} on the returned object when finished.
#'
#' @export
setMethod(
  f = "as_h5",
  signature = signature(object = "NeuroVol"),
  definition = function(object, file = NULL, data_type = "FLOAT",
                         chunk_dim = NULL, compression = 6) { 
    
    # --- Determine output file path --- 
    out_file <- file
    if (is.null(out_file)) {
        out_file <- tempfile(fileext = ".h5")
        message("Output file not specified, using temp file: ", out_file)
    }
    
    # --- Map data_type string to HDF5 type object --- 
    h5dtype_obj <- switch(toupper(data_type),
                           "FLOAT"   = hdf5r::h5types$H5T_NATIVE_FLOAT,
                           "DOUBLE"  = hdf5r::h5types$H5T_NATIVE_DOUBLE,
                           "INT"     = hdf5r::h5types$H5T_NATIVE_INT32,
                           "INTEGER" = hdf5r::h5types$H5T_NATIVE_INT32,
                           "SHORT"   = hdf5r::h5types$H5T_NATIVE_INT16,
                           "CHAR"    = hdf5r::h5types$H5T_NATIVE_CHAR,
                           "UINT8"   = hdf5r::h5types$H5T_NATIVE_UCHAR,
                           stop("Unsupported data_type: ", data_type)
                          )
    
    # --- Write Phase: Open, Write, Explicitly Close --- 
    fh_write <- NULL 
    h5_write_obj <- NULL
    write_success <- FALSE
    tryCatch({
        fh_write <- open_h5(out_file, mode = "w")
        h5_write_obj <- fh_write$h5
        
        # Write Root Attribute
        hdf5r::h5attr(h5_write_obj, "rtype") <- "DenseNeuroVol"
        
        # Write Space Group 
        sp <- space(object)
        sp_dims <- dim(sp)
        if (length(sp_dims) != 3) stop("Input NeuroVol must be 3-dimensional.")
        
        # --- Debug: Print trans matrix before writing ---
        current_trans <- trans(sp)
        #message("DEBUG: trans matrix BEFORE writing:")
        print(current_trans)
        
        h5_write(h5_write_obj, "/space/dim", as.integer(sp_dims), overwrite = TRUE)
        h5_write(h5_write_obj, "/space/origin", as.double(origin(sp)), overwrite = TRUE)
        h5_write(h5_write_obj, "/space/spacing", as.double(spacing(sp)), overwrite = TRUE)
        h5_write(h5_write_obj, "/space/trans", current_trans, overwrite = TRUE) 
        
        # Write Data Group
        data_arr <- as.array(object) 
        final_chunk_dim <- chunk_dim
        if (is.null(final_chunk_dim) && compression > 0) {
             final_chunk_dim <- pmin(sp_dims, c(32L, 32L, 32L))
        }
        h5_write(h5_write_obj, "/data/elements", data_arr,
                 dtype = h5dtype_obj,
                 chunk_dims = final_chunk_dim, 
                 compression = compression,
                 overwrite = TRUE)
        
        write_success <- TRUE # Mark success if we reached here
        message("Successfully wrote NeuroVol to: ", out_file)
        
    }, error = function(e) {
        # Error occurred during writing
        stop("Failed during HDF5 write phase for ", out_file, ": ", e$message)
    }, finally = {
        # Always attempt to close the write handle
        if (!is.null(fh_write) && fh_write$owns) {
             safe_h5_close(h5_write_obj)
        }
    })
    
    # If writing didn't complete successfully, stop before reopening
    if (!write_success) {
        stop("HDF5 write failed for ", out_file, ", cannot create H5NeuroVol object.")
    }
    
    # --- Reopen in Read Mode and Return H5NeuroVol --- 
    h5_read_obj <- NULL
    tryCatch({
        # Write handle is now closed. Reopen the file in read mode.
        h5_read_obj <- hdf5r::H5File$new(out_file, mode = "r")
        
        # Read space info for the new object
        trans_data <- h5_read_obj[["/space/trans"]]$read()
       
        spacing <- diag(trans_data)[1:3]
        sp_read <- NeuroSpace(
            dim    = h5_read_obj[["/space/dim"]]$read(),
            spacing= h5_read_obj[["/space/spacing"]]$read(),   
            origin = h5_read_obj[["/space/origin"]]$read(), 
            trans  = trans_data
        )
      
        
        # Return the H5NeuroVol with the NEW, OPEN, read-mode handle
        new("H5NeuroVol", space=sp_read, h5obj=h5_read_obj)
        
    }, error = function(e) {
        # If reopening or reading space fails, ensure the read handle is closed if it opened
        if (!is.null(h5_read_obj) && h5_read_obj$is_valid) safe_h5_close(h5_read_obj)
        stop("Failed to reopen/read HDF5 file ", out_file, " to create H5NeuroVol: ", e$message)
    })
  }
)


#' @export
#' @rdname linear_access-methods
setMethod(
  f = "linear_access",
  signature = signature(x="H5NeuroVol", i="numeric"),
  definition = function(x, i) {

    # 1) Check range
    n_vox <- prod(dim(x))  # total number of voxels in 3D
    if (any(i < 1 | i > n_vox)) {
      stop("Some linear indices are out of range 1..", n_vox)
    }
    # If you also want an error if i==0 is found, the above condition will catch it.

    # 2) If i is empty => return empty
    if (length(i) == 0) {
      return(numeric(0))
    }

    # 3) Convert linear -> (x,y,z)
    coords <- arrayInd(i, dim(x))  # Nx3

    # 4) bounding box
    minx <- min(coords[,1]); maxx <- max(coords[,1])
    miny <- min(coords[,2]); maxy <- max(coords[,2])
    minz <- min(coords[,3]); maxz <- max(coords[,3])

    # 5) Read that bounding box from the dataset
    subvol <- with_h5_dataset(x@h5obj, "data/elements", function(ds) {
        if (is.null(ds)) stop("Could not open dataset '/data/elements' for linear_access")
        ds[minx:maxx, miny:maxy, minz:maxz, drop = FALSE]
    })
    # shape => (maxx - minx + 1) x (maxy - miny + 1) x (maxz - minz + 1)

    # 6) Offset coords to index subvol
    off_coords <- cbind(coords[,1] - minx + 1,
                        coords[,2] - miny + 1,
                        coords[,3] - minz + 1)

    # 7) Gather values
    n <- nrow(coords)
    out_vals <- numeric(n)
    for (k in seq_len(n)) {
      out_vals[k] <- subvol[ off_coords[k,1],
                             off_coords[k,2],
                             off_coords[k,3] ]
    }
    out_vals
  }
)

#' @export
setMethod(
  f = "linear_access",
  signature = signature(x="H5NeuroVol", i="integer"),
  definition = function(x, i) {
    callGeneric(x, as.numeric(i))  # passes off to the numeric method above
  }
)



#' 3D bracket subsetting for H5NeuroVol (handles partial arguments)
#'
#' @description
#' Allows \code{h5vol[i, j, k]} where each of \code{i,j,k} may be missing.
#' Missing arguments default to the entire range in that dimension.
#' Zero-length arguments immediately yield an empty array of the correct shape.
#'
#' @param x An \code{H5NeuroVol} instance
#' @param i,j,k Numeric (or integer) index vectors for each dimension. If missing,
#'   we take the full range in that dimension.
#' @param drop Logical: whether to drop dimensions of size 1. Default \code{TRUE}.
#' @param ... Unused
#'
#' @return A numeric \code{array} of shape \code{c(length(i), length(j), length(k))},
#'   or fewer dims if \code{drop=TRUE}.
#'
#' @export
setMethod(
  f = "[",
  signature = signature(x="H5NeuroVol"),
  definition = function(x, i, j, k, ..., drop=TRUE) {

    # 1) Determine dimension of the underlying volume
    dimx <- dim(x)  # c(X, Y, Z)
    if (length(dimx) != 3) {
      stop("H5NeuroVol is not 3D? Found dim=", paste(dimx, collapse="x"))
    }

    # 2) If i, j, k are missing, default them
    if (missing(i)) {
      i <- seq_len(dimx[1])
    }
    if (missing(j)) {
      j <- seq_len(dimx[2])
    }
    if (missing(k)) {
      k <- seq_len(dimx[3])
    }

    # Convert to numeric in case user gave integer
    i <- as.numeric(i)
    j <- as.numeric(j)
    k <- as.numeric(k)

    # 3) If any index has length=0 => return an empty array right away
    if (length(i)==0 || length(j)==0 || length(k)==0) {
      outdim <- c(length(i), length(j), length(k))
      empty_arr <- array(numeric(0), dim=outdim)
      if (drop) empty_arr <- drop(empty_arr)
      return(empty_arr)
    }

    # 4) Check out-of-range
    if (min(i)<1 || max(i)>dimx[1]) {
      stop("Subscript 'i' out of range for dimension 1")
    }
    if (min(j)<1 || max(j)>dimx[2]) {
      stop("Subscript 'j' out of range for dimension 2")
    }
    if (min(k)<1 || max(k)>dimx[3]) {
      stop("Subscript 'k' out of range for dimension 3")
    }

    # 5) Determine bounding box
    minI <- floor(min(i)); maxI <- ceiling(max(i))
    minJ <- floor(min(j)); maxJ <- ceiling(max(j))
    minK <- floor(min(k)); maxK <- ceiling(max(k))

    # 6) Read the bounding box from the dataset
    subvol <- with_h5_dataset(x@h5obj, "data/elements", function(ds) {
        if (is.null(ds)) stop("Could not open dataset '/data/elements'")
        ds[minI:maxI, minJ:maxJ, minK:maxK, drop = FALSE]
    })
    # shape => c((maxI-minI+1), (maxJ-minJ+1), (maxK-minK+1))

    # 7) We then re-map i,j,k into local sub-box coords
    i_off <- i - minI + 1
    j_off <- j - minJ + 1
    k_off <- k - minK + 1

    subdimI <- maxI - minI + 1
    subdimJ <- maxJ - minJ + 1

    # We'll build the output array
    out_dim <- c(length(i), length(j), length(k))
    out_vals <- numeric(prod(out_dim))

    # Flatten subvol
    subvol_vec <- as.vector(subvol)

    # Build a systematic index
    N  <- length(i) * length(j) * length(k)
    ix_i <- rep(seq_along(i), times = length(j)*length(k))
    ix_j <- rep(rep(seq_along(j), each=length(i)), times=length(k))
    ix_k <- rep(seq_along(k), each = length(i)*length(j))

    loc_i <- i_off[ix_i]
    loc_j <- j_off[ix_j]
    loc_k <- k_off[ix_k]

    # local 3D => linear index in subvol
    sub_lin_idx <- loc_i +
      (loc_j-1)* subdimI +
      (loc_k-1)* subdimI * subdimJ

    out_vals <- subvol_vec[sub_lin_idx]
    arr_out  <- array(out_vals, dim=out_dim)

    # 8) drop dims if requested
    if (drop) {
      arr_out <- drop(arr_out)
    }
    arr_out
  }
)

#' Convert DenseNeuroVol to H5NeuroVol
#'
#' @description
#' Converts a \code{DenseNeuroVol} to \code{H5NeuroVol} by writing it to disk in HDF5 format.
#'
#' @param from A \code{DenseNeuroVol} object.
#' @return An \code{H5NeuroVol} referencing the resulting HDF5 file.
#'
#' @keywords internal
#' @name DenseNeuroVol,H5NeuroVol
setAs(
  from = "DenseNeuroVol",
  to   = "H5NeuroVol",
  def  = function(from) {
    to_nih5_vol(from, file_name=NULL, data_type="FLOAT")
  }
)

#' Show Method for H5NeuroVol
#' 
#' Displays a summary of the H5NeuroVol object, including dimensions,
#' spacing, origin, and HDF5 file information, without reading voxel data.
#' 
#' @param object The H5NeuroVol object to display.
#' @importFrom crayon bold blue silver yellow green italic
#' @export
setMethod(
  f = "show",
  signature = "H5NeuroVol",
  definition = function(object) {
    cat("\n", crayon::bold(crayon::blue("H5NeuroVol Object")), "\n")
    cat(crayon::silver("═══════════════════\n"))
    
    # Display spatial info from the space slot
    sp <- object@space
    dims <- dim(sp)
    spacing_str <- paste(round(neuroim2::spacing(sp), 2), collapse = " × ")
    origin_str <- paste(round(neuroim2::origin(sp), 2), collapse = " × ")
    
    cat(crayon::yellow("Dimensions:"), crayon::green(paste(dims, collapse = " × ")), "\n")
    cat(crayon::yellow("Spacing:"), crayon::green(spacing_str), "\n")
    cat(crayon::yellow("Origin:"), crayon::green(origin_str), "\n")
    
    # Display HDF5 file info
    h5f <- object@h5obj
    file_status <- "<Invalid Handle>"
    file_path <- "<unknown>"
    
    if (!is.null(h5f) && inherits(h5f, "H5File")) {
      is_valid <- tryCatch(h5f$is_valid, error = function(e) FALSE)
      if (is_valid) {
        file_path <- tryCatch(h5f$get_filename(), error = function(e) "<error getting path>")
        file_status <- paste(crayon::green("VALID handle"), "for file:", crayon::italic(file_path))
      } else {
        # Try to get filename even if closed
        file_path <- tryCatch(h5f$get_filename(), error = function(e) "<unknown path>")
        file_status <- paste(crayon::red("CLOSED handle"), "for file:", crayon::italic(file_path))
      }
    } else {
       file_status <- crayon::red("INVALID H5File object slot")
    }
    
    cat(crayon::yellow("HDF5 Source:"), file_status, "\n")
    cat(crayon::silver("═══════════════════\n"))
    cat("Access data using standard array indexing (e.g., object[1:10, 1:10, 1])\n")
    invisible(NULL)
  }
)

#' Close the HDF5 file associated with an H5NeuroVol
#'
#' This method manually closes the HDF5 file handle stored within the
#' H5NeuroVol object. It uses the \code{safe_h5_close} helper.
#'
#' @param con An \code{H5NeuroVol} object.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns \code{NULL}.
#' @rdname close
#' @export
setMethod("close", "H5NeuroVol", function(con, ...) {
  if (!is.null(con@h5obj)) {
    safe_h5_close(con@h5obj)
    # Nulling out the reference is problematic as the slot expects H5File object.
    # The hdf5r object itself will become invalid after closing.
    # con@h5obj <- NULL # This line causes S4 validation error
  }
  invisible(NULL)
})
