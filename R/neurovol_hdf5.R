#' @include all_class.R
NULL

#' H5NeuroVol Constructor
#'
#' @description
#' Constructs an \code{\link{H5NeuroVol}} object representing a 3D brain volume
#' stored in HDF5 format.
#'
#' @param file_name A \code{character} string giving the path to a 3D HDF5 file.
#' @return A new \code{\link{H5NeuroVol-class}} instance.
#'
#' @seealso
#' \code{\link[neuroim2]{NeuroVol-class}} for the base 3D brain volume class.
#'
#' @importFrom assertthat assert_that
#' @importFrom neuroim2 NeuroSpace
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

#' A simple bounding-box approach for linear indexing of H5NeuroVol
#'
#' @description
#' Given a set of linear voxel indices in a 3D volume, this method maps them to
#' \code{(x,y,z)} coordinates, reads the minimal bounding box from HDF5 once,
#' and extracts the requested values in the correct order.
#'
#' @param x An \code{\link{H5NeuroVol}} object.
#' @param i A numeric vector of linear voxel indices (1-based).
#'
#' @return A numeric vector of length \code{length(i)} containing the voxel values
#'   in the same order as \code{i}.
#'
#' @export
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
    dset <- x@h5obj[["data/elements"]]
    subvol <- dset[minx:maxx, miny:maxy, minz:maxz, drop=FALSE]
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
    dset   <- x@h5obj[["data/elements"]]
    subvol <- dset[minI:maxI, minJ:maxJ, minK:maxK, drop=FALSE]
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

#' Convert a NeuroVol to HDF5 Format
#'
#' @description
#' Saves a \code{\link[neuroim2]{NeuroVol}} to an HDF5 file with minimal necessary
#' metadata to reconstruct an \code{H5NeuroVol}.
#'
#' @param vol A \code{NeuroVol} object (3D).
#' @param file_name \code{character} path to the output file (if NULL, uses tempfile).
#' @param data_type \code{character}: "FLOAT", "DOUBLE", "INT", etc.
#' @param chunk_dim \code{numeric} vector specifying chunk sizes.
#' @param nbit \code{logical}; if TRUE, applies N-bit compression.
#' @param compression \code{integer} [1..9], default 6.
#'
#' @return A new \code{\link{H5NeuroVol}} referencing the written file.
#'
#' @export
to_nih5_vol <- function(vol, file_name=NULL, data_type="FLOAT",
                        chunk_dim=NULL, nbit=FALSE, compression=6)
{
  if (is.null(file_name)) {
    file_name <- paste0(tempfile(), ".h5")
  }

  # 1) Open HDF5 file
  h5obj <- hdf5r::H5File$new(file_name, mode="w")

  # 2) Create dataset
  space_ds <- hdf5r::H5S$new(dims=dim(vol), maxdims=dim(vol))
  dtype_pl  <- hdf5r::H5P_DATASET_CREATE$new()

  # match data_type
  h5dtype <- switch(data_type,
                    "BINARY"  = hdf5r::h5types$H5T_NATIVE_HBOOL,
                    "SHORT"   = hdf5r::h5types$H5T_NATIVE_SHORT,
                    "INT"     = hdf5r::h5types$H5T_NATIVE_INT,
                    "INTEGER" = hdf5r::h5types$H5T_NATIVE_INT,
                    "FLOAT"   = hdf5r::h5types$H5T_NATIVE_FLOAT,
                    "DOUBLE"  = hdf5r::h5types$H5T_NATIVE_DOUBLE,
                    "LONG"    = hdf5r::h5types$H5T_NATIVE_LONG,
                    NULL
  )
  if (is.null(h5dtype)) {
    stop("Unsupported data_type: ", data_type)
  }

  if (is.null(chunk_dim)) {
    # default chunking: entire X,Y, but chunk Z by 1
    chunk_dim <- c(dim(vol)[1], dim(vol)[2], 1)
  }
  dtype_pl <- dtype_pl$set_chunk(chunk_dim)$set_fill_value(h5dtype, 0)
  if (compression>0) {
    dtype_pl <- dtype_pl$set_deflate(compression)
  }
  if (nbit) {
    dtype_pl <- dtype_pl$set_nbit()
  }

  # 3) Create /data
  dgroup <- h5obj$create_group("data")

  # 4) Write data
  dset <- dgroup$create_dataset(
    name             = "elements",
    space            = space_ds,
    dtype            = h5dtype,
    dataset_create_pl= dtype_pl,
    chunk_dims       = chunk_dim
  )
  dset[,,] <- vol@.Data

  # 5) Mark file attribute for H5NeuroVol
  hdf5r::h5attr(h5obj, "rtype") <- "DenseNeuroVol"

  # 6) Store space info
  space_grp <- h5obj$create_group("space")
  space_grp[["dim"]]     <- dim(vol)
  space_grp[["origin"]]  <- origin(space(vol))
  space_grp[["spacing"]] <- spacing(space(vol))
  space_grp[["trans"]]   <- trans(space(vol))
  space_grp[["axes"]]    <- c("x","y","z")  # minimal

  # 7) Return the new H5NeuroVol
  new("H5NeuroVol", space=space(vol), h5obj=h5obj)
}

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
