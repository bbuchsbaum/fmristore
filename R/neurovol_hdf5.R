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
#' @param x A \code{\link{H5NeuroVol}} object.
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
    # 1) Convert linear -> (x,y,z)
    coords <- arrayInd(i, dim(x))  # Nx3

    # 2) bounding box
    minx <- min(coords[,1]); maxx <- max(coords[,1])
    miny <- min(coords[,2]); maxy <- max(coords[,2])
    minz <- min(coords[,3]); maxz <- max(coords[,3])

    # 3) Read that bounding box
    # Adjust path if needed ("/data/elements" vs "/data"). We'll assume "/data/elements".
    dset <- x@h5obj[["data/elements"]]
    subvol <- dset[minx:maxx, miny:maxy, minz:maxz, drop=FALSE]

    # 4) Offset coords for subvol indexing
    off_coords <- cbind(coords[,1] - minx + 1,
                        coords[,2] - miny + 1,
                        coords[,3] - minz + 1)

    # 5) Gather values
    n <- nrow(coords)
    out_vals <- numeric(n)
    for (k in seq_len(n)) {
      cx <- off_coords[k,1]
      cy <- off_coords[k,2]
      cz <- off_coords[k,3]
      out_vals[k] <- subvol[cx, cy, cz]
    }
    out_vals
  }
)

#' @export
setMethod(
  f = "linear_access",
  signature = signature(x="H5NeuroVol", i="integer"),
  definition = function(x, i) {
    callGeneric(x, as.numeric(i))
  }
)

#' Subset a 3D H5NeuroVol via bounding box (vectorized indexing)
#'
#' @description
#' Handles calls like \code{x[i, j, k]} where \code{i,j,k} are numeric vectors (potentially
#' non-contiguous). Reads the minimal bounding box from HDF5 in one go, then picks out
#' the requested values in correct order (without triple nested loops in R).
#'
#' @param x A \code{\link{H5NeuroVol}} object (3D).
#' @param i Numeric indices for dimension 1.
#' @param j Numeric indices for dimension 2.
#' @param k Numeric indices for dimension 3.
#' @param drop Logical; if TRUE, drops dimensions of size 1.
#' @param ... Ignored.
#'
#' @return An \code{array} of size \code{length(i) × length(j) × length(k)}, or fewer
#'   dimensions if \code{drop=TRUE}.
#'
#' @export
setMethod(
  f = "[",
  signature = signature(x="H5NeuroVol", i="numeric", j="numeric", drop="ANY"),
  definition = function(x, i, j, k, ..., drop=TRUE) {
    # Convert i,j,k
    i <- as.numeric(i)
    j <- as.numeric(j)
    k <- as.numeric(k)

    # Basic bounds
    dimx <- dim(x)
    if (any(i<1 | i>dimx[1]) || any(j<1|j>dimx[2]) || any(k<1|k>dimx[3])) {
      stop("Subscript out of range for H5NeuroVol dims.")
    }

    # bounding box
    minI <- min(i); maxI <- max(i)
    minJ <- min(j); maxJ <- max(j)
    minK <- min(k); maxK <- max(k)

    # read bounding box
    dset <- x@h5obj[["data/elements"]]
    subvol <- dset[minI:maxI, minJ:maxJ, minK:maxK, drop=FALSE]

    # offsets
    i_off <- i - minI + 1
    j_off <- j - minJ + 1
    k_off <- k - minK + 1

    subdimI <- maxI - minI + 1
    subdimJ <- maxJ - minJ + 1

    # result shape
    nI <- length(i)
    nJ <- length(j)
    nK <- length(k)
    out_dim <- c(nI, nJ, nK)

    # flatten subvol
    subvol_vec <- as.vector(subvol)

    # build linear indices
    N <- nI * nJ * nK
    ix_i <- rep(seq_len(nI), times=nJ*nK)
    ix_j <- rep(rep(seq_len(nJ), each=nI), times=nK)
    ix_k <- rep(seq_len(nK), each=nI*nJ)

    loc_i <- i_off[ix_i]
    loc_j <- j_off[ix_j]
    loc_k <- k_off[ix_k]

    sub_lin_idx <- loc_i +
      (loc_j-1)* subdimI +
      (loc_k-1)* subdimI * subdimJ

    out_vals <- subvol_vec[sub_lin_idx]
    arr_out <- array(out_vals, dim=out_dim)
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
