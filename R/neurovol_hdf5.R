#' @include all_class.R
NULL

#' H5NeuroVol Constructor
#'
#' @description
#' This function constructs an H5NeuroVol object, which represents a three-dimensional
#' brain image stored in the HDF5 format.
#'
#' @param file_name The name of the 3-dimensional image file in HDF5 format.
#'
#' @return An instance of the \code{\link{H5NeuroVol-class}} class.
#'
#' @seealso
#' \code{\link[neuroim2]{NeuroVol-class}} for the base 3D brain volume class.
#'
#' @importFrom assertthat assert_that
#' @importFrom neuroim2 NeuroSpaceQ
#' @export
H5NeuroVol <- function(file_name) {
  assert_that(is.character(file_name))
  assert_that(file.exists(file_name))

  h5obj <- hdf5r::H5File$new(file_name)

  rtype <- try(hdf5r::h5attr(h5obj, which="rtype"))
  if (! (rtype == "DenseNeuroVol")) {
    stop("invalid h5 file: ", file_name)
  }

  if (length(h5obj[["space/dim"]][]) != 3) {
    stop(paste("cannot H5NeuroVol: must have 3 dimensions: ", paste(h5obj[["space/dim"]][], collapse=" ")))
  }

  sp <- NeuroSpace(dim=h5obj[["space/dim"]][], origin=h5obj[["space/origin"]][],
                   trans=h5obj[["space/trans"]][,])

  new("H5NeuroVol", space=sp, h5obj=h5obj)
}


# # Methods for H5NeuroVol
# setMethod("[", signature(x="H5NeuroVol", i="numeric", j="numeric", k="numeric"),
#           def=function(x, i, j, k) {
#             assertthat::assert_that(max(i) <= dim(x)[1])
#             assertthat::assert_that(max(j) <= dim(x)[2])
#             assertthat::assert_that(max(k) <= dim(x)[3])
#             assertthat::assert_that(min(i) >= 1)
#             assertthat::assert_that(min(j) >= 1)
#             assertthat::assert_that(min(k) >= 1)
#             x@h5obj[["data"]][i, j, k, drop=FALSE]
#           })


#' A simple bounding-box approach for linear indexing of H5NeuroVol
#'
#' @param x A H5NeuroVol object
#' @param i Numeric vector of linear voxel indices
#'
#' @return A numeric vector containing values at those voxel indices
#' @export
setMethod(
  f = "linear_access",
  signature = signature(x = "H5NeuroVol", i = "numeric"),
  definition = function(x, i) {
    # 1) Convert linear -> (x,y,z) coords
    coords <- arrayInd(i, dim(x))    # coords is N x 3

    # 2) Find bounding box of requested coords
    minx <- min(coords[,1]); maxx <- max(coords[,1])
    miny <- min(coords[,2]); maxy <- max(coords[,2])
    minz <- min(coords[,3]); maxz <- max(coords[,3])

    # 3) Read that bounding box from HDF5
    #    Adjust path if needed ("/data" vs. "/data/elements", etc.).
    dset <- x@h5obj[["data/elements"]]
    subvol <- dset[minx:maxx, miny:maxy, minz:maxz, drop=FALSE]

    # 4) Convert coords so we can index into subvol
    #    subvol has dimensions c( maxx-minx+1, maxy-miny+1, maxz-minz+1 )
    #    So subtract offsets from coords
    offset_coords <- cbind(coords[,1] - minx + 1,
                           coords[,2] - miny + 1,
                           coords[,3] - minz + 1)

    # 5) Extract the requested voxels from subvol in the user-specified order
    #    We'll do this in a loop or apply; a loop is simplest but for large i,
    #    you might want a vectorized approach.
    n <- nrow(coords)
    out_vals <- numeric(n)
    for (k in seq_len(n)) {
      cx <- offset_coords[k,1]
      cy <- offset_coords[k,2]
      cz <- offset_coords[k,3]
      out_vals[k] <- subvol[cx, cy, cz]
    }

    # 6) Return in the same order as i was given
    return(out_vals)
  }
)


#' @export
setMethod(
  f = "linear_access",
  signature = signature(x = "H5NeuroVol", i = "integer"),
  definition = function(x, i) {
    callGeneric(x, as.numeric(i))
  }
)

#' Subset a 3D H5NeuroVol with arbitrary numeric i, j, k via bounding box (vectorized)
#'
#' @description
#' Handles calls like \code{x[i, j, k]} where \code{i}, \code{j}, and \code{k} are
#' numeric vectors (possibly non-contiguous). We read the minimal bounding box from
#' HDF5 in one go, then pick out the requested points in the correct order *without*
#' a triple nested loop.
#'
#' @param x A \code{H5NeuroVol} object (3D).
#' @param i Numeric vector of indices for the first dimension.
#' @param j Numeric vector of indices for the second dimension.
#' @param k Numeric vector for the third dimension.
#' @param drop Logical whether to drop dimensions of size 1.
#' @param ... Additional args (ignored).
#'
#' @return An R array (or vector if dimensions are dropped) of size
#'   \code{length(i) x length(j) x length(k)} (then dropped if \code{drop=TRUE}).
#'
#' @export
setMethod(
  f = "[",
  signature = signature(x="H5NeuroVol", i="numeric", j="numeric", drop="ANY"),
  definition = function(x, i, j, k, ..., drop=TRUE) {

    # 1) Convert i, j, k to numeric if not already
    i <- as.numeric(i)
    j <- as.numeric(j)
    k <- as.numeric(k)

    # 2) Basic bounds checking
    dimx <- dim(x)  # c(nx, ny, nz)
    if (any(i < 1 | i > dimx[1]) ||
        any(j < 1 | j > dimx[2]) ||
        any(k < 1 | k > dimx[3])) {
      stop("Subscript out of range for one of the dimensions in H5NeuroVol.")
    }

    # 3) Identify bounding box
    minI <- min(i); maxI <- max(i)
    minJ <- min(j); maxJ <- max(j)
    minK <- min(k); maxK <- max(k)

    # 4) Read bounding box from HDF5 in a single call
    dset <- x@h5obj[["data/elements"]]
    subvol <- dset[minI:maxI, minJ:maxJ, minK:maxK, drop=FALSE]
    # subvol dims: (maxI-minI+1) x (maxJ-minJ+1) x (maxK-minK+1)

    # 5) Offsets to map subvol( minI,minJ,minK ) => local index (1,1,1)
    i_off <- i - minI + 1
    j_off <- j - minJ + 1
    k_off <- k - minK + 1

    subdimI <- maxI - minI + 1
    subdimJ <- maxJ - minJ + 1
    # subdimK <- maxK - minK + 1  # Not strictly needed for the formula.

    # 6) We want an output array of size length(i) x length(j) x length(k)
    nI <- length(i)
    nJ <- length(j)
    nK <- length(k)
    out_dim <- c(nI, nJ, nK)

    # Flatten subvol into a 1D vector in R's column-major order
    # so subvol_vec[ a + (b-1)*subdimI + (c-1)*subdimI*subdimJ ]
    # corresponds to subvol[a,b,c].
    subvol_vec <- as.vector(subvol)

    # 7) Build a vector of linear indices for each (x_i, y_j, z_k)
    #    R is column-major, so first dimension changes fastest => i is fastest, then j, then k.
    #    We'll do a single pass to construct them all.

    # a) Expand i, j, k in the correct order:
    # i cycles fastest, j next, k slowest => typical R array layout.
    # E.g. if we had i=2, j=3, k=4, the total length is 2*3*4=24. We'll replicate indexes accordingly.
    # xi = 1,2,1,2,1,2, ...
    # yj = 1,1,2,2,3,3, ...
    # zk changes last.
    n <- nI * nJ * nK
    # We'll define vectors of the same length n:
    xi <- rep(seq_len(nI), times = nJ*nK)
    yj <- rep(rep(seq_len(nJ), each = nI), times = nK)
    zk <- rep(seq_len(nK), each = nI*nJ)

    # b) For each triple (xi[p], yj[p], zk[p]) in out array space,
    #    map to subvol index: i_off, j_off, k_off
    #    Then convert to 1D offset in subvol_vec
    #    subvol's linear index =
    #    = ( i_off[xi[p]] ) +
    #      ( j_off[yj[p]] - 1 )* subdimI +
    #      ( k_off[zk[p]] - 1 )* subdimI*subdimJ
    #
    # We do -1 because R indexing starts at 1, so for dimension offsets we do (val - 1).

    # c) Actually build them:
    li_i <- i_off[xi]  # local i coords
    li_j <- j_off[yj]  # local j coords
    li_k <- k_off[zk]  # local k coords

    subvol_lin_index <-
      li_i +
      (li_j - 1L)* subdimI +
      (li_k - 1L)* subdimI * subdimJ

    # 8) Extract from subvol_vec
    out_vals <- subvol_vec[subvol_lin_index]

    # 9) Convert 'out_vals' into an array of dim (nI, nJ, nK)
    arr_out <- array(out_vals, dim=out_dim)

    # 10) Possibly drop dimensions
    if (drop) {
      arr_out <- drop(arr_out)
    }

    return(arr_out)
  }
)


#' Convert NeuroVol to HDF5 Format (Improved)
#'
#' @description
#' Converts a NeuroVol object to HDF5 format and saves it to a file,
#' storing enough spatial info to reconstruct a valid NeuroSpace in H5NeuroVol.
#'
#' @param vol A \code{\link[neuroim2]{NeuroVol-class}} object to convert
#' @param file_name Character string specifying the output file path
#' @param data_type String: "FLOAT", "DOUBLE", "INT", etc.
#' @param chunk_dim Optional numeric vector specifying chunk dimensions
#' @param nbit Logical; if TRUE, uses N-bit compression (default: FALSE)
#' @param compression Integer from 1-9 specifying compression level (default: 6)
#'
#' @return A \code{\link{H5NeuroVol-class}} object representing the converted volume
#' @export
to_nih5_vol <- function(vol, file_name=NULL, data_type="FLOAT",
                        chunk_dim=NULL, nbit=FALSE, compression=6)
{
  if (is.null(file_name)) {
    file_name <- paste0(tempfile(), ".h5")
  }

  # 1) Open or create HDF5 file
  h5obj <- hdf5r::H5File$new(file_name, mode="w")

  # 2) Write the main 3D dataset
  space_ds <- hdf5r::H5S$new(dims = dim(vol), maxdims = dim(vol))
  dtype <- hdf5r::H5P_DATASET_CREATE$new()

  # Match the requested data_type
  h5dtype <- switch(data_type,
                    "BINARY"=hdf5r::h5types$H5T_NATIVE_HBOOL,
                    "SHORT"=hdf5r::h5types$H5T_NATIVE_SHORT,
                    "INT"=hdf5r::h5types$H5T_NATIVE_INT,
                    "INTEGER"=hdf5r::h5types$H5T_NATIVE_INT,
                    "FLOAT"=hdf5r::h5types$H5T_NATIVE_FLOAT,
                    "DOUBLE"=hdf5r::h5types$H5T_NATIVE_DOUBLE,
                    "LONG"=hdf5r::h5types$H5T_NATIVE_LONG,
                    NULL
  )
  if (is.null(h5dtype)) {
    stop(paste("unsupported data_type:", data_type))
  }

  # Optionally define chunk size
  if (is.null(chunk_dim)) {
    # By default, chunk in Z dimension or something minimal
    chunk_dim <- c(dim(vol)[1], dim(vol)[2], 1)
  }

  dtype <- dtype$set_chunk(chunk_dim)$set_fill_value(h5dtype, 0)
  if (compression > 0) {
    # If you want compression, do so here
    dtype <- dtype$set_deflate(compression)
  }
  if (nbit) {
    # If you want N-bit, do so here
    dtype <- dtype$set_nbit()
  }

  # 3) Create group for image data
  dgroup <- h5obj$create_group("data")

  # 4) Actually write the array data
  dset <- dgroup$create_dataset(
    name             = "elements",
    space            = space_ds,
    dtype            = h5dtype,
    dataset_create_pl= dtype,
    chunk_dims       = chunk_dim
  )
  dset[,,] <- vol@.Data

  # 5) Store minimal "rtype" attribute so H5NeuroVol knows the file type
  hdf5r::h5attr(h5obj, "rtype") <- "DenseNeuroVol"

  # 6) Store the NeuroSpace info in a "space" group
  #    (dim, origin, spacing, transform, plus a marker for "3D axes")
  space_grp <- h5obj$create_group("space")
  space_grp[["dim"]]     <- dim(vol)
  space_grp[["origin"]]  <- origin(space(vol))
  space_grp[["spacing"]] <- spacing(space(vol))
  space_grp[["trans"]]   <- trans(space(vol))
  # We'll store a simple 'axes' marker
  space_grp[["axes"]]    <- c("x","y","z")

  # (You could store more info about the direction of axes if needed.)

  # 7) Return new H5NeuroVol object
  new("H5NeuroVol",
      space = space(vol),  # or a placeholder; typically replaced in constructor
      h5obj = h5obj)
}


#' Convert DenseNeuroVol to H5NeuroVol
#'
#' This function converts a DenseNeuroVol object to an H5NeuroVol object.
#'
#' @param from A DenseNeuroVol object.
#'
#' @return An H5NeuroVol object resulting from the conversion.
#'
#' @keywords internal
#' @name DenseNeuroVol,H5NeuroVol
setAs(from="DenseNeuroVol", to="H5NeuroVol", def=function(from) {
  to_nih5_vol(from, file_name=NULL, data_type="FLOAT")
})
