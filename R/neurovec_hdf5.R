#' @include all_class.R
NULL

#' H5NeuroVecSource
#'
#' @description
#' Creates a source object for building \code{\link{H5NeuroVec}} instances from an HDF5 file.
#'
#' @param file_name A \code{character} string specifying the HDF5 file name.
#' @return A new \code{H5NeuroVecSource} object (internal).
#' @keywords internal
#' @noRd
H5NeuroVecSource <- function(file_name) {
  new("H5NeuroVecSource", file_name = file_name)
}

#' H5NeuroVec Constructor
#'
#' @description
#' Constructs an \code{\link{H5NeuroVec}} object, which represents a 4D brain image
#' stored in HDF5 format.
#'
#' @param file_name The path to a 4D HDF5 file.
#'
#' @return A new \code{\link{H5NeuroVec-class}} instance.
#'
#' @seealso
#' \code{\link[neuroim2]{NeuroVec-class}} for details on 4D brain images.
#'
#' @importFrom assertthat assert_that
#' @importFrom neuroim2 NeuroSpace
#' @export
H5NeuroVec <- function(file_name) {
  assert_that(is.character(file_name))
  assert_that(file.exists(file_name))

  h5obj <- hdf5r::H5File$new(file_name)

  # Check the 'rtype' attribute
  rtype <- try(hdf5r::h5attr(h5obj, which="rtype"), silent=TRUE)
  if (!is.character(rtype) || rtype != "DenseNeuroVec") {
    stop("Invalid HDF5 file for H5NeuroVec: ", file_name)
  }

  # Check dimension count
  if (length(h5obj[["space/dim"]][]) != 4) {
    stop(
      "Cannot create H5NeuroVec: file must have 4 dimensions; found: ",
      paste(h5obj[["space/dim"]][], collapse=" ")
    )
  }

  # Build NeuroSpace
  sp <- NeuroSpace(
    dim    = h5obj[["space/dim"]][],
    origin = h5obj[["space/origin"]][],
    trans  = h5obj[["space/trans"]][,]
  )

  new("H5NeuroVec", space=sp, obj=h5obj)
}

#' Convert DenseNeuroVec to H5NeuroVec
#'
#' @description
#' Converts a \code{DenseNeuroVec} object to an \code{H5NeuroVec} by writing it to HDF5.
#'
#' @param from A \code{DenseNeuroVec} object.
#' @return An \code{H5NeuroVec} referencing the newly created HDF5 file.
#'
#' @keywords internal
#' @name DenseNeuroVec,H5NeuroVec
setAs(
  from = "DenseNeuroVec",
  to   = "H5NeuroVec",
  def  = function(from) {
    to_nih5_vec(from, file_name=NULL, data_type="FLOAT")
  }
)

#' @export
setMethod(
  f   = "series",
  signature = signature(x="H5NeuroVec", i="integer"),
  definition = function(x, i) {
    assertthat::assert_that(max(i) <= dim(x)[4])
    assertthat::assert_that(min(i) >= 1)
    x@obj[["data"]][,,, i, drop=FALSE]
  }
)

#' @export
setMethod(
  f = "series",
  signature = signature(x="H5NeuroVec", i="numeric"),
  definition = function(x, i) {
    callGeneric(x, as.integer(i))
  }
)

setMethod(
  f = "series",
  signature = signature(x="H5NeuroVec", i="matrix"),
  definition = function(x, i) {
    assertthat::assert_that(ncol(i) == 3)
    assertthat::assert_that(max(i) <= prod(dim(x)[1:3]))
    assertthat::assert_that(min(i) >= 1)
    x@obj[["data"]][i]
  }
)

#' @param drop Whether to drop dimensions of length 1
#' @rdname series-methods
#' @export
setMethod(
  f = "series",
  signature = signature(x="H5NeuroVec", i="integer"),
  definition = function(x, i, j, k, drop=TRUE) {
    if (missing(j) && missing(k)) {
      # i is a linear index; convert to (x,y,z)
      grid <- indexToGridCpp(i, dim(x)[1:3])
      callGeneric(x, grid)
    } else {
      # Possibly expand.grid approach
      assertthat::assert_that(
        length(i)==1 && length(j)==1 && length(k)==1,
        msg="Expecting single-voxel indices for i,j,k"
      )
      ret <- x@obj[["data/elements"]][i,j,k,]
      if (drop) drop(ret) else ret
    }
  }
)

#' @rdname series-methods
#' @export
setMethod(
  f = "series",
  signature = signature(x="H5NeuroVec", i="numeric"),
  definition = function(x, i, j, k) {
    if (missing(j) && missing(k)) {
      callGeneric(x, as.integer(i))
    } else {
      callGeneric(x, as.integer(i), as.integer(j), as.integer(k))
    }
  }
)

#' @rdname series-methods
#' @export
setMethod(
  f = "series",
  signature = signature(x="H5NeuroVec", i="matrix"),
  definition = function(x, i) {
    assertthat::assert_that(ncol(i) == 3)
    d4 <- dim(x)[4]

    # Build bounding box for i
    ir <- lapply(seq_len(ncol(i)), function(j) seq(min(i[,j]), max(i[,j])))

    # e.g. sub-block
    ret <- x@obj[["data/elements"]][
      ir[[1]][1]:ir[[1]][length(ir[[1]])],
      ir[[2]][1]:ir[[2]][length(ir[[2]])],
      ir[[3]][1]:ir[[3]][length(ir[[3]])],
      , drop=FALSE
    ]

    # flatten
    ret2 <- t(array(ret, c(prod(dim(ret)[1:3]), dim(ret)[4])))
    # check if shape matches nrow(i)
    if (ncol(ret2) != nrow(i)) {
      i2 <- apply(i, 2, function(ind) {
        ind - min(ind) + 1
      })
      i3 <- gridToIndex3DCpp(dim(ret)[1:3], i2)
      ret2[, i3, drop=FALSE]
    } else {
      ret2
    }
  }
)

#' A bounding-box approach for linear indexing of a 4D H5NeuroVec
#'
#' @description
#' Converts linear indices \code{i} into 4D subscripts \code{(x,y,z,t)} in
#' \code{dim(x)}, then reads the minimal bounding region from HDF5.
#'
#' @param x A \code{H5NeuroVec} object (4D).
#' @param i Numeric vector of linear indices into the 4D array.
#' @return A numeric vector of length \code{length(i)}, in the same order as \code{i}.
#'
#' @export
setMethod(
  f = "linear_access",
  signature = signature(x="H5NeuroVec", i="numeric"),
  definition = function(x, i) {
    # 1) Convert linear => (x,y,z,t)
    coords <- arrayInd(i, dim(x))
    n <- nrow(coords)

    # 2) bounding box
    minX <- min(coords[,1]); maxX <- max(coords[,1])
    minY <- min(coords[,2]); maxY <- max(coords[,2])
    minZ <- min(coords[,3]); maxZ <- max(coords[,3])
    minT <- min(coords[,4]); maxT <- max(coords[,4])

    # 3) Read sub-block
    dset <- x@obj[["data/elements"]]
    sub4d <- dset[minX:maxX, minY:maxY, minZ:maxZ, minT:maxT, drop=FALSE]

    # 4) Flatten sub4d
    sub4d_vec <- as.vector(sub4d)

    # 5) Compute local offsets
    offX <- coords[,1] - minX + 1
    offY <- coords[,2] - minY + 1
    offZ <- coords[,3] - minZ + 1
    offT <- coords[,4] - minT + 1

    lenX <- maxX - minX + 1
    lenY <- maxY - minY + 1
    lenZ <- maxZ - minZ + 1

    sub_lin_idx <- offX +
      (offY - 1L)* lenX +
      (offZ - 1L)* lenX*lenY +
      (offT - 1L)* lenX*lenY*lenZ

    out_vals <- sub4d_vec[sub_lin_idx]
    out_vals
  }
)

#' @export
setMethod(
  f = "linear_access",
  signature = signature(x="H5NeuroVec", i="integer"),
  definition = function(x, i) {
    callGeneric(x, as.numeric(i))
  }
)

#' Subset a 4D H5NeuroVec with arbitrary numeric i, j, k, l (bounding box)
#'
#' @description
#' Handles \code{x[i,j,k,l]} subsetting. We read a bounding hyper-slab from HDF5 once,
#' then pick out the requested points in correct order, avoiding loops in R.
#'
#' @param x \code{H5NeuroVec} (4D).
#' @param i,j,k,l Numeric vectors of indices (possibly non-contiguous).
#' @param drop Logical; drop dimensions of size 1?
#' @param ... Ignored.
#' @return An array of dimension \code{length(i) × length(j) × length(k) × length(l)},
#'   dropped if \code{drop=TRUE}.
#'
#' @export
setMethod(
  f = "[",
  signature = signature(x="H5NeuroVec", i="numeric", j="numeric", drop="ANY"),
  definition = function(x, i, j, k, l, ..., drop=TRUE) {
    i <- as.numeric(i)
    j <- as.numeric(j)

    if (missing(k)) {
      k <- seq_len(dim(x)[3])
    }
    if (missing(l)) {
      l <- seq_len(dim(x)[4])
    }

    dims_x <- dim(x)
    if (any(i<1|i>dims_x[1]) || any(j<1|j>dims_x[2]) ||
        any(k<1|k>dims_x[3]) || any(l<1|l>dims_x[4])) {
      stop("Subscript out of range in H5NeuroVec dimensions.")
    }

    minI <- min(i); maxI <- max(i)
    minJ <- min(j); maxJ <- max(j)
    minK <- min(k); maxK <- max(k)
    minL <- min(l); maxL <- max(l)

    dset <- x@obj[["data/elements"]]
    subvol <- dset[minI:maxI, minJ:maxJ, minK:maxK, minL:maxL, drop=FALSE]

    i_off <- i - minI + 1
    j_off <- j - minJ + 1
    k_off <- k - minK + 1
    l_off <- l - minL + 1

    subdimI <- maxI - minI + 1
    subdimJ <- maxJ - minJ + 1
    subdimK <- maxK - minK + 1

    subvol_vec <- as.vector(subvol)

    nI <- length(i)
    nJ <- length(j)
    nK <- length(k)
    nL <- length(l)

    out_dim <- c(nI, nJ, nK, nL)
    N <- nI * nJ * nK * nL

    idx_i <- rep(seq_len(nI), times=nJ*nK*nL)
    idx_j <- rep(rep(seq_len(nJ), each=nI), times=nK*nL)
    idx_k <- rep(rep(seq_len(nK), each=nI*nJ), times=nL)
    idx_l <- rep(seq_len(nL), each=nI*nJ*nK)

    loc_i <- i_off[idx_i]
    loc_j <- j_off[idx_j]
    loc_k <- k_off[idx_k]
    loc_l <- l_off[idx_l]

    sub_lin_idx <- loc_i +
      (loc_j-1)* subdimI +
      (loc_k-1)* subdimI*subdimJ +
      (loc_l-1)* subdimI*subdimJ*subdimK

    out_vals <- subvol_vec[sub_lin_idx]

    arr_out <- array(out_vals, dim=out_dim)
    if (drop) {
      arr_out <- drop(arr_out)
    }
    arr_out
  }
)

#' Load data from an H5NeuroVecSource object
#'
#' @description
#' Loads a \code{\link{H5NeuroVec}} from an \code{H5NeuroVecSource},
#' effectively calling \code{H5NeuroVec()} on the stored file path.
#'
#' @param x A \code{H5NeuroVecSource} with \code{file_name}.
#' @return A new \code{H5NeuroVec}.
#'
#' @seealso \code{\link{H5NeuroVecSource}}, \code{\link{H5NeuroVec}}
#'
#' @noRd
setMethod(
  f = "load_data",
  signature = c("H5NeuroVecSource"),
  definition = function(x) {
    H5NeuroVec(x@file_name)
  }
)

#' Convert NeuroVec to HDF5 Format, returning an H5NeuroVec
#'
#' @description
#' Creates an HDF5 file from a 4D \code{\link[neuroim2]{NeuroVec}} object, storing
#' data & spatial info. Returns an \code{\link{H5NeuroVec}} referencing the file
#' for on-demand access.
#'
#' @param vec A \code{NeuroVec} (4D) to convert (coerced internally to \code{DenseNeuroVec}).
#' @param file_name Path for the output HDF5 file (if NULL, a temp file is used).
#' @param data_type Storage type (e.g., "FLOAT"). Default "FLOAT".
#' @param chunk_dim Chunk dimensions. Default c(4,4,4,dim(vec)[4]).
#' @param nbit \code{logical}: use N-bit filter? Default \code{FALSE}.
#' @param compression \code{integer} [0..9], default 6.
#'
#' @return An \code{H5NeuroVec} referencing the new HDF5 file.
#'
#' @details
#' 1. Coerces \code{vec} to \code{DenseNeuroVec}
#' 2. Writes data to \code{/data/elements} in HDF5
#' 3. Writes dimension/origin/spacing/transform in \code{/space}
#' 4. Tags file attribute "rtype" => "DenseNeuroVec"
#' 5. Returns \code{H5NeuroVec}
#'
#' @keywords internal
to_nih5_vec <- function(vec,
                        file_name   = NULL,
                        data_type   = "FLOAT",
                        chunk_dim   = c(4,4,4, dim(vec)[4]),
                        nbit        = FALSE,
                        compression = 6)
{
  if (!requireNamespace("hdf5r", quietly=TRUE)) {
    stop("Package 'hdf5r' must be installed for HDF5 I/O.", call.=FALSE)
  }

  # 1) Coerce to DenseNeuroVec
  vec <- as(vec, "DenseNeuroVec")

  # 2) Check compression
  assert_that(compression >= 0 && compression <= 9)

  # 3) If no file_name, use temp
  if (is.null(file_name)) {
    file_name <- tempfile(fileext=".h5")
  }

  # 4) Create or overwrite
  h5obj <- hdf5r::H5File$new(file_name, mode="w")

  # 5) Setup dimension
  space_ds <- hdf5r::H5S$new(dims=dim(vec), maxdims=dim(vec))
  dtype_pl  <- hdf5r::H5P_DATASET_CREATE$new()

  # 6) Data type
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
    stop("Unsupported 'data_type': ", data_type)
  }

  # 7) Chunk/fill/compression
  dtype_pl$set_chunk(chunk_dim)$set_fill_value(h5dtype,0)$set_deflate(compression)
  if (nbit && compression>0) {
    dtype_pl$set_nbit()
  }

  # 8) Tag the file
  hdf5r::h5attr(h5obj, "rtype") <- "DenseNeuroVec"

  # 9) /data & /space
  dgroup <- h5obj$create_group("data")
  sgroup <- h5obj$create_group("space")

  # 10) dataset "elements"
  dset <- dgroup$create_dataset(
    name             = "elements",
    space            = space_ds,
    dtype            = h5dtype,
    dataset_create_pl= dtype_pl,
    chunk_dims       = chunk_dim,
    gzip_level       = compression
  )

  # 11) space metadata
  sgroup[["dim"]]     <- dim(vec)
  sgroup[["spacing"]] <- spacing(vec)
  sgroup[["origin"]]  <- origin(space(vec))
  sgroup[["trans"]]   <- trans(space(vec))

  # 12) Write data
  dset[,,,] <- as.array(vec)

  # 13) Return new H5NeuroVec
  new("H5NeuroVec", space=space(vec), obj=h5obj)
}
