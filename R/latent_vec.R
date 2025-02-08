#' @importFrom neuroim2 matricized_access
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
  if (!inherits(mask, "LogicalNeuroVol")) {
    mspace <- NeuroSpace(dim(space)[1:3],
                         spacing(space),
                         origin(space),
                         axes(space),
                         trans(space))
    mask <- LogicalNeuroVol(as.logical(mask), mspace)
  }

  cardinality <- sum(mask)

  # Handle offset
  if (is.null(offset)) {
    offset <- rep(0, nrow(loadings))
  } else if (length(offset) != nrow(loadings)) {
    stop("'offset' length must match number of voxels in 'loadings'")
  }

  # Dimension checks
  if (nrow(loadings) != cardinality) {
    stop("'loadings' must have ", cardinality, " rows (i.e. #nonzero in mask)")
  }
  if (ncol(loadings) != ncol(basis)) {
    stop("'basis' and 'loadings' must have the same number of columns")
  }
  if (nrow(basis) != dim(space)[4]) {
    stop("'basis' must have ", dim(space)[4], " rows (the 4th dimension of space)")
  }

  # Convert to sparse if needed
  basis    <- if (is.matrix(basis)) Matrix(basis) else basis
  loadings <- if (is.matrix(loadings)) Matrix(loadings) else loadings

  # Create the object
  new("LatentNeuroVec",
      basis    = basis,
      loadings = loadings,
      space    = space,
      mask     = mask,
      map      = IndexLookupVol(space(mask), as.integer(which(mask))),
      offset   = offset,
      label    = label)
}

#' Write LatentNeuroVec to HDF5 File
#'
#' @description
#' Writes a \code{LatentNeuroVec} to an HDF5 file with optional compression.
#'
#' @param x A \code{\link[neuroim2]{LatentNeuroVec-class}} to write.
#' @param file_name \code{character} file path to the output HDF5.
#' @param nbit \code{logical}; if TRUE, uses N-bit compression (default: FALSE).
#' @param compression \code{integer} in [1..9] specifying compression level (default: 9).
#' @param chunk_dim Optional numeric vector specifying chunk dimensions.
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
#' # Write a LatentNeuroVec
#' write_vec(lvec, "brain_components.h5", compression=6)
#' }
#'
#' @seealso
#' \code{\link[neuroim2]{read_vec}} for reading back the file.
#'
#' @rdname write_vec-methods
#' @export
setMethod(
  f = "write_vec",
  signature = signature(x="LatentNeuroVec", file_name="character", format="missing", data_type="missing"),
  definition = function(x, file_name, nbit=FALSE, compression=9, chunk_dim=NULL) {
    if (!is.character(file_name) || length(file_name) != 1) {
      stop("'file_name' must be a single character string")
    }
    if (!is.numeric(compression) || compression < 1 || compression > 9) {
      stop("'compression' must be an integer between 1 and 9")
    }

    obj <- to_h5_latentvec(
      vec        = x,
      file_name  = file_name,
      chunk_dim  = chunk_dim,
      nbit       = nbit,
      compression= compression
    )
    obj$close()
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
    if (!is.numeric(i) || ncol(i) != 2) {
      stop("Index matrix 'i' must be numeric with 2 columns")
    }
    b1 <- x@basis[as.numeric(i[,1]),, drop=FALSE]
    b2 <- x@loadings[as.numeric(i[,2]),, drop=FALSE]
    rowSums(b1 * b2) + x@offset[i[,2]]
  }
)

#' @keywords internal
#' @noRd
#' @importFrom neuroim2 matricized_access
setMethod(
  f = "matricized_access",
  signature = signature(x="LatentNeuroVec", i="integer"),
  definition = function(x, i) {
    if (any(i < 1) || any(i > nrow(x@loadings))) {
      stop("Index out of bounds for 'loadings'")
    }
    b1 <- x@basis
    b2 <- x@loadings[as.numeric(i),, drop=FALSE]
    out <- tcrossprod(b1, b2)
    # Add offsets
    as.matrix(sweep(out, 2, x@offset[i], "+"))
  }
)

#' @keywords internal
#' @noRd
setMethod(
  f = "matricized_access",
  signature = signature(x="LatentNeuroVec", i="numeric"),
  definition = function(x, i) {
    callGeneric(x, as.integer(i))
  }
)

#' Internal Linear Access Method for LatentNeuroVec
#'
#' @description
#' Internal method providing linear access to elements.
#'
#' @param x A \code{LatentNeuroVec}.
#' @param i A numeric vector of indices.
#'
#' @return Computed values at those linear indices.
#' @keywords internal
#' @noRd
setMethod(
  f = "linear_access",
  signature = signature(x="LatentNeuroVec", i="numeric"),
  definition = function(x, i) {
    callGeneric(x, i)
  }
)

#' Internal Linear Access Method for LatentNeuroVec
#'
#' @description
#' Internal method providing linear access to elements.
#'
#' @param x A \code{LatentNeuroVec}.
#' @param i A numeric vector of indices.
#'
#' @return Computed values
#' @keywords internal
#' @noRd
setMethod(
  f = "linear_access",
  signature = signature(x="LatentNeuroVec", i="integer"),
  definition = function(x, i) {
    browser()  # presumably for debugging
    if (!is.numeric(i) || any(is.na(i))) {
      stop("Index 'i' must be numeric without NA values")
    }
    nels <- prod(dim(x)[1:3])
    if (any(i < 1) || any(i > nels * dim(x)[4])) {
      stop("Index out of bounds for 4D volume")
    }

    n <- ceiling(i / nels)
    offset <- i %% nels
    offset[offset == 0] <- nels

    ll <- lookup(x@map, offset)
    nz <- which(ll > 0)
    if (length(nz) == 0) {
      return(numeric(length(i)))
    }

    idx2d <- cbind(n[nz], ll[nz])
    b1 <- x@basis[idx2d[,1],, drop=FALSE]
    b2 <- x@loadings[idx2d[,2],, drop=FALSE]
    vals <- rowSums(b1 * b2) + x@offset[idx2d[,2]]

    ovals <- numeric(length(i))
    ovals[nz] <- vals
    ovals
  }
)

#' Extract a Single Volume from LatentNeuroVec
#'
#' @description
#' Extracts a single volume from a \code{LatentNeuroVec} as a \code{SparseNeuroVol}.
#'
#' @param x A \code{\link[neuroim2]{LatentNeuroVec-class}} object.
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
#' @rdname extract-methods
#' @export
setMethod(
  f = "[[",
  signature = signature(x="LatentNeuroVec", i="numeric"),
  definition = function(x, i) {
    if (length(i) != 1) {
      stop("Index must be a single number")
    }
    if (i < 1 || i > dim(x)[4]) {
      stop("Index out of bounds")
    }

    xs <- space(x)
    dat <- (tcrossprod(x@basis[i,,drop=FALSE], x@loadings))[1,]
    dat <- dat + x@offset

    newdim <- dim(xs)[1:3]
    bspace <- NeuroSpace(newdim,
                         spacing=spacing(xs),
                         origin=origin(xs),
                         axes=axes(xs),
                         trans=trans(xs))

    SparseNeuroVol(dat, bspace, indices=indices(x))
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
#' @rdname concat-methods
#' @export
setMethod(
  f = "concat",
  signature = signature(x="LatentNeuroVec", y="LatentNeuroVec"),
  definition = function(x, y, ...) {
    do.call(NeuroVecSeq, list(x, y, ...))
  }
)

#' Convert LatentNeuroVec to HDF5 Format
#'
#' @description
#' Internal function to convert a \code{LatentNeuroVec} to HDF5 format.
#'
#' @param vec A \code{LatentNeuroVec} object to be saved.
#' @param file_name Path for the output HDF5 file. If NULL, a temp file is used.
#' @param data_type String specifying data type (e.g., "FLOAT"). Default is "FLOAT".
#' @param chunk_dim Numeric vector for chunking. If NULL, no chunking.
#' @param nbit Logical; if TRUE, use N-bit filter (default FALSE).
#' @param compression Integer in [0..9], default 6.
#'
#' @return An \code{hdf5r} \code{H5File} object.
#'
#' @details
#' Creates an HDF5 file with \code{basis}, \code{loadings}, \code{offset},
#' plus minimal space/mask metadata.
#'
#' @keywords internal
#' @noRd
to_h5_latentvec <- function(vec, file_name=NULL, data_type="FLOAT",
                            chunk_dim=NULL, nbit=FALSE, compression=6) {

  assert_that(inherits(vec, "LatentNeuroVec"))

  if (!is.null(file_name) && !endsWith(file_name, ".lv.h5")) {
    file_name <- paste0(file_name, ".lv.h5")
  } else if (is.null(file_name)) {
    file_name <- tempfile(fileext = ".lv.h5")
  }

  h5obj <- hdf5r::H5File$new(file_name)
  hdf5r::h5attr(h5obj, "rtype") <- class(vec)

  dgroup <- h5obj$create_group("data")
  dspace <- h5obj$create_group("space")

  dspace[["dim"]]     <- dim(vec)
  dspace[["spacing"]] <- spacing(vec)
  dspace[["origin"]]  <- origin(space(vec))
  dspace[["trans"]]   <- trans(space(vec))

  dtype <- hdf5r::H5P_DATASET_CREATE$new()
  if (nbit && compression > 0) {
    dtype <- dtype$set_nbit()
  }

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

  if (!is.null(chunk_dim)) {
    dtype <- dtype$set_chunk(chunk_dim)$set_fill_value(h5dtype, 0)$set_deflate(compression)
  } else {
    dtype <- dtype$set_fill_value(h5dtype, 0)$set_deflate(compression)
  }

  basis_ds    <- hdf5r::H5S$new(dims=dim(vec@basis),    maxdims=dim(vec@basis))
  loadings_ds <- hdf5r::H5S$new(dims=dim(vec@loadings), maxdims=dim(vec@loadings))

  basis_dset <- dgroup$create_dataset(
    name             = "basis",
    space            = basis_ds,
    dtype            = h5dtype,
    dataset_create_pl= dtype,
    gzip_level       = compression
  )

  loadings_dset <- dgroup$create_dataset(
    name             = "loadings",
    space            = loadings_ds,
    dtype            = h5dtype,
    dataset_create_pl= dtype,
    gzip_level       = compression
  )

  basis_dset[,]    <- as.matrix(vec@basis)
  loadings_dset[,] <- as.matrix(vec@loadings)
  dgroup[["offset"]]  <- vec@offset
  dgroup[["indices"]] <- as.integer(vec@map@indices)

  h5obj
}

#' Load data from a LatentNeuroVecSource object
#'
#' @description
#' Constructs a \code{LatentNeuroVec} from a \code{LatentNeuroVecSource} by reading
#' basis/loadings/offset info in an HDF5 file.
#'
#' @param x A \code{LatentNeuroVecSource} specifying the HDF5 file path.
#'
#' @return A \code{LatentNeuroVec} with loaded data.
#'
#' @details
#' 1. Opens file
#' 2. Reads basis, loadings, offset, indices
#' 3. Builds a \code{NeuroSpace}
#' 4. Creates a mask from \code{indices}
#' 5. Returns \code{LatentNeuroVec}
#'
#' @note Requires \pkg{hdf5r} for file I/O.
#' @seealso \code{\link{LatentNeuroVecSource}}, \code{\link{LatentNeuroVec}}
#'
#' @importFrom hdf5r H5File
#' @noRd
setMethod(
  f = "load_data",
  signature = c("LatentNeuroVecSource"),
  definition = function(x) {
    h5obj <- hdf5r::H5File$new(x@file_name)
    basis    <- h5obj[["data/basis"]][,]
    loadings <- h5obj[["data/loadings"]][,]
    offset   <- h5obj[["data/offset"]][]
    indices  <- h5obj[["data/indices"]][]

    sp <- NeuroSpace(
      dim     = h5obj[["space/dim"]][],
      spacing = h5obj[["space/spacing"]][],
      origin  = h5obj[["space/origin"]][],
      trans   = h5obj[["space/trans"]][,]
    )

    mask_vol <- NeuroVol(array(0, dim(sp)[1:3]), drop_dim(sp))
    mask_vol[indices] <- 1
    mask_vol <- as.logical(mask_vol)

    h5obj$close_all()

    LatentNeuroVec(basis, loadings, space=sp, mask=mask_vol, offset=offset,
                   label=basename(x@file_name))
  }
)

#' Show Method for LatentNeuroVec Objects
#'
#' @description
#' Displays a formatted summary of a \code{LatentNeuroVec} with colored output.
#'
#' @param object A \code{LatentNeuroVec} to display.
#'
#' @details
#' Uses the \pkg{crayon} package for colored terminal output, showing:
#' \itemize{
#'   \item Dimensions & components
#'   \item Memory usage
#'   \item Sparsity info
#'   \item Space metadata
#' }
#'
#' @importFrom crayon bold red green blue yellow silver
#' @importFrom utils object.size
#' @keywords internal
setMethod(
  f = "show",
  signature = "LatentNeuroVec",
  definition = function(object) {
    # Header
    cat("\n", crayon::bold(crayon::blue("LatentNeuroVec Object")), "\n")
    cat(crayon::silver("══════════════════════\n"))

    # Dimensions
    dims <- dim(object)
    spatial_dims <- paste(dims[1:3], collapse=" × ")
    cat("\n", crayon::yellow("Dimensions:"), "\n")
    cat(" ", crayon::silver("•"), " Spatial: ", crayon::green(spatial_dims), "\n")
    cat(" ", crayon::silver("•"), " Temporal: ", crayon::green(dims[4]), "\n")

    # Components
    n_components <- ncol(object@basis)
    variance_explained <- format(object@basis[1:min(5, nrow(object@basis)), 1], digits=3)
    cat("\n", crayon::yellow("Components:"), "\n")
    cat(" ", crayon::silver("•"), " Number: ", crayon::green(n_components), "\n")
    cat(" ", crayon::silver("•"), " First coefficients: ",
        crayon::green(paste(variance_explained, collapse=", ")),
        if (length(variance_explained) < 5) "" else "...", "\n")

    # Memory Usage
    basis_size    <- format(object.size(object@basis),    units="auto")
    loadings_size <- format(object.size(object@loadings), units="auto")
    total_size    <- format(object.size(object),          units="auto")

    cat("\n", crayon::yellow("Memory Usage:"), "\n")
    cat(" ", crayon::silver("•"), " Basis: ",    crayon::green(basis_size),    "\n")
    cat(" ", crayon::silver("•"), " Loadings: ", crayon::green(loadings_size), "\n")
    cat(" ", crayon::silver("•"), " Total: ",    crayon::green(total_size),    "\n")

    # Sparsity
    n_nonzero <- sum(object@mask)
    sparsity <- round(100 * n_nonzero / prod(dims[1:3]), 2)
    cat("\n", crayon::yellow("Sparsity:"), "\n")
    cat(" ", crayon::silver("•"), " Non-zero voxels: ", crayon::green(n_nonzero), "\n")
    cat(" ", crayon::silver("•"), " Coverage: ",       crayon::green(sparsity), "%\n")

    # Space Info
    sp <- space(object)
    spacing_str <- paste(round(spacing(sp), 2), collapse=" × ")
    origin_str  <- paste(round(origin(sp), 2), collapse=" × ")
    cat("\n", crayon::yellow("Space Information:"), "\n")
    cat(" ", crayon::silver("•"), " Spacing: ", crayon::green(spacing_str), "\n")
    cat(" ", crayon::silver("•"), " Origin:  ", crayon::green(origin_str),  "\n")

    # Footer
    cat("\n", crayon::bold("Data Access:"), "\n")
    cat("\n", crayon::yellow("Reconstructed Space Access:"), "\n")
    cat(" ", crayon::silver("•"), " Extract volume: ",
        crayon::blue("object[[i]]"),
        crayon::silver("  # 3D volume at timepoint i\n"))
    cat(" ", crayon::silver("•"), " Get value: ",
        crayon::blue("object[i]"),
        crayon::silver("  # Value at linear index i\n"))
    cat(" ", crayon::silver("•"), " Subset: ",
        crayon::blue("object[mask]"),
        crayon::silver("  # Values at mask positions\n"))

    cat("\n", crayon::yellow("Latent Space Access:"), "\n")
    cat(" ", crayon::silver("•"), " Basis vectors: ",
        crayon::blue("basis(object)"),
        crayon::silver("  # n×k temporal basis\n"))
    cat(" ", crayon::silver("•"), " Loadings: ",
        crayon::blue("loadings(object)"),
        crayon::silver("  # p×k spatial loadings\n"))
    cat(" ", crayon::silver("•"), " Components: ",
        crayon::blue("components(object)"),
        crayon::silver("  # List of k component volumes\n"))

    cat("\n", crayon::yellow("Conversions:"), "\n")
    cat(" ", crayon::silver("•"), " as.array(object): ",
        crayon::silver("4D reconstruction\n"))
    cat(" ", crayon::silver("•"), " as.matrix(object): ",
        crayon::silver("n×p matrix of reconstructed values\n"))

    cat("\n", crayon::silver("Note: All access methods reconstruct data (X = B × L^T + offset)"),
        "\n", crayon::silver("unless you're explicitly accessing latent space."), "\n\n")
  }
)
