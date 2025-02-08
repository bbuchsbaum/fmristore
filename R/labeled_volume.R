#' Write a NeuroVec to an HDF5 file with NIfTI-like quaternions
#'
#' @description
#' Creates an HDF5 file following a NIfTI-like header layout, storing:
#' \itemize{
#'   \item \code{/header/dim} => \code{[4, X, Y, Z, nVols, 1,1,1]}
#'   \item \code{/header/pixdim} => \code{[qfac, dx, dy, dz, ...]}
#'   \item \code{/header/quatern_b,c,d} and \code{qoffset_x,y,z}
#'   \item \code{/mask} => 3D dataset \code{[X, Y, Z]} (0/1)
#'   \item \code{/header/labels} => array of label strings
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
#' We also gather voxel spacing (dx,dy,dz) from \code{spacing(space(vec))} and
#' the origin from \code{origin(space(vec))}.
#'
#' We then store a minimal set of NIfTI-like header fields in the
#' \code{/header} group. The user can supply \code{header_values} (a named
#' list) to override or augment additional fields (e.g., \code{qform_code=1L}).
#'
#' @param vec A 4D \code{\link[neuroim2]{NeuroVec}} with dimension \code{[X,Y,Z,nVols]}.
#' @param mask A \code{\link[neuroim2]{LogicalNeuroVol}} of shape \code{[X,Y,Z]}
#'   (the same 3D shape as \code{vec}).
#' @param labels A character vector of length \code{nVols}, labeling each 4D sub-volume.
#' @param file Either a character path to the HDF5 file to create or
#'   an open \code{\link[hdf5r]{H5File}} in write mode.
#' @param compression Integer \code{0-9} for gzip level; default \code{4}.
#' @param dtype An HDF5 data type object (e.g., \code{h5types$H5T_NATIVE_FLOAT})
#'   or a list of types (one per volume). Default is \code{h5types$H5T_NATIVE_DOUBLE}.
#' @param chunk_size If non-NULL, the chunk dimension for the 1D datasets. Default is \code{1024}.
#' @param header_values A named list of optional overrides for fields in the header
#'   (e.g., \code{list(qform_code=1L, sform_code=2L)}).
#'
#' @return Invisibly returns the \code{\link[hdf5r]{H5File}} object, which now
#'   has all the data and header info stored in it.
#'
#' @seealso
#' \code{\link[neuroim2]{matrixToQuatern}} for how the quaternion is derived,
#' \code{\link[neuroim2]{quaternToMatrix}} for reconstructing the 4×4,
#' \code{\link{read_labeled_vec}} for reading the file back in.
#'
#' @import hdf5r
#' @importFrom neuroim2 spacing space origin trans matrixToQuatern
#' @export
write_labeled_vec <- function(vec,
                              mask,
                              labels,
                              file,
                              compression = 4,
                              dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE,
                              chunk_size = 1024,
                              header_values = list())
{
  newly_opened <- FALSE
  if (is.character(file)) {
    file_path <- file
    file <- hdf5r::H5File$new(file_path, mode = "w")
    newly_opened <- TRUE
    on.exit({
      if (newly_opened && file$is_valid) file$close_all()
    }, add = TRUE)
  }
  stopifnot(inherits(file, "H5File"))

  # Validate shapes
  nd <- dim(vec)  # [X, Y, Z, nVols]
  stopifnot(length(nd) == 4)
  nVols <- nd[4]
  if (length(labels) != nVols) {
    stop("Length of 'labels' must match the 4th dimension of 'vec'.")
  }
  stopifnot(inherits(mask, "LogicalNeuroVol"))
  stopifnot(all(dim(mask) == nd[1:3]))

  X <- nd[1]; Y <- nd[2]; Z <- nd[3]

  # Extract 4×4 transformation matrix from NeuroSpace
  tmat <- trans(space(vec))          # e.g. a 4×4
  # Convert to quaternion + qfac
  q    <- matrixToQuatern(tmat)      # => list(quaternion=c(b,c,d), qfac=±1)

  # Gather spacing & origin
  sp  <- spacing(space(vec))         # c(dx, dy, dz)
  org <- origin(space(vec))          # c(ox, oy, oz)
  if (length(sp) < 3) {
    sp  <- c(sp, rep(1, 3 - length(sp)))
  }
  if (length(org) < 3) {
    org <- c(org, rep(0, 3 - length(org)))
  }

  # Build minimal NIfTI-like header fields
  hdr_default <- list(
    dim         = c(4L, X, Y, Z, nVols, 1L, 1L, 1L),
    pixdim      = c(q$qfac, sp[1], sp[2], sp[3], 0, 0, 0, 0),
    quatern_b   = q$quaternion[1],
    quatern_c   = q$quaternion[2],
    quatern_d   = q$quaternion[3],
    qoffset_x   = org[1],
    qoffset_y   = org[2],
    qoffset_z   = org[3],
    sizeof_hdr  = 348L,
    magic       = "n+1"
    # Additional fields if desired: qform_code, sform_code, etc.
  )

  # Merge user overrides from header_values
  for (nm in names(header_values)) {
    hdr_default[[nm]] <- header_values[[nm]]
  }

  # 1) Write these fields into /header
  hdr_grp <- file$create_group("header")
  .write_nifti_header_fields(hdr_grp, hdr_default)

  # 2) Write /mask => shape [X, Y, Z]
  mask_arr <- array(as.integer(mask@.Data), dim = c(X, Y, Z))
  dset_mask <- file$create_dataset("mask",
                                   dims  = dim(mask_arr),
                                   dtype = hdf5r::h5types$H5T_STD_U8LE)
  dset_mask[,,] <- mask_arr

  # 3) Write labels => /header/labels (variable-length string array)
  lbl_space <- hdf5r::H5S$new(dims = length(labels))
  str_vlen  <- hdf5r::H5T_STRING$new(size = Inf)
  lbl_dset  <- hdr_grp$create_dataset("labels", space = lbl_space, dtype = str_vlen)
  lbl_dset[] <- labels

  # 4) /data => subdatasets for each volume, storing masked data
  data_grp <- file$create_group("data")
  idx_nonzero <- which(mask_arr == 1)
  n_nonzero   <- length(idx_nonzero)

  # unify dtype => list
  if (!is.list(dtype)) {
    dtypes <- rep(list(dtype), nVols)
  } else {
    dtypes <- dtype
    stopifnot(length(dtypes) == nVols,
              "If 'dtype' is a list, it must match the #vols.")
  }

  chunk_dim <- if (!is.null(chunk_size)) min(chunk_size, n_nonzero) else NULL

  for (i in seq_len(nVols)) {
    message("Writing label: ", labels[i])
    vol_3d <- as.array(vec[,,, i])
    vol_1d <- vol_3d[idx_nonzero]

    vol_space <- hdf5r::H5S$new(dims = n_nonzero)
    ds <- data_grp$create_dataset(name  = labels[i],
                                  space = vol_space,
                                  dtype = dtypes[[i]],
                                  chunk = chunk_dim,
                                  gzip_level = compression)
    ds[] <- vol_1d
  }

  invisible(file)
}


#' Internal helper to write each field from a list into the HDF5 header
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
    # store single string dataset
    stype <- hdf5r::H5T_STRING$new(size = Inf)
    ds    <- hdr_grp$create_dataset(nm, dims=1, dtype=stype)
    ds[]  <- if (length(val) > 1) paste(val, collapse=" ") else val
  } else if (is.integer(val)) {
    ds <- hdr_grp$create_dataset(nm, dims=length(val),
                                 dtype=hdf5r::h5types$H5T_NATIVE_INT32)
    ds[] <- val
  } else if (is.numeric(val)) {
    ds <- hdr_grp$create_dataset(nm, dims=length(val),
                                 dtype=hdf5r::h5types$H5T_NATIVE_DOUBLE)
    ds[] <- val
  } else {
    # fallback => store as string
    stype <- hdf5r::H5T_STRING$new(size = Inf)
    ds    <- hdr_grp$create_dataset(nm, dims=1, dtype=stype)
    ds[]  <- as.character(val)
  }
}

#' Internal utility used inside write_labeled_vec (and read_labeled_vec) to store integer dataset
#' @noRd
.write_ds_int <- function(hdr_grp, name, vec) {
  ds <- hdr_grp$create_dataset(name,
                               dims=length(vec),
                               dtype=hdf5r::h5types$H5T_NATIVE_INT32)
  ds[] <- as.integer(vec)
}

#' Read a labeled volume set from an HDF5 file, restoring NeuroSpace from quaternion
#'
#' @description
#' Reads an HDF5 file conforming to the "Brain Image HDF5" spec (similar to NIfTI).
#' It reconstructs a minimal \code{\link[neuroim2]{NeuroSpace}} from the stored
#' quaternion parameters, builds a \code{LogicalNeuroVol} mask, and returns an S4 object
#' (e.g. \code{LabeledVolumeSet}) for lazy access to each 4D sub-volume.
#'
#' @details
#' Expected datasets:
#' \itemize{
#'   \item \code{/header/dim} => e.g. \code{[4, X, Y, Z, nVols, 1,1,1]}
#'   \item \code{/header/pixdim} => e.g. \code{[qfac, dx, dy, dz, ...]}
#'   \item \code{/header/quatern_b}, \code{quatern_c}, \code{quatern_d}, \code{qoffset_x}, etc.
#'   \item \code{/mask} => 3D [X,Y,Z]
#'   \item \code{/header/labels} => array of strings
#'   \item \code{/data/<label>} => 1D array (length = #nonzero in mask)
#' }
#'
#' The quaternion expansions are done via \code{\link[neuroim2]{quaternToMatrix}}.
#'
#' @param file Character path or open \code{\link[hdf5r]{H5File}} in read mode.
#' @param memoise Logical. If TRUE, caches volumes in memory after first load.
#' @return A \code{\link{LabeledVolumeSet}} referencing the open file. It contains:
#' \itemize{
#'   \item \code{mask}: A \code{LogicalNeuroVol}
#'   \item \code{labels}: The sub-volume labels
#'   \item \code{space}: A \code{NeuroSpace} capturing orientation
#'   \item An internal \code{environment} for lazy load
#' }
#'
#' @import hdf5r
#' @importFrom neuroim2 quaternToMatrix
#' @export
read_labeled_vec <- function(file, memoise=FALSE) {
  if (is.character(file)) {
    file <- hdf5r::H5File$new(file, mode="r")
  }
  stopifnot(inherits(file, "H5File"))

  hdr_grp <- file[["header"]]
  if (is.null(hdr_grp)) stop("No '/header' group found in file")

  # Helper to read dataset if present
  .rd <- function(nm) {
    d <- hdr_grp[[nm]]
    if (!is.null(d)) d[] else NULL
  }

  dims       <- .rd("dim")         # c(4, X, Y, Z, nVols, 1,1,1)
  pixdim     <- .rd("pixdim")      # c(qfac, dx, dy, dz, ...)
  qb         <- .rd("quatern_b")
  qc         <- .rd("quatern_c")
  qd         <- .rd("quatern_d")
  qx         <- .rd("qoffset_x")
  qy         <- .rd("qoffset_y")
  qz         <- .rd("qoffset_z")
  labels_arr <- .rd("labels")

  if (is.null(dims) || length(dims)<5 || dims[1]!=4) {
    stop("Invalid or missing 'dim' in /header/dim")
  }
  X <- dims[2]; Y <- dims[3]; Z <- dims[4]
  nVols <- dims[5]
  if (is.null(labels_arr) || length(labels_arr) != nVols) {
    warning("Mismatch: #labels != nVols. Possibly incomplete 'labels' dataset?")
  }

  # read /mask => 3D
  mask_dset <- file[["mask"]]
  if (is.null(mask_dset)) stop("No '/mask' dataset found.")
  mask_arr <- mask_dset$read()
  stopifnot(all(dim(mask_arr) == c(X,Y,Z)))

  # Rebuild 4x4 transform from quaternion
  if (!is.null(pixdim) && length(pixdim) >= 4) {
    qfac <- pixdim[1]
    dx   <- pixdim[2]
    dy   <- pixdim[3]
    dz   <- pixdim[4]
  } else {
    qfac <- 1; dx <- 1; dy <- 1; dz <- 1
  }

  if (!is.null(qb) && !is.null(qc) && !is.null(qd) &&
      !is.null(qx) && !is.null(qy) && !is.null(qz)) {
    mat <- neuroim2::quaternToMatrix(
      quat     = c(qb,qc,qd),
      origin   = c(qx,qy,qz),
      stepSize = c(dx,dy,dz),
      qfac     = qfac
    )
  } else {
    # fallback => identity
    mat <- diag(4)
    mat[1,1] <- dx
    mat[2,2] <- dy
    mat[3,3] <- dz
  }

  # build space => just do a 3D NeuroSpace
  spc <- NeuroSpace(dim=c(X,Y,Z), spacing=c(dx,dy,dz), trans=mat)

  # build mask
  mask_vol <- LogicalNeuroVol(as.logical(mask_arr), space=spc)

  # environment for lazy loading
  load_env <- new.env(parent=emptyenv())
  load_env$file     <- file
  load_env$data_grp <- file[["data"]]
  load_env$mask_idx <- which(mask_arr==1)
  load_env$dims     <- c(X,Y,Z)
  load_env$labels   <- labels_arr
  load_env$nVols    <- nVols
  load_env$space    <- spc

  # a loader function
  loader <- function(i) {
    lab  <- labels_arr[i]
    ds   <- load_env$data_grp[[lab]]
    val1 <- ds[]
    vol  <- array(0, dim=load_env$dims)
    vol[load_env$mask_idx] <- val1
    DenseNeuroVol(vol, space=spc)
  }

  if (memoise) {
    cache <- vector("list", nVols)
    load_env$load_vol <- function(i) {
      if (!is.null(cache[[i]])) return(cache[[i]])
      v <- loader(i)
      cache[[i]] <- v
      v
    }
  } else {
    load_env$load_vol <- loader
  }

  # Return the S4 object (LabeledVolumeSet)
  out_obj <- new("LabeledVolumeSet",
                 obj      = file,
                 mask     = mask_vol,
                 labels   = labels_arr,
                 load_env = load_env)
  out_obj
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
  signature = signature(x="LabeledVolumeSet", i="numeric", j="numeric", drop="ANY"),
  definition = function(x, i, j, k, l, ..., drop=TRUE) {

    # 1) Figure out any missing dims => use full range
    dims_3d <- dim(space(x@mask))
    nVols   <- length(x@labels)

    if (missing(k)) {
      k <- seq_len(dims_3d[3])
    }
    if (missing(l)) {
      l <- seq_len(nVols)
    }

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
    mask_arr <- as.array(x@mask)
    loader   <- x@load_env$load_vol

    out_l_pos <- 1
    for (lv in l) {
      vol_3d <- loader(lv)
      vol_arr <- as.array(vol_3d)
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

#' linear_access for LabeledVolumeSet
#'
#' @param x A \code{LabeledVolumeSet} object.
#' @param i numeric vector of 1D indices in the full [X * Y * Z * nVols] space.
#'
#' @return A numeric vector of length \code{length(i)}, in the same order as \code{i}.
#' @export
setMethod(
  f = "linear_access",
  signature = signature(x="LabeledVolumeSet", i="numeric"),
  definition = function(x, i) {

    dims_3d <- dim(space(x@mask))
    nVols   <- length(x@labels)
    bigDim  <- c(dims_3d, nVols)
    total   <- prod(bigDim)

    i <- as.integer(i)
    if (any(i < 1 | i > total)) {
      stop("Some indices out of range 1..", total)
    }

    sub_4d <- arrayInd(i, dim=bigDim)
    vol_groups <- split(seq_len(nrow(sub_4d)), sub_4d[,4])
    out <- numeric(length(i))

    loader <- x@load_env$load_vol

    for (v_str in names(vol_groups)) {
      v_idx <- as.integer(v_str)
      these_rows <- vol_groups[[v_str]]
      coords <- sub_4d[these_rows, , drop=FALSE]

      vol_3d <- loader(v_idx)
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
    cat("  ", crayon::yellow("HDF5 File"), "    : ",
        if (object@obj$is_valid) object@obj$get_filename() else "CLOSED", "\n", sep="")
    cat("  ", crayon::yellow("Data Path"), "    : /data/<label>\n", sep="")
    cat("  ", crayon::yellow("Mask Path"), "    : /mask\n", sep="")

    cat("\n")
  }
)
