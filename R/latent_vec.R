#' @importFrom neuroim2 NeuroSpace LogicalNeuroVol drop_dim space spacing origin trans axes
#' @importFrom neuroim2 quaternToMatrix matrixToQuatern write_vec
#' @importFrom fmrilatent LatentNeuroVec basis loadings offset
#' @importFrom Matrix Matrix t nnzero
#' @importFrom hdf5r H5File h5attr h5types H5T_STRING H5P_DATASET_CREATE
#' @importFrom assertthat assert_that
#' @importFrom withr defer
#' @importFrom methods setMethod signature new is as
NULL

#' HDF5 I/O for LatentNeuroVec
#'
#' @description
#' This file provides HDF5 reading and writing functionality for
#' \code{LatentNeuroVec} objects.
#'
#' The class definition and core methods are in the \code{fmrilatent} package.
#' This file provides:
#' \itemize{
#'   \item \code{write_vec} method for writing to HDF5
#'   \item \code{load_data} method for reading from HDF5
#'   \item \code{validate_latent_file} for file validation
#'   \item Internal HDF5 helper functions
#' }
#'
#' @seealso \code{LatentNeuroVec} (from \pkg{fmrilatent})
#' @name LatentNeuroVec-io
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
  new("LatentNeuroVecSource", file_name = file_name)
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
#'   \item \code{basis} matrix (as embedding)
#'   \item \code{loadings} matrix (as spatial basis)
#'   \item \code{offset} vector
#'   \item Spatial metadata
#'   \item Mask information
#' }
#' all inside an HDF5 file for future loading.
#'
#' @examples
#' \dontrun{
#' library(fmrilatent)
#' library(fmristore)
#'
#' # Create a LatentNeuroVec
#' lvec <- LatentNeuroVec(basis, loadings, space, mask)
#'
#' # Write to HDF5
#' temp_file <- tempfile(fileext = ".h5")
#' write_vec(lvec, temp_file, compression = 6)
#' }
#'
#' @seealso \code{LatentNeuroVec} (from \pkg{fmrilatent}) for the class definition.
#'
#' @importFrom neuroim2 write_vec
#' @rdname write_vec-methods
#' @export
setMethod(
  f = "write_vec",
  signature = signature(x = "LatentNeuroVec", file_name = "character", format = "missing", data_type = "missing"),
  definition = function(x, file_name, compression = 9, nbit = FALSE) {
    if (!is.character(file_name) || length(file_name) != 1) {
      stop("'file_name' must be a single character string")
    }
    if (!is.numeric(compression) || compression < 1 || compression > 9) {
      stop("'compression' must be an integer between 1 and 9")
    }

    h5obj <- to_h5_latentvec(
      vec = x,
      file_name = file_name,
      compression = compression,
      nbit = nbit
    )

    invisible(NULL)
  }
)

#' Convert LatentNeuroVec to HDF5 Format (Spec Compliant)
#' @keywords internal
#' @noRd
to_h5_latentvec <- function(vec, file_name = NULL, data_type = "FLOAT",
                            compression = 6, nbit = FALSE) {
  assert_that(inherits(vec, "LatentNeuroVec"))

  if (!is.null(file_name) && !endsWith(file_name, ".lv.h5")) {
    file_name <- paste0(file_name, ".lv.h5")
  } else if (is.null(file_name)) {
    file_name <- tempfile(fileext = ".lv.h5")
  }

  fh <- open_h5(file_name, mode = "w")
  h5obj <- fh$h5
  defer(if (fh$owns) h5obj$close_all(), envir = parent.frame())

  tryCatch(
    {
      verbose <- isTRUE(getOption("fmristore.verbose"))
      if (verbose) message("[to_h5_latentvec] Writing LatentNeuroVec to: ", file_name)

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
          list(quaternion = c(0, 0, 0), qoffset = c(0, 0, 0), qfac = 1)
        }
      )
      sp_spacing <- spacing(sp)
      if (length(sp_spacing) < 3) sp_spacing <- c(sp_spacing, rep(1, 3 - length(sp_spacing)))
      TR <- attr(sp, "TR") %||% 0.0

      h5dtype_internal <- switch(toupper(data_type),
        "FLOAT"   = hdf5r::h5types$H5T_NATIVE_FLOAT,
        "DOUBLE"  = hdf5r::h5types$H5T_NATIVE_DOUBLE,
        hdf5r::h5types$H5T_NATIVE_FLOAT
      )
      if (is.null(h5dtype_internal)) stop(paste0("Invalid data_type specified: ", data_type))

      # Use centralized map_dtype instead of inline maps
      databit <- map_dtype(h5dtype_internal)
      nifti_datatype_code <- databit[1]
      nifti_bitpix <- databit[2]
      if (nifti_datatype_code == 0L) {
        warning("Could not map data_type ", data_type, " to NIfTI codes.")
      }

      # Use shared NIfTI header builder
      hdr_fields <- build_nifti_header(
        dims = c(X, Y, Z, T_vec),
        spacing = sp_spacing,
        quat = q,
        tmat = tmat,
        datatype_code = nifti_datatype_code,
        bitpix = nifti_bitpix,
        descrip = paste("LatentNeuroVec data:", vec@label %||% "(no label)"),
        overrides = list(
          pixdim = c(q$qfac %||% 1.0, sp_spacing[1], sp_spacing[2], sp_spacing[3], TR, 0.0, 0.0, 0.0),
          slice_end = as.integer(Z - 1),
          xyzt_units = 10L
        )
      )

      .write_header(h5obj, hdr_fields, q$qfac %||% 1.0, verbose = verbose)

      nvox <- .write_mask(h5obj, vec@mask, compression, verbose = verbose)

      basis_info <- .write_basis(
        h5obj, vec@loadings, h5dtype_internal,
        compression, vec@offset, nvox, verbose = verbose
      )

      scan_name <- vec@label %||% tools::file_path_sans_ext(basename(file_name))
      .write_scans(
        h5obj, as.matrix(vec@basis), scan_name, TR,
        basis_info$k, T_vec, h5dtype_internal, compression, verbose = verbose
      )

      if (verbose) message("[to_h5_latentvec] HDF5 write SUCCESSFUL.")
      return(h5obj)
    },
    error = function(e) {
      warning(paste0("[to_h5_latentvec] ERROR during HDF5 write: ", e$message))
      return(NULL)
    }
  )
}

#' @keywords internal
#' @noRd
.write_header <- function(h5, hdr_fields, qfac, verbose = FALSE) {
  if (verbose) message("[to_h5_latentvec] Writing Header Group...")
  for (nm in names(hdr_fields)) {
    h5_write(h5, file.path("/header", nm), hdr_fields[[nm]],
      dtype = guess_h5_type(hdr_fields[[nm]]), overwrite = TRUE
    )
  }
  h5_write(h5, "/header/qfac", qfac,
    dtype = h5types$H5T_NATIVE_DOUBLE, overwrite = TRUE
  )
  if (verbose) message("[to_h5_latentvec] Header Group DONE.")
}

#' @keywords internal
#' @noRd
.write_mask <- function(h5, mask_vol, compression, verbose = FALSE) {
  if (verbose) message("[to_h5_latentvec] Writing Mask Dataset...")
  if (!inherits(mask_vol, "LogicalNeuroVol")) {
    stop("vec@mask is not a LogicalNeuroVol")
  }
  mask_arr <- array(as.integer(mask_vol@.Data), dim = dim(mask_vol))
  dims <- dim(mask_vol)
  mask_chunk_dim <- c(min(32, dims[1]), min(32, dims[2]), min(32, dims[3]))
  h5_write(h5, "/mask", mask_arr,
    dtype = hdf5r::h5types$H5T_NATIVE_UCHAR,
    chunk_dims = mask_chunk_dim, compression = compression,
    overwrite = TRUE
  )
  nvox <- sum(mask_vol)

  if (verbose) message("[to_h5_latentvec] Writing Voxel Coords...")
  mask_indices <- which(mask_arr == 1L)
  coords_to_write <- if (length(mask_indices) > 0) {
    coords <- arrayInd(mask_indices, .dim = dim(mask_arr))
    matrix(as.integer(coords - 1), nrow = length(mask_indices), ncol = 3)
  } else {
    matrix(integer(), nrow = 0, ncol = 3)
  }
  if (nrow(coords_to_write) != nvox) {
    stop("Internal dimension mismatch: voxel_coords rows != non-zero count in mask")
  }
  coord_chunk <- if (nrow(coords_to_write) > 0) c(min(1024, nrow(coords_to_write)), 3) else NULL
  h5_write(h5, "/voxel_coords", coords_to_write,
    dtype = hdf5r::h5types$H5T_NATIVE_INT32,
    chunk_dims = coord_chunk, compression = compression,
    overwrite = TRUE
  )
  return(nvox)
}

#' @keywords internal
#' @noRd
.write_basis <- function(h5, loadings, h5dtype, compression, offset, nvox_mask,
                         verbose = FALSE) {
  if (verbose) message("[to_h5_latentvec] Writing Basis Group (spatial components)...")
  h5$create_group("/basis")
  t_loadings <- Matrix::t(loadings)
  k <- nrow(t_loadings)
  nvox <- ncol(t_loadings)
  if (nvox != nvox_mask) {
    stop(paste0(
      "Internal dimension mismatch: spatial basis columns (", nvox,
      ") != nVox in mask (", nvox_mask, ")"
    ))
  }
  density <- Matrix::nnzero(t_loadings) / length(t_loadings)
  write_sparse <- density < 0.30
  if (verbose) {
    message(paste0(
      "  Spatial basis density: ", round(density * 100, 2),
      "%. Writing as ", if (write_sparse) "SPARSE" else "DENSE", "."
    ))
  }
  if (!write_sparse) {
    dense_chunk <- c(k, min(1024, nvox))
    h5_write(h5, "/basis/basis_matrix", as.matrix(t_loadings),
      dtype = h5dtype,
      chunk_dims = dense_chunk, compression = compression,
      overwrite = TRUE
    )
  } else {
    sparse_grp <- h5$create_group("/basis/basis_matrix_sparse")
    hdf5r::h5attr(sparse_grp, "storage") <- "csc"
    hdf5r::h5attr(sparse_grp, "shape") <- dim(t_loadings)
    if (!inherits(t_loadings, "dgCMatrix")) {
      t_loadings <- methods::as(t_loadings, "CsparseMatrix")
    }
    nnz <- length(t_loadings@x)
    target_min_chunk <- 128 * 1024
    element_size <- h5dtype$get_size()
    chunk_len <- max(1L, floor(target_min_chunk / element_size))
    if (nnz > 0 && chunk_len > nnz) chunk_len <- nnz else if (nnz == 0) chunk_len <- NULL
    h5_write(h5, "/basis/basis_matrix_sparse/data", t_loadings@x, h5dtype,
      chunk_dims = chunk_len, compression = compression,
      overwrite = TRUE
    )
    h5_write(h5, "/basis/basis_matrix_sparse/indices", as.integer(t_loadings@i),
      hdf5r::h5types$H5T_NATIVE_INT32,
      chunk_dims = chunk_len, compression = compression,
      overwrite = TRUE
    )
    h5_write(h5, "/basis/basis_matrix_sparse/indptr", as.integer(t_loadings@p),
      hdf5r::h5types$H5T_NATIVE_INT32,
      chunk_dims = NULL, compression = 0,
      overwrite = TRUE
    )
  }
  if (verbose) message("[to_h5_latentvec] Basis Group DONE.")
  if (verbose) message("[to_h5_latentvec] Writing Offset...")
  if (length(offset) > 0) {
    if (length(offset) != nvox) {
      stop(paste0("Offset length (", length(offset), ") does not match spatial basis nVox (", nvox, ")"))
    }
    h5_write(h5, "/offset", offset,
      dtype = h5dtype,
      chunk_dims = NULL, compression = 0, overwrite = TRUE
    )
  } else {
    if (verbose) message("  Offset is empty, skipping write.")
  }
  invisible(list(k = k, nVox = nvox))
}

#' @keywords internal
#' @noRd
.write_scans <- function(h5, embedding_matrix, scan_name, TR, k, T_vec,
                         h5dtype, compression, verbose = FALSE) {
  if (verbose) message("[to_h5_latentvec] Writing Scans Group...")
  h5$create_group("/scans")
  scan_name <- gsub("[^a-zA-Z0-9_.-]", "_", scan_name)
  if (!nzchar(scan_name)) scan_name <- "scan_1"

  meta_path <- file.path("/scans", scan_name, "metadata")
  h5$create_group(file.path("/scans", scan_name))
  h5$create_group(meta_path)
  run_length <- nrow(embedding_matrix)
  h5_write(h5, file.path(meta_path, "run_length"), as.integer(run_length), overwrite = TRUE)
  if (TR > 0) h5_write(h5, file.path(meta_path, "TR"), as.double(TR), overwrite = TRUE)
  h5_write(h5, file.path(meta_path, "subject_id"), "", dtype = hdf5r::H5T_STRING$new(size = Inf), overwrite = TRUE)
  h5_write(h5, file.path(meta_path, "task"), "", dtype = hdf5r::H5T_STRING$new(size = Inf), overwrite = TRUE)
  h5_write(h5, file.path(meta_path, "session"), "", dtype = hdf5r::H5T_STRING$new(size = Inf), overwrite = TRUE)

  if (verbose) message("  Writing embedding matrix...")
  embed_dims <- dim(embedding_matrix)
  if (embed_dims[1] != T_vec) {
    warning(paste0("Embedding time points (", embed_dims[1], ") does not match header time dim (", T_vec, "). Using embedding dim."))
  }
  if (embed_dims[2] != k) {
    stop(paste0("Embedding components (", embed_dims[2], ") mismatch spatial basis components (", k, ")"))
  }
  embed_chunk <- c(min(128, embed_dims[1]), k)
  h5_write(h5, file.path("/scans", scan_name, "embedding"), embedding_matrix,
    dtype = h5dtype, chunk_dims = embed_chunk,
    compression = compression, overwrite = TRUE
  )
  if (verbose) message("[to_h5_latentvec] Scans Group DONE.")
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
#' @importFrom fmrilatent LatentNeuroVec
#' @noRd
setMethod(
  f = "load_data",
  signature = c("LatentNeuroVecSource"),
  definition = function(x, scan_name = NULL) {
    fh <- open_h5(x@file_name, mode = "r")
    h5obj <- fh$h5
    defer(if (fh$owns) h5obj$close_all(), envir = parent.frame())

    # Check Version Attribute
    version_attr <- NULL
    if (h5obj$attr_exists("latent_spec_version")) {
      version_attr <- hdf5r::h5attr(h5obj, "latent_spec_version")
    }
    if (is.null(version_attr)) {
      warning(
        "[load_data,LatentNeuroVecSource] HDF5 file '", x@file_name,
        "' is missing the 'latent_spec_version' attribute. Assuming version 1.0 format."
      )
    } else if (version_attr != "1.0") {
      warning(
        "[load_data,LatentNeuroVecSource] HDF5 file '", x@file_name,
        "' has unexpected 'latent_spec_version'='", version_attr,
        "'. Attempting to load assuming version 1.0 format."
      )
    }

    # Read Header and Reconstruct Space
    .rd_hdr <- function(nm, required = TRUE) {
      h5_read(h5obj, file.path("/header", nm), missing_ok = !required)
    }
    dims_hdr <- .rd_hdr("dim", required = TRUE)
    pixdim_hdr <- .rd_hdr("pixdim", required = FALSE)
    qb <- .rd_hdr("quatern_b", required = FALSE)
    qc <- .rd_hdr("quatern_c", required = FALSE)
    qd <- .rd_hdr("quatern_d", required = FALSE)
    qx <- .rd_hdr("qoffset_x", required = FALSE)
    qy <- .rd_hdr("qoffset_y", required = FALSE)
    qz <- .rd_hdr("qoffset_z", required = FALSE)
    qfac <- .rd_hdr("qfac", required = FALSE)

    if (length(dims_hdr) < 5 || dims_hdr[1] != 4) {
      stop("Invalid '/header/dim' dimensions found.")
    }
    dims_3d <- dims_hdr[2:4]
    nTime_hdr <- dims_hdr[5]

    qfac_val <- if (is.null(qfac)) 1.0 else qfac
    spacing_3d <- if (!is.null(pixdim_hdr) && length(pixdim_hdr) >= 4) pixdim_hdr[2:4] else c(1, 1, 1)
    origin_3d <- c(
      if (is.null(qx)) 0 else qx,
      if (is.null(qy)) 0 else qy,
      if (is.null(qz)) 0 else qz
    )

    if (!all(sapply(list(qb, qc, qd), function(q) !is.null(q) && is.numeric(q)))) {
      warning("Missing or non-numeric quaternion b,c,d parameters. Using default orientation.")
      trans_mat <- diag(4)
      trans_mat[1, 1] <- spacing_3d[1]
      trans_mat[2, 2] <- spacing_3d[2]
      trans_mat[3, 3] <- spacing_3d[3]
      trans_mat[1:3, 4] <- origin_3d
    } else {
      trans_mat <- tryCatch(
        neuroim2::quaternToMatrix(
          quat     = c(qb, qc, qd),
          origin   = origin_3d,
          stepSize = spacing_3d,
          qfac     = qfac_val
        ),
        error = function(e) {
          warning("Error calling quaternToMatrix: ", e$message, ". Using default orientation.")
          mat_fallback <- diag(4)
          mat_fallback[1, 1] <- spacing_3d[1]
          mat_fallback[2, 2] <- spacing_3d[2]
          mat_fallback[3, 3] <- spacing_3d[3]
          mat_fallback[1:3, 4] <- origin_3d
          mat_fallback
        }
      )
    }

    full_space <- NeuroSpace(
      dim = dims_hdr[2:5],
      spacing = spacing_3d,
      origin = origin_3d,
      trans = trans_mat
    )
    space_3d <- drop_dim(full_space)

    # Read Mask
    mask_arr <- h5_read(h5obj, "/mask", missing_ok = FALSE)
    check_same_dims(dim(mask_arr), dims_3d,
      dims_to_compare = 1:3,
      msg = "load_data: Dimensions of /mask mismatch header dims"
    )

    mask_vol <- LogicalNeuroVol(as.logical(mask_arr), space = space_3d)
    nVox_mask <- sum(mask_vol)

    # Read Basis Matrix -> object@loadings
    basis_grp_exists <- h5obj$exists("/basis")
    if (!basis_grp_exists) stop("Mandatory group '/basis' not found.")
    internal_loadings <- NULL
    k_basis <- NULL
    nVox_basis <- NULL
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
      internal_loadings <- Matrix::Matrix(t(basis_matrix))
      message("[load_data] Dense spatial basis loaded.")
    } else if (has_sparse_grp) {
      sparse_grp <- h5obj[["/basis/basis_matrix_sparse"]]
      on.exit(if (!is.null(sparse_grp) && sparse_grp$is_valid) sparse_grp$close(), add = TRUE)

      storage_fmt <- "csc"
      if ("storage" %in% names(hdf5r::h5attributes(sparse_grp))) {
        storage_fmt <- tolower(hdf5r::h5attr(sparse_grp, "storage"))
      }
      shape_attr <- hdf5r::h5attr(sparse_grp, "shape")
      k_basis <- shape_attr[1]
      nVox_basis <- shape_attr[2]

      data_data <- h5_read(h5obj, "/basis/basis_matrix_sparse/data", missing_ok = FALSE)
      indices_data <- h5_read(h5obj, "/basis/basis_matrix_sparse/indices", missing_ok = FALSE)
      indptr_data <- h5_read(h5obj, "/basis/basis_matrix_sparse/indptr", missing_ok = FALSE)

      message("[load_data] Sparse triplet datasets read.")

      if (storage_fmt == "csc") {
        expected_indptr_len <- shape_attr[2] + 1
        if (length(indptr_data) != expected_indptr_len) stop("CSC indptr length mismatch.")
        reconstructed_k_p <- Matrix::sparseMatrix(
          i = indices_data, p = indptr_data, x = data_data,
          dims = shape_attr, index1 = FALSE
        )
        internal_loadings <- Matrix::t(reconstructed_k_p)
      } else if (storage_fmt == "csr") {
        expected_indptr_len <- shape_attr[1] + 1
        if (length(indptr_data) != expected_indptr_len) stop("CSR indptr length mismatch.")
        message("  Reconstructing CSR matrix (p x k) from triplet...")
        reconstructed_p_k <- Matrix::sparseMatrix(
          j = indices_data,
          p = indptr_data,
          x = data_data,
          dims = rev(shape_attr),
          index1 = FALSE
        )
        internal_loadings <- reconstructed_p_k
      } else {
        stop(paste0("Unsupported sparse storage format: ", storage_fmt))
      }

      message("[load_data] Sparse spatial basis loaded.")
    } else {
      stop("No basis matrix found in '/basis'.")
    }

    # Read Offset (Optional)
    offset_vec <- h5_read(h5obj, "/offset", missing_ok = TRUE)
    if (!is.null(offset_vec) && length(offset_vec) != nVox_basis) {
      warning("Offset length (", length(offset_vec), ") mismatch basis nVox (", nVox_basis, "). Using zero offset.")
      offset_vec <- NULL
    }
    if (is.null(offset_vec)) {
      offset_vec <- numeric(0)
    }

    # Select Scan and Read Embedding -> object@basis
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
      stop(
        "Requested scan_name '", target_scan_name, "' not found in file. Available: ",
        paste(available_scans, collapse = ", ")
      )
    }

    embedding_path <- file.path("/scans", target_scan_name, "embedding")
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

    # Construct LatentNeuroVec using constructor from fmrilatent
    LatentNeuroVec(
      basis = Matrix::Matrix(internal_basis),
      loadings = internal_loadings,
      space = full_space,
      mask = mask_vol,
      offset = offset_vec,
      label = target_scan_name
    )
  }
)


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
#'   \item Dimension consistency between header, mask, basis, and embedding.
#' }
#'
#' @importFrom hdf5r H5File h5attr h5attributes
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
  fh <- NULL

  tryCatch(
    {
      fh <- open_h5(file_path, mode = "r")
      h5obj <- fh$h5
      defer(if (!is.null(fh) && fh$owns) try(h5obj$close_all(), silent = TRUE), envir = parent.frame())

      # Check for mandatory groups/datasets
      required_groups <- c("header", "basis", "scans")
      required_datasets <- c("mask")
      required_header_datasets <- c("dim")

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
      for (ds_name in required_header_datasets) {
        hdr_ds_path <- file.path("/header", ds_name)
        if (!h5obj$exists(hdr_ds_path)) {
          stop(paste0("Mandatory dataset '", hdr_ds_path, "' not found."))
        }
      }

      # Check basis dataset/group
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
      scan_names <- tryCatch(h5obj[["scans"]]$ls()$name, finally = try(h5obj[["scans"]]$close(), silent = TRUE))
      if (length(scan_names) == 0) {
        stop("No scan subgroups found under '/scans'.")
      }
      first_scan_name <- scan_names[1]
      embed_path <- file.path("/scans", first_scan_name, "embedding")
      if (!h5obj$exists(embed_path)) {
        stop(paste0("Mandatory dataset '", embed_path, "' not found."))
      }

      # Read dimensions for consistency checks
      header_dim <- h5_read(h5obj, "/header/dim", missing_ok = FALSE)
      mask_data <- h5_read(h5obj, "/mask", missing_ok = FALSE)
      mask_dim <- dim(mask_data)

      basis_dim <- NULL
      k_basis <- NA_integer_
      nVox_basis <- NA_integer_
      if (has_dense_basis) {
        basis_dset <- NULL
        tryCatch(
          {
            basis_dset <- h5obj[[basis_dense_path]]
            basis_dim <- basis_dset$dims
          },
          finally = if (!is.null(basis_dset) && basis_dset$is_valid) try(basis_dset$close(), silent = TRUE)
        )
      } else {
        sparse_grp <- NULL
        tryCatch(
          {
            sparse_grp <- h5obj[[basis_sparse_path]]
            if (!"shape" %in% names(hdf5r::h5attributes(sparse_grp))) stop("Sparse basis group missing 'shape' attribute.")
            basis_dim <- hdf5r::h5attr(sparse_grp, "shape")
          },
          finally = if (!is.null(sparse_grp) && sparse_grp$is_valid) try(sparse_grp$close(), silent = TRUE)
        )
      }
      if (is.null(basis_dim) || length(basis_dim) != 2) stop("Failed to read valid 2D dimensions for basis.")
      k_basis <- basis_dim[1]
      nVox_basis <- basis_dim[2]

      embed_dset <- NULL
      embed_dim <- NULL
      tryCatch(
        {
          embed_dset <- h5obj[[embed_path]]
          embed_dim <- embed_dset$dims
        },
        finally = if (!is.null(embed_dset) && embed_dset$is_valid) try(embed_dset$close(), silent = TRUE)
      )
      if (is.null(embed_dim) || length(embed_dim) != 2) stop("Failed to read valid 2D dimensions for embedding.")
      T_emb <- embed_dim[1]
      k_embed <- embed_dim[2]

      # Perform Consistency Checks
      valid_checks <- list()

      if (length(header_dim) < 5 || header_dim[1] != 4) {
        valid_checks$hdr_dim0 <- paste0("'/header/dim' should start with 4, but starts with ", header_dim[1])
      }
      if (!identical(header_dim[2:4], as.integer(mask_dim))) {
        valid_checks$hdr_mask_dims <- paste0(
          "'/header/dim[2:4]' (", paste(header_dim[2:4], collapse = ","),
          ") mismatch '/mask' dims (", paste(mask_dim, collapse = ","), ")"
        )
      }
      if (k_basis != k_embed) {
        valid_checks$k_mismatch <- paste0(
          "Component dimension mismatch: basis k (", k_basis,
          ") != embedding k (", k_embed, ")"
        )
      }
      T_hdr <- header_dim[5]
      if (T_hdr != T_emb) {
        valid_checks$time_mismatch <- paste0(
          "Time dimension mismatch: header dim[5] (", T_hdr,
          ") != embedding rows (", T_emb, ")"
        )
      }
      nVox_mask <- sum(mask_data > 0)
      if (nVox_basis != nVox_mask) {
        valid_checks$basis_nvox_mask <- paste0(
          "Basis nVox mismatch: basis nVox (", nVox_basis,
          ") != non-zero count in /mask (", nVox_mask, ")"
        )
      }
      offset_val <- h5_read(h5obj, "/offset", missing_ok = TRUE)
      if (!is.null(offset_val)) {
        offset_len <- length(offset_val)
        if (offset_len != nVox_basis) {
          valid_checks$offset_len <- paste0(
            "Offset length mismatch: offset length (", offset_len,
            ") != basis nVox (", nVox_basis, ")"
          )
        }
      }

      if (length(valid_checks) > 0) {
        is_valid <- FALSE
        warning(
          "[validate_latent_file] Validation failed for '", file_path, "' with issues:\n",
          paste("  - ", unlist(valid_checks), collapse = "\n")
        )
      }
    },
    error = function(e) {
      is_valid <<- FALSE
      error_message <<- e$message
    }
  )

  if (!is.null(error_message)) {
    stop(error_message)
  }

  return(is_valid)
}

#' LatentNeuroVec constructor (re-exported from fmrilatent)
#'
#' @description
#' Creates a latent representation of neuroimaging data using basis functions
#' and loadings. This function is re-exported from the `fmrilatent` package
#' for convenience.
#'
#' @param basis A matrix of temporal embeddings (time x components).
#' @param loadings A matrix of spatial loadings (voxels x components).
#' @param space A `NeuroSpace` object defining the 4D space.
#' @param mask A `LogicalNeuroVol` specifying the brain mask.
#' @param offset Optional numeric vector of voxel offsets.
#' @param label Optional character label for the dataset.
#' @param meta Optional named list of additional metadata.
#'
#' @return A `LatentNeuroVec` object.
#'
#' @name LatentNeuroVec
#' @seealso \code{LatentNeuroVec} in package \pkg{fmrilatent} for full documentation.
#' @export
#' @importFrom fmrilatent LatentNeuroVec
LatentNeuroVec <- fmrilatent::LatentNeuroVec

#' @rdname mask-methods
#' @export
setMethod("mask", "LatentNeuroVec", function(x) x@mask)

#' Convert LatentNeuroVec to Array
#'
#' @description
#' Converts a `LatentNeuroVec` to a standard R 4D array by reconstructing
#' the full data from the latent representation.
#'
#' @param x A `LatentNeuroVec` object to convert.
#' @param ... Not used.
#'
#' @return A 4D array representing the full reconstructed data.
#'
#' @details
#' This method reconstructs the full 4D array from the latent representation,
#' which can be memory-intensive for large datasets. It calculates:
#' \deqn{array[,,, t] = Map(basis[t,] * t(loadings)) + offset}
#' for each time point, where values outside the mask are set to zero.
#'
#' @export
as.array.LatentNeuroVec <- function(x, ...) {
  dims <- dim(neuroim2::space(x))
  basis_mat <- as.matrix(x@basis)
  loadings_mat <- as.matrix(x@loadings)
  offset_vec <- x@offset
  mask_arr <- as.array(x@mask)
  mask_idx <- which(as.logical(mask_arr))

  result <- array(0, dim = dims)

  for (t in seq_len(dims[4])) {
    vals <- as.vector(tcrossprod(basis_mat[t, , drop = FALSE], loadings_mat))
    if (length(offset_vec) > 0) {
      vals <- vals + offset_vec
    }
    vol <- array(0, dim = dims[1:3])
    vol[mask_idx] <- vals
    result[, , , t] <- vol
  }

  result
}
