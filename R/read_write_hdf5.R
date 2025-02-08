
#' Read or Write Cluster-Based Time-Series Datasets from HDF5
#'
#' @description
#' Two main functions for handling cluster-based time-series in HDF5:
#'
#' \itemize{
#'   \item \code{read_clustered_dataset()}: Reads a dataset and decides whether to return
#'         \code{H5ReducedClusteredVecSeq} or \code{H5ClusteredVecSeq} based on the
#'         HDF5 attribute \code{/scans/summary_only} and a user preference.
#'
#'   \item \code{write_clustered_dataset()}: Writes a cluster-based dataset to a new (or
#'         overwritten) HDF5 file, optionally storing only cluster-level summaries (if
#'         \code{summary_only=TRUE}) or the full voxel-level data.
#' }
#'
#' @import hdf5r
#' @importFrom neuroim2 space
#' @importFrom neuroim2 LogicalNeuroVol
#' @importFrom neuroim2 ClusteredNeuroVol
#' @importFrom neuroim2 NeuroSpace
#' @export
NULL

#' Read a cluster-based time-series dataset from HDF5
#'
#' @description
#' Reads the dataset under \code{/scans}, checking the \code{summary_only} attribute:
#' \itemize{
#'   \item If \code{summary_only=TRUE}, returns an \code{H5ReducedClusteredVecSeq}.
#'   \item Otherwise, returns an \code{H5ClusteredVecSeq} by default.
#' }
#' Additionally, if \code{prefer_reduced=TRUE}, attempts to read summary data even if
#' \code{summary_only=FALSE}, falling back to the full data if the summary data is absent.
#'
#' @param file A path (character) or an open \code{\link[hdf5r]{H5File}} in read mode.
#' @param scan_names An optional character vector of scan names to load. If \code{NULL},
#'   loads all scans found under \code{/scans}.
#' @param prefer_reduced Logical; if TRUE, tries to load reduced data even if
#'   \code{summary_only=FALSE} (falling back if unavailable).
#'
#' @return Either an \code{\link{H5ReducedClusteredVecSeq}} or an \code{\link{H5ClusteredVecSeq}}.
#'
#' @export
read_clustered_dataset <- function(file, scan_names=NULL, prefer_reduced=FALSE) {
  if (is.character(file)) {
    file <- hdf5r::H5File$new(file, mode="r")
  }

  # 1) Parse /header
  hdr_dim <- file$read("header/dim")  # c(4, X, Y, Z, nTimes, 1,1,1)
  dims3d  <- hdr_dim[2:4]
  nTimes  <- hdr_dim[5]

  pixdim  <- file$read("header/pixdim")
  spacing <- pixdim[2:4]

  # 2) Build mask
  mask_arr <- file$read("mask")
  spc      <- NeuroSpace(dim=dims3d, spacing=spacing)
  mask_vol <- LogicalNeuroVol(as.logical(mask_arr), spc)

  # 3) cluster_map, voxel_coords
  cluster_map <- file$read("cluster_map")
  tmp         <- file$read("voxel_coords")  # optional, if needed

  # 4) Global /clusters => cluster_ids, cluster_meta
  clus_grp     <- file$get("clusters")
  cluster_ids  <- clus_grp$read("cluster_ids")
  cluster_metadata <- data.frame()

  if ("cluster_meta" %in% names(clus_grp)) {
    meta_grp <- clus_grp$get("cluster_meta")
    fields   <- meta_grp$names()
    meta_list <- list()
    for (f in fields) {
      meta_list[[f]] <- meta_grp$read(f)
    }
    # if consistent length => as.data.frame
    lens <- sapply(meta_list, length)
    if (length(unique(lens)) == 1) {
      cluster_metadata <- as.data.frame(meta_list)
    } else {
      cluster_metadata <- meta_list
    }
  }

  label_map <- NULL
  if ("cluster_id" %in% names(cluster_metadata) &&
      "description" %in% names(cluster_metadata)) {
    ids  <- cluster_metadata$cluster_id
    desc <- cluster_metadata$description
    label_map <- setNames(as.list(ids), desc)
  }
  cVol <- ClusteredNeuroVol(mask_vol, cluster_map, label_map)

  # 5) /scans => attribute summary_only => bool
  scans_grp <- file$get("scans")
  sum_attr  <- scans_grp$attr_open("summary_only")$read()
  sum_attr  <- isTRUE(sum_attr)

  all_scans <- scans_grp$names()
  if (is.null(scan_names)) {
    scan_names <- all_scans
  } else {
    missing <- setdiff(scan_names, all_scans)
    if (length(missing)) {
      stop("Some requested scans not found: ", paste(missing, collapse=","))
    }
  }

  # read metadata for each scan
  out_scan_md <- list()
  for (scn in scan_names) {
    sgrp <- scans_grp$get(scn)
    mg   <- sgrp$get("metadata")
    flds <- mg$names()
    s_meta <- list()
    for (ff in flds) {
      s_meta[[ff]] <- mg$read(ff)
    }
    out_scan_md[[scn]] <- s_meta
  }

  # Decision logic
  # Case A: sum_attr=TRUE => definitely H5Reduced
  # Case B: sum_attr=FALSE => build H5Clustered or maybe reduced if prefer_reduced=TRUE
  if (sum_attr) {
    # summary_only => H5Reduced
    rseq <- new("H5ReducedClusteredVecSeq",
                obj              = file,
                scan_names       = scan_names,
                mask             = mask_vol,
                clusters         = cVol,
                scan_metadata    = out_scan_md,
                cluster_metadata = cluster_metadata,
                cluster_names    = character(length(cluster_ids)),
                cluster_ids      = as.integer(cluster_ids),
                n_time           = rep(nTimes, length(scan_names)),
                space            = spc)
    return(rseq)
  } else {
    # summary_only=FALSE => full data available
    if (prefer_reduced) {
      # check if each scan has clusters_summary/summary_data
      can_reduced <- TRUE
      for (scn in scan_names) {
        scn_path <- paste0("/scans/", scn, "/clusters_summary")
        if (!file$exists(scn_path)) {
          can_reduced <- FALSE
          break
        }
      }
      if (can_reduced) {
        # build H5Reduced
        rseq <- new("H5ReducedClusteredVecSeq",
                    obj              = file,
                    scan_names       = scan_names,
                    mask             = mask_vol,
                    clusters         = cVol,
                    scan_metadata    = out_scan_md,
                    cluster_metadata = cluster_metadata,
                    cluster_names    = character(length(cluster_ids)),
                    cluster_ids      = as.integer(cluster_ids),
                    n_time           = rep(nTimes, length(scan_names)),
                    space            = spc)
        return(rseq)
      } else {
        message("Reduced summary data not found for all scans; falling back to full data.")
      }
    }
    # fallback => build H5ClusteredVecSeq
    cseq <- new("H5ClusteredVecSeq",
                obj              = file,
                scan_names       = scan_names,
                mask             = mask_vol,
                clusters         = cVol,
                scan_metadata    = out_scan_md,
                cluster_metadata = cluster_metadata,
                space            = spc)
    return(cseq)
  }
}

#' Write a cluster-based time-series dataset to an HDF5 file
#'
#' @description
#' Creates or overwrites an HDF5 file with:
#' \itemize{
#'   \item NIfTI-like header in \code{/header}
#'   \item 3D mask in \code{/mask}
#'   \item 1D \code{cluster_map} in \code{/cluster_map}
#'   \item \code{voxel_coords} at \code{/voxel_coords}
#'   \item Global \code{/clusters} group for \code{cluster_ids} & optional metadata
#'   \item A \code{/scans} group with attribute \code{summary_only}
#' }
#'
#' If \code{summary_only=TRUE}, stores *only* cluster-level summary data (one 2D dataset
#' per scan). If \code{summary_only=FALSE}, stores full voxel-level data in \code{/clusters}.
#' (You could store both, but here we do either/or for simplicity.)
#'
#' @param file A path (character) or an open HDF5 file in write mode.
#' @param vecs If \code{summary_only=FALSE}, a list of \code{\link[neuroim2]{NeuroVec}} objects.
#'   If \code{summary_only=TRUE}, a list of numeric matrices (\eqn{nTime \times nClusters}).
#' @param scan_names A character vector naming each scan (same length as \code{vecs}).
#' @param mask A \code{\link[neuroim2]{LogicalNeuroVol}} for 3D geometry.
#' @param clusters A \code{\link[neuroim2]{ClusteredNeuroVol}} with \code{clusters@clusters} IDs.
#' @param scan_metadata A list of metadata lists, one per scan.
#' @param cluster_metadata An optional \code{data.frame} or \code{list} describing clusters.
#' @param summary_only Logical; if TRUE, store only summary data; otherwise voxel-level data.
#' @param compression \code{integer} 0..9 specifying gzip level (default=4).
#' @param chunk_size The chunk dimension for 2D writes (default=1024).
#'
#' @return An \code{H5ClusteredVecSeq} (for full data) or \code{H5ReducedClusteredVecSeq}
#'   (for summary-only).
#'
#' @export
write_clustered_dataset <- function(file,
                                    vecs,
                                    scan_names,
                                    mask,
                                    clusters,
                                    scan_metadata,
                                    cluster_metadata = NULL,
                                    summary_only    = FALSE,
                                    compression     = 4,
                                    chunk_size      = 1024)
{
  newly_opened <- FALSE
  if (is.character(file)) {
    file <- hdf5r::H5File$new(file, mode="w")
    newly_opened <- TRUE
  }

  # Basic checks
  stopifnot(length(vecs) == length(scan_names),
            length(vecs) == length(scan_metadata),
            inherits(mask, "LogicalNeuroVol"),
            inherits(clusters, "ClusteredNeuroVol"))

  # Decide # timepoints
  if (!summary_only) {
    # Using NeuroVec => check consistent time dimension
    nTimes_each <- sapply(vecs, function(v) dim(v)[4])
    if (length(unique(nTimes_each)) != 1) {
      stop("All NeuroVec objects must share the same number of timepoints.")
    }
    nTimes <- nTimes_each[1]
  } else {
    # Summaries => [nTime, nClusters]
    rowCounts <- sapply(vecs, nrow)
    if (length(unique(rowCounts)) != 1) {
      stop("All summary matrices must share the same # of rows (timepoints).")
    }
    nTimes <- rowCounts[1]
  }

  # 1) Write minimal NIfTI-like header
  hdr_grp <- file$create_group("header")
  dims_3d <- dim(mask)
  hdr_dim <- c(4L, dims_3d, nTimes, 1L, 1L, 1L)
  hdr_grp$create_dataset(
    "dim",
    dims  = 8,
    dtype = hdf5r::h5types$H5T_NATIVE_INT16
  )$write(hdr_dim)

  sp <- space(mask)
  # If user metadata for first scan has "TR", use it; else default 2.0
  pixdim <- c(
    0,
    sp@spacing[1],
    sp@spacing[2],
    sp@spacing[3],
    scan_metadata[[1]]$TR %||% 2.0,
    0, 0, 0
  )
  hdr_grp$create_dataset(
    "pixdim",
    dims  = 8,
    dtype = hdf5r::h5types$H5T_NATIVE_FLOAT
  )$write(as.numeric(pixdim))

  # 2) /mask => [X,Y,Z]
  mask_arr <- as.array(mask)
  file$create_dataset(
    "mask",
    dims  = dims_3d,
    dtype = hdf5r::h5types$H5T_NATIVE_UINT8
  )$write(as.integer(mask_arr))

  # 3) /cluster_map => length sum(mask)
  nVox <- sum(mask)
  file$create_dataset(
    "cluster_map",
    dims  = nVox,
    dtype = hdf5r::h5types$H5T_NATIVE_INT32
  )$write(as.integer(clusters@clusters))

  # 4) /voxel_coords => [nVox, 3]
  vox_coords <- which(mask, arr.ind=TRUE)
  file$create_dataset(
    "voxel_coords",
    dims  = c(nVox, 3),
    dtype = hdf5r::h5types$H5T_NATIVE_INT32
  )$write(vox_coords)

  # 5) /clusters => cluster_ids + optional metadata
  clus_grp    <- file$create_group("clusters")
  cluster_ids <- sort(unique(clusters@clusters))
  clus_grp$create_dataset(
    "cluster_ids",
    dims  = length(cluster_ids),
    dtype = hdf5r::h5types$H5T_NATIVE_INT32
  )$write(as.integer(cluster_ids))

  if (!is.null(cluster_metadata)) {
    meta_grp <- clus_grp$create_group("cluster_meta")
    if (is.data.frame(cluster_metadata)) {
      for (coln in names(cluster_metadata)) {
        dat <- cluster_metadata[[coln]]
        meta_grp$create_dataset(
          coln,
          dims  = length(dat),
          dtype = guess_h5_type(dat)
        )$write(dat)
      }
    } else if (is.list(cluster_metadata)) {
      for (fld in names(cluster_metadata)) {
        dat <- cluster_metadata[[fld]]
        meta_grp$create_dataset(
          fld,
          dims  = length(dat),
          dtype = guess_h5_type(dat)
        )$write(dat)
      }
    }
  }

  # 6) /scans => summary_only attribute
  scans_grp <- file$create_group("scans")
  scans_grp$create_attribute(
    "summary_only",
    space = hdf5r::H5S$new(dims=1),
    dtype = hdf5r::h5types$H5T_NATIVE_HBOOL
  )$write(as.logical(summary_only))

  # 7) Write each scan
  for (i in seq_along(vecs)) {
    sname   <- scan_names[i]
    scn_grp <- scans_grp$create_group(sname)

    # /scans/<scan_name>/metadata
    meta_subgrp <- scn_grp$create_group("metadata")
    smeta       <- scan_metadata[[i]]
    for (fld in names(smeta)) {
      dat <- smeta[[fld]]
      meta_subgrp$create_dataset(
        fld,
        dims  = length(dat),
        dtype = guess_h5_type(dat)
      )$write(dat)
    }

    if (!summary_only) {
      # Full voxel-level data => /clusters
      cl_subgrp <- scn_grp$create_group("clusters")
      # build [nVox, nTimes]
      mat_data <- as.matrix(vecs[[i]])[mask, ]  # => [nVox, nTimes]
      for (cid in cluster_ids) {
        idx   <- which(clusters@clusters == cid)
        cdata <- mat_data[idx, , drop=FALSE]
        dname <- paste0("cluster_", cid)
        chunk_dims <- c(min(nrow(cdata), chunk_size),
                        min(nTimes,      chunk_size))
        cl_subgrp$create_dataset(
          dname,
          dims      = c(nrow(cdata), nTimes),
          dtype     = hdf5r::h5types$H5T_NATIVE_FLOAT,
          chunk     = chunk_dims,
          gzip_level= compression
        )$write(cdata)
      }
    } else {
      # Summary-only => single dataset summary_data => shape [nTime, nClusters]
      cl_subgrp <- scn_grp$create_group("clusters_summary")
      summ_mat  <- vecs[[i]]  # user gave [nTime, nClusters]
      if (ncol(summ_mat) != length(cluster_ids)) {
        stop("Summary matrix has # of columns != # of cluster IDs!")
      }
      chunk_dims <- c(min(nrow(summ_mat), chunk_size),
                      min(ncol(summ_mat), chunk_size))
      dset <- cl_subgrp$create_dataset(
        "summary_data",
        dims       = dim(summ_mat),
        dtype      = hdf5r::h5types$H5T_NATIVE_FLOAT,
        chunk      = chunk_dims,
        gzip_level = compression
      )
      dset$write(summ_mat)
    }
  }

  # 8) Return appropriate object
  out_clust_md <- data.frame()
  if (!is.null(cluster_metadata)) {
    if (is.data.frame(cluster_metadata)) {
      out_clust_md <- cluster_metadata
    } else {
      out_clust_md <- as.data.frame(cluster_metadata)
    }
  }

  sp <- space(mask)
  if (!summary_only) {
    # Return H5ClusteredVecSeq
    obj <- new("H5ClusteredVecSeq",
               obj              = file,
               scan_names       = scan_names,
               mask             = mask,
               clusters         = clusters,
               scan_metadata    = scan_metadata,
               cluster_metadata = out_clust_md,
               space            = sp)
    return(obj)
  } else {
    # Return H5ReducedClusteredVecSeq
    obj <- new("H5ReducedClusteredVecSeq",
               obj              = file,
               scan_names       = scan_names,
               mask             = mask,
               clusters         = clusters,
               scan_metadata    = scan_metadata,
               cluster_metadata = out_clust_md,
               cluster_names    = character(length(cluster_ids)),
               cluster_ids      = as.integer(cluster_ids),
               n_time           = rep(nTimes, length(scan_names)),
               space            = sp)
    return(obj)
  }
}

#' Internal helper to guess the HDF5 type from an R vector
#'
#' @keywords internal
guess_h5_type <- function(x) {
  if (is.character(x)) {
    return(hdf5r::h5types$H5T_STRING)
  } else if (is.integer(x)) {
    return(hdf5r::h5types$H5T_NATIVE_INT32)
  } else if (is.numeric(x)) {
    return(hdf5r::h5types$H5T_NATIVE_DOUBLE)
  } else if (is.logical(x)) {
    return(hdf5r::h5types$H5T_NATIVE_HBOOL)
  }
  stop("Unsupported R data type for guess_h5_type().")
}
