
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
read_clustered_dataset <- function(file, scan_names = NULL, prefer_reduced = FALSE) {
  # 1) Open HDF5 file if 'file' is a character path
  if (is.character(file)) {
    file <- hdf5r::H5File$new(file, mode = "r")
  }

  ## 2) Read the header group data
  #    "header/dim" => a 1D dataset of length 8, e.g. c(4, X, Y, Z, nTimes, 1,1,1)
  hdr_dim <- file[["header/dim"]][]
  if (length(hdr_dim) != 8) {
    stop("Unexpected 'dim' dataset shape: should be length 8.")
  }
  # We assume 'hdr_dim' = [4, X, Y, Z, nTimes, 1, 1, 1]
  dims_3d <- hdr_dim[2:4]
  nTimes  <- hdr_dim[5]

  # "header/pixdim" => also length 8, e.g. c(0, dx, dy, dz, TR, 0,0,0)
  pixdim <- file[["header/pixdim"]][]
  if (length(pixdim) != 8) {
    stop("Unexpected 'pixdim' dataset shape: should be length 8.")
  }

  # Spacing = dx, dy, dz
  spacing <- pixdim[2:4]

  ## 3) Read the 3D mask => /mask
  # shape: (X, Y, Z)
  mask_arr <- file[["mask"]][]
  if (!all(dim(mask_arr) == dims_3d)) {
    stop("Mismatch between 'mask' shape and 'header/dim' 3D portion.")
  }
  spc <- NeuroSpace(dim = dims_3d, spacing = spacing)
  mask_vol <- LogicalNeuroVol(as.logical(mask_arr), spc)

  ## 4) cluster_map => /cluster_map
  # shape: (nVox) 1D
  cluster_map <- file[["cluster_map"]][]
  nVox <- sum(mask_vol)

  if (length(cluster_map) != nVox) {
    stop("Length of 'cluster_map' does not match the number of TRUE voxels in 'mask'.")
  }

  ## 5) voxel_coords => /voxel_coords
  # shape: (nVox, 3)
  voxel_coords <- file[["voxel_coords"]][, ]
  if (!all(dim(voxel_coords) == c(nVox, 3))) {
    stop("Shape of 'voxel_coords' is not [nVox, 3].")
  }

  ## 6) The '/clusters' group => cluster_ids + optional cluster_meta
  clus_grp <- file[["clusters"]]
  if (is.null(clus_grp)) {
    stop("File is missing the '/clusters' group.")
  }
  cluster_ids <- clus_grp[["cluster_ids"]][]
  # cluster_ids is 1D, length == number_of_active_clusters_in_the_data
  # possibly not the same as 1:10 but the actual set in the data

  # cluster_metadata => stored in /clusters/cluster_meta if present
  cluster_metadata <- data.frame()
  if ("cluster_meta" %in% names(clus_grp)) {
    meta_grp <- clus_grp[["cluster_meta"]]
    fields   <- names(meta_grp)
    meta_list <- list()
    for (f in fields) {
      meta_list[[f]] <- meta_grp[[f]][]  # read entire dataset
    }
    # If all have the same length, coerce to data frame
    lens <- sapply(meta_list, length)
    if (length(unique(lens)) == 1) {
      cluster_metadata <- as.data.frame(meta_list)
    } else {
      cluster_metadata <- meta_list
    }
  }

  # Optionally build a label_map if 'cluster_id' and 'description' columns exist
  label_map <- NULL
  if ("cluster_id" %in% names(cluster_metadata) &&
      "description" %in% names(cluster_metadata)) {
    # associate each cluster_id with a textual description
    ids  <- cluster_metadata[["cluster_id"]]
    desc <- cluster_metadata[["description"]]
    label_map <- setNames(as.list(ids), desc)
  }

  # Construct the ClusteredNeuroVol
  cVol <- ClusteredNeuroVol(mask_vol, cluster_map, label_map)

  ## 7) The '/scans' group => check 'summary_only' attribute
  scans_grp <- file[["scans"]]
  if (is.null(scans_grp)) {
    stop("File is missing the '/scans' group.")
  }
  # The hdf5r recommended way to get an attribute is via h5attr()
  sum_attr <- h5attr(scans_grp, "summary_only")
  summary_only <- isTRUE(sum_attr)

  # all scan subgroups in /scans
  all_scans <- names(scans_grp)
  # user might specify a subset
  if (is.null(scan_names)) {
    scan_names <- all_scans
  } else {
    missing_scans <- setdiff(scan_names, all_scans)
    if (length(missing_scans) > 0) {
      stop("Some requested scans not found: ", paste(missing_scans, collapse = ", "))
    }
  }

  ## 8) For each scan, read the "metadata" group
  #    We'll store that in out_scan_md[[scanName]]
  out_scan_md <- list()
  for (scn in scan_names) {
    sgrp <- scans_grp[[scn]]
    if (is.null(sgrp)) {
      stop("Scan subgroup not found: ", scn)
    }
    mg <- sgrp[["metadata"]]
    if (is.null(mg)) {
      # No metadata group => we can store an empty list
      out_scan_md[[scn]] <- list()
    } else {
      flds <- names(mg)
      s_meta <- list()
      for (ff in flds) {
        s_meta[[ff]] <- mg[[ff]][]
      }
      out_scan_md[[scn]] <- s_meta
    }
  }

  ## 9) Decide on what type of object to return
  #    If summary_only => return H5ReducedClusteredVecSeq
  #    Otherwise => H5ClusteredVecSeq
  #    But if (prefer_reduced) and summary data is present => do that.

  if (summary_only) {
    # summary_only => definitely H5Reduced
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
    # summary_only=FALSE => we wrote full voxel-level data
    # If user set prefer_reduced=TRUE, check if we actually have /clusters_summary for each scan
    if (prefer_reduced) {
      can_reduced <- TRUE
      for (scn in scan_names) {
        scn_path <- paste0("/scans/", scn, "/clusters_summary")
        if (!file$exists(scn_path)) {
          can_reduced <- FALSE
          break
        }
      }
      if (can_reduced) {
        # We can build an H5ReducedClusteredVecSeq
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
        message("Reduced summary data not found for all scans; returning full data object.")
      }
    }
    # Fallback => full data => H5ClusteredVecSeq
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
                                    chunk_size      = 1024) {
  newly_opened <- FALSE
  if (is.character(file)) {
    file <- hdf5r::H5File$new(file, mode = "w")
    newly_opened <- TRUE
  }

  # Basic checks
  stopifnot(length(vecs) == length(scan_names),
            length(vecs) == length(scan_metadata),
            inherits(mask, "LogicalNeuroVol"),
            inherits(clusters, "ClusteredNeuroVol"))

  # Decide on number of timepoints
  if (!summary_only) {
    # For full NeuroVecs, ensure all have the same 4th dimension
    nTimes_each <- sapply(vecs, function(v) dim(v)[4])
    if (length(unique(nTimes_each)) != 1) {
      stop("All NeuroVec objects must share the same number of timepoints.")
    }
    nTimes <- nTimes_each[1]
  } else {
    # For summaries: each summary matrix must have the same number of rows (timepoints)
    rowCounts <- sapply(vecs, nrow)
    if (length(unique(rowCounts)) != 1) {
      stop("All summary matrices must share the same number of rows (timepoints).")
    }
    nTimes <- rowCounts[1]
  }

  ## 1) Write a minimal NIfTI-like header in group "header"
  hdr_grp <- file$create_group("header")
  dims_3d <- dim(mask)
  hdr_dim <- c(4L, dims_3d, nTimes, 1L, 1L, 1L)
  hdr_dim_arr <- array(hdr_dim, dim = c(8))
  # Instead of using a non-existent create_attribute(), assign the dataset:
  hdr_grp[["dim"]] <- hdr_dim_arr

  sp <- space(mask)
  # Use the TR value from the first scan’s metadata (or use a default)
  pixdim <- c(0,
              sp@spacing[1],
              sp@spacing[2],
              sp@spacing[3],
              scan_metadata[[1]]$TR,
              0, 0, 0)
  hdr_grp[["pixdim"]] <- pixdim

  ## 2) Write the mask dataset (/mask)
  mask_arr <- as.array(mask)
  dset <- file$create_dataset(
    "mask",
    robj = mask_arr@.Data,
    dtype = hdf5r::h5types$H5T_NATIVE_UINT8
  )
  dset[,,] <- mask_arr

  ## 3) Write the cluster map (/cluster_map)
  nVox <- sum(mask)
  dset <- file$create_dataset(
    "cluster_map",
    dims  = nVox,
    dtype = hdf5r::h5types$H5T_NATIVE_INT32
  )
  dset[] <- clusters@clusters

  ## 4) Write voxel coordinates (/voxel_coords)
  vox_coords <- which(mask, arr.ind = TRUE)
  dset <- file$create_dataset(
    "voxel_coords",
    dims  = c(nVox, 3),
    dtype = hdf5r::h5types$H5T_NATIVE_INT32
  )
  dset[,] <- vox_coords

  ## 5) Write cluster IDs and (optionally) cluster metadata in group "/clusters"
  clus_grp <- file$create_group("clusters")
  cluster_ids <- sort(unique(clusters@clusters))
  dset <- clus_grp$create_dataset(
    "cluster_ids",
    dims  = length(cluster_ids),
    dtype = hdf5r::h5types$H5T_NATIVE_INT32
  )
  dset[] <- cluster_ids

  if (!is.null(cluster_metadata)) {
    meta_grp <- clus_grp$create_group("cluster_meta")
    if (is.data.frame(cluster_metadata)) {
      for (coln in names(cluster_metadata)) {
        dat <- cluster_metadata[[coln]]
        dset <- meta_grp$create_dataset(
          coln,
          dims  = length(dat),
          dtype = guess_h5_type(dat)
        )
        dset[] <- dat
      }
    } else if (is.list(cluster_metadata)) {
      for (fld in names(cluster_metadata)) {
        dat <- cluster_metadata[[fld]]
        dset <- meta_grp$create_dataset(
          fld,
          dims  = length(dat),
          dtype = guess_h5_type(dat)
        )
        dset[] <- dat
      }
    }
  }

  ## 6) Create the "scans" group and assign the "summary_only" attribute
  scans_grp <- file$create_group("scans")
  h5attr(scans_grp, "summary_only") <- as.logical(summary_only)

  ## 7) Write each scan’s data and metadata
  for (i in seq_along(vecs)) {
    sname <- scan_names[i]
    scn_grp <- scans_grp$create_group(sname)

    ## Write scan metadata into /scans/<scan_name>/metadata
    meta_subgrp <- scn_grp$create_group("metadata")
    smeta <- scan_metadata[[i]]
    for (fld in names(smeta)) {
      dat <- smeta[[fld]]
      dset <- meta_subgrp$create_dataset(
        fld,
        dims  = length(dat),
        dtype = guess_h5_type(dat)
      )
      dset[] <- dat
    }

    if (!summary_only) {
      ## For full voxel-level data: write into /scans/<scan_name>/clusters
      cl_subgrp <- scn_grp$create_group("clusters")
      # Extract data as a matrix [nVox, nTimes] (only for voxels in mask)
      mat_data <- as.matrix(vecs[[i]])[mask, ]
      for (cid in cluster_ids) {
        idx <- which(clusters@clusters == cid)
        cdata <- mat_data[idx, , drop = FALSE]
        dname <- paste0("cluster_", cid)
        chunk_dims <- c(min(nrow(cdata), chunk_size),
                        min(nTimes, chunk_size))
        dset <- cl_subgrp$create_dataset(
          dname,
          dims = c(nrow(cdata), nTimes),
          dtype = hdf5r::h5types$H5T_NATIVE_FLOAT,
          chunk = chunk_dims,
          gzip_level = compression
        )
        dset[] <- cdata
      }
    } else {
      ## For summary-only data: each scan is stored as a single dataset "summary_data"
      ## in group /scans/<scan_name>/clusters_summary. Its shape must be [nTime, nClusters].
      cl_subgrp <- scn_grp$create_group("clusters_summary")
      summ_mat <- vecs[[i]]  # user–provided summary matrix
      # IMPORTANT: The number of columns in summ_mat must equal the number of active clusters
      if (ncol(summ_mat) != length(cluster_ids)) {
        stop("Summary matrix has number of columns not equal to number of cluster IDs!")
      }
      chunk_dims <- c(min(nrow(summ_mat), chunk_size),
                      min(ncol(summ_mat), chunk_size))
      dset <- cl_subgrp$create_dataset(
        "summary_data",
        dims = dim(summ_mat),
        dtype = hdf5r::h5types$H5T_NATIVE_FLOAT,
        chunk = chunk_dims,
        gzip_level = compression
      )
      dset[,] <- summ_mat
    }
  }

  ## 8) Return the appropriate object
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
    obj <- new("H5ClusteredVecSeq",
               obj = file,
               scan_names = scan_names,
               mask = mask,
               clusters = clusters,
               scan_metadata = scan_metadata,
               cluster_metadata = out_clust_md,
               space = sp)
    return(obj)
  } else {
    obj <- new("H5ReducedClusteredVecSeq",
               obj = file,
               scan_names = scan_names,
               mask = mask,
               clusters = clusters,
               scan_metadata = scan_metadata,
               cluster_metadata = out_clust_md,
               cluster_names = character(length(cluster_ids)),
               cluster_ids = as.integer(cluster_ids),
               n_time = rep(nTimes, length(scan_names)),
               space = sp)
    return(obj)
  }
}

#' Internal helper to guess the HDF5 type from an R vector
#'
#' @keywords internal
guess_h5_type <- function(x) {
  if (is.character(x)) {
    return(hdf5r::H5T_STRING$new(size = Inf))
  } else if (is.integer(x)) {
    return(hdf5r::h5types$H5T_NATIVE_INT32)
  } else if (is.numeric(x)) {
    return(hdf5r::h5types$H5T_NATIVE_DOUBLE)
  } else if (is.logical(x)) {
    return(hdf5r::h5types$H5T_NATIVE_HBOOL)
  }
  stop("Unsupported R data type for guess_h5_type().")
}
