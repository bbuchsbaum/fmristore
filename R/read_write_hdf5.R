
#' Read a cluster-based time-series dataset from HDF5 (Full or Summaries)
#'
#' @description
#' Checks the /scans attribute 'summary_only':
#'   * If TRUE => returns H5ReducedClusteredVecSeq
#'   * Else => returns H5ClusteredVecSeq by default
#' Additionally, if prefer_reduced=TRUE, we try to read summary data even if
#' summary_only=FALSE. If that is absent, fallback to full data.
#'
#' @param file Path or open H5File
#' @param scan_names optional vector of scans
#' @param prefer_reduced boolean: if TRUE, attempt to read /clusters_summary over /clusters
#'
#' @return H5ClusteredVecSeq or H5ReducedClusteredVecSeq
#' @export
read_clustered_dataset <- function(file, scan_names=NULL, prefer_reduced=FALSE) {
  if (is.character(file)) {
    file <- hdf5r::H5File$new(file, mode="r")
  }

  # 1) parse /header
  hdr_dim <- file$read("header/dim")
  dims3d <- hdr_dim[2:4]
  nTimes <- hdr_dim[5]

  pixdim <- file$read("header/pixdim")
  spacing <- pixdim[2:4]

  # 2) build mask
  mask_arr <- file$read("mask")
  spc <- NeuroSpace(dim=dims3d, spacing=spacing)
  mask_vol <- LogicalNeuroVol(as.logical(mask_arr), spc)

  # 3) cluster_map, voxel_coords
  cluster_map <- file$read("cluster_map")
  tmp <- file$read("voxel_coords")  # we don't necessarily need them for lazy read, but let's store if needed

  # 4) global /clusters => cluster_ids, cluster_meta
  clus_grp <- file$get("clusters")
  cluster_ids <- clus_grp$read("cluster_ids")

  cluster_metadata <- data.frame()
  if ("cluster_meta" %in% names(clus_grp)) {
    meta_grp <- clus_grp$get("cluster_meta")
    fields <- meta_grp$names()
    meta_list <- list()
    for (f in fields) {
      meta_list[[f]] <- meta_grp$read(f)
    }
    # if consistent length => as.data.frame
    lens <- sapply(meta_list, length)
    if (length(unique(lens))==1) {
      cluster_metadata <- as.data.frame(meta_list)
    } else {
      cluster_metadata <- meta_list
    }
  }

  label_map <- NULL
  if ("cluster_id" %in% names(cluster_metadata) &&
      "description" %in% names(cluster_metadata)) {
    ids <- cluster_metadata$cluster_id
    desc <- cluster_metadata$description
    label_map <- setNames(as.list(ids), desc)
  }
  cVol <- ClusteredNeuroVol(mask_vol, cluster_map, label_map)

  # 5) /scans => attribute summary_only => bool
  scans_grp <- file$get("scans")
  sum_attr <- scans_grp$attr_open("summary_only")$read()
  sum_attr <- isTRUE(sum_attr)

  # read available scan subgroups
  all_scans <- scans_grp$names()
  if (is.null(scan_names)) {
    scan_names <- all_scans
  } else {
    missing <- setdiff(scan_names, all_scans)
    if (length(missing)) {
      stop("Some requested scans not found: ", paste(missing, collapse=","))
    }
  }

  # read metadata
  out_scan_md <- list()
  for (scn in scan_names) {
    sgrp <- scans_grp$get(scn)
    mg <- sgrp$get("metadata")
    flds <- mg$names()
    s_meta <- list()
    for (ff in flds) {
      s_meta[[ff]] <- mg$read(ff)
    }
    out_scan_md[[scn]] <- s_meta
  }

  # Decision logic:
  # Case A: If sum_attr=TRUE => definitely summary only => build H5ReducedClusteredVecSeq
  # Case B: If sum_attr=FALSE => we either build H5ClusteredVecSeq or
  #          if prefer_reduced=TRUE and /clusters_summary is present for each scan, build H5ReducedClusteredVecSeq
  # If the user asked prefer_reduced=TRUE but the data doesn't exist => fallback to full.

  # check if user wants reduced data
  if (sum_attr) {
    # summary_only => definitely H5Reduced
    rseq <- new("H5ReducedClusteredVecSeq",
                obj = file,
                scan_names = scan_names,
                mask = mask_vol,
                clusters = cVol,
                scan_metadata = out_scan_md,
                cluster_metadata = cluster_metadata,
                cluster_names = character(length(cluster_ids)),
                cluster_ids = as.integer(cluster_ids),
                n_time = rep(nTimes, length(scan_names)),
                space = spc)
    return(rseq)
  } else {
    # summary_only=FALSE => we have full data
    # maybe also have summary => if user asked prefer_reduced => check presence
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
        # so we build H5ReducedClusteredVecSeq
        rseq <- new("H5ReducedClusteredVecSeq",
                    obj = file,
                    scan_names = scan_names,
                    mask = mask_vol,
                    clusters = cVol,
                    scan_metadata = out_scan_md,
                    cluster_metadata = cluster_metadata,
                    cluster_names = character(length(cluster_ids)),
                    cluster_ids = as.integer(cluster_ids),
                    n_time = rep(nTimes, length(scan_names)),
                    space = spc)
        return(rseq)
      } else {
        message("Reduced summary data not found for all scans, falling back to full data.")
      }
    }
    # else or fallback => build H5ClusteredVecSeq
    cseq <- new("H5ClusteredVecSeq",
                obj = file,
                scan_names = scan_names,
                mask = mask_vol,
                clusters = cVol,
                scan_metadata = out_scan_md,
                cluster_metadata = cluster_metadata,
                space = spc)
    return(cseq)
  }
}

#' Write a cluster-based time-series dataset to an HDF5 file
#' (Optionally summary-only or both)
#'
#' @description
#' Creates/overwrites an HDF5 file with:
#'   * NIfTI-like header at /header
#'   * 3D mask at /mask
#'   * 1D cluster_map at /cluster_map
#'   * voxel_coords at /voxel_coords
#'   * global 'clusters' group for cluster_ids/metadata
#'   * a /scans group with attribute 'summary_only'
#'
#' If \code{summary_only=TRUE}, we store only cluster-level summary data
#' (a single 2D dataset per scan). If \code{summary_only=FALSE}, we store
#' full voxel-level data in /clusters. (Optionally you can store both if you like
#' by adding code for \emph{both} blocks, but below we do either/or for clarity.)
#'
#' Returns either \code{H5ClusteredVecSeq} (for full data) or
#' \code{H5ReducedClusteredVecSeq} (for summary-only).
#'
#' @param file character path or open H5File (write mode)
#' @param vecs If \code{summary_only=FALSE}, a list of NeuroVec objects.
#'   If \code{summary_only=TRUE}, a list of numeric matrices [nTime, nClusters].
#' @param scan_names names for each scan
#' @param mask a LogicalNeuroVol for 3D geometry
#' @param clusters a ClusteredNeuroVol with \code{clusters@clusters} IDs
#' @param scan_metadata list of list metadata for each scan
#' @param cluster_metadata optional data.frame or list describing cluster IDs
#' @param summary_only boolean: store only summary data (TRUE) or voxel-level (FALSE)
#' @param compression 0-9 gzip level
#' @param chunk_size chunk dimension for 2D writes
#'
#' @return H5ClusteredVecSeq or H5ReducedClusteredVecSeq
#'
#' @export
write_clustered_dataset <- function(file,
                                    vecs,
                                    scan_names,
                                    mask,
                                    clusters,
                                    scan_metadata,
                                    cluster_metadata = NULL,
                                    summary_only = FALSE,
                                    compression = 4,
                                    chunk_size = 1024)
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

  # Decide # timepoints from either NeuroVec or matrix row dimension
  if (!summary_only) {
    # vecs => a list of NeuroVec, check same time dim
    nTimes_each <- sapply(vecs, function(v) dim(v)[4])
    if (length(unique(nTimes_each)) != 1) {
      stop("All NeuroVec objects must share same # of timepoints.")
    }
    nTimes <- nTimes_each[1]
  } else {
    # summary => each vecs[[i]] is [nTime, nClusters]
    rowCounts <- sapply(vecs, nrow)
    if (length(unique(rowCounts)) != 1) {
      stop("All summary matrices must share same # of rows (timepoints).")
    }
    nTimes <- rowCounts[1]
  }

  # 1) Write minimal NIfTI-like header
  hdr_grp <- file$create_group("header")
  dims_3d <- dim(mask)
  hdr_dim <- c(4L, dims_3d, nTimes, 1L, 1L, 1L)
  hdr_grp$create_dataset("dim", dims=8,
                         dtype=hdf5r::h5types$H5T_NATIVE_INT16)$write(hdr_dim)

  sp <- space(mask)
  pixdim <- c(0, sp@spacing[1], sp@spacing[2], sp@spacing[3],
              scan_metadata[[1]]$TR %||% 2.0,
              0,0,0)
  hdr_grp$create_dataset("pixdim", dims=8,
                         dtype=hdf5r::h5types$H5T_NATIVE_FLOAT)$write(as.numeric(pixdim))

  # 2) /mask
  dims_3d <- dim(mask)
  file$create_dataset("mask", dims=dims_3d,
                      dtype=hdf5r::h5types$H5T_NATIVE_UINT8)$write(as.integer(as.array(mask)))

  # 3) /cluster_map => length sum(mask)
  nVox <- sum(mask)
  file$create_dataset("cluster_map",
                      dims=nVox,
                      dtype=hdf5r::h5types$H5T_NATIVE_INT32)$write(as.integer(clusters@clusters))

  # 4) /voxel_coords => shape [nVox, 3]
  vox_coords <- which(mask, arr.ind=TRUE)
  file$create_dataset("voxel_coords",
                      dims=c(nVox,3),
                      dtype=hdf5r::h5types$H5T_NATIVE_INT32)$write(vox_coords)

  # 5) /clusters group => cluster_ids + optional metadata
  clus_grp <- file$create_group("clusters")
  cluster_ids <- sort(unique(clusters@clusters))
  clus_grp$create_dataset("cluster_ids",
                          dims=length(cluster_ids),
                          dtype=hdf5r::h5types$H5T_NATIVE_INT32)$write(as.integer(cluster_ids))

  if (!is.null(cluster_metadata)) {
    meta_grp <- clus_grp$create_group("cluster_meta")
    if (is.data.frame(cluster_metadata)) {
      for (coln in names(cluster_metadata)) {
        dat <- cluster_metadata[[coln]]
        meta_grp$create_dataset(coln,
                                dims=length(dat),
                                dtype=guess_h5_type(dat))$write(dat)
      }
    } else if (is.list(cluster_metadata)) {
      for (fld in names(cluster_metadata)) {
        dat <- cluster_metadata[[fld]]
        meta_grp$create_dataset(fld,
                                dims=length(dat),
                                dtype=guess_h5_type(dat))$write(dat)
      }
    }
  }

  # 6) /scans => summary_only attribute
  scans_grp <- file$create_group("scans")
  scans_grp$create_attribute("summary_only",
                             space=hdf5r::H5S$new(dims=1),
                             dtype=hdf5r::h5types$H5T_NATIVE_HBOOL)$write(as.logical(summary_only))

  # 7) Write each scan
  for (i in seq_along(vecs)) {
    sname <- scan_names[i]
    scn_grp <- scans_grp$create_group(sname)

    # /scans/<scan_name>/metadata
    meta_subgrp <- scn_grp$create_group("metadata")
    smeta <- scan_metadata[[i]]
    for (fld in names(smeta)) {
      dat <- smeta[[fld]]
      meta_subgrp$create_dataset(fld,
                                 dims=length(dat),
                                 dtype=guess_h5_type(dat))$write(dat)
    }

    if (!summary_only) {
      # Full voxel-level data => /clusters
      cl_subgrp <- scn_grp$create_group("clusters")
      # build [nVox, nTimes]
      mat_data <- as.matrix(vecs[[i]])[mask, ] # => [nVox, nTimes]
      for (cid in cluster_ids) {
        idx <- which(clusters@clusters == cid)
        cdata <- mat_data[idx, , drop=FALSE]
        dset_name <- paste0("cluster_", cid)
        chunk_dims <- c(min(nrow(cdata), chunk_size),
                        min(nTimes, chunk_size))
        cl_subgrp$create_dataset(dset_name,
                                 dims=c(nrow(cdata), nTimes),
                                 dtype=hdf5r::h5types$H5T_NATIVE_FLOAT,
                                 chunk=chunk_dims,
                                 gzip_level=compression)$write(cdata)
      }
    } else {
      # Summary-only => a single dataset "summary_data" => shape [nTime, nClusters]
      cl_subgrp <- scn_grp$create_group("clusters_summary")
      summ_mat <- vecs[[i]]  # user gave [nTime, nClusters]
      if (ncol(summ_mat) != length(cluster_ids)) {
        stop("Summary matrix has # columns != # cluster IDs!")
      }
      chunk_dims <- c(min(nrow(summ_mat), chunk_size),
                      min(ncol(summ_mat), chunk_size))
      dset <- cl_subgrp$create_dataset("summary_data",
                                       dims=dim(summ_mat),
                                       dtype=hdf5r::h5types$H5T_NATIVE_FLOAT,
                                       chunk=chunk_dims,
                                       gzip_level=compression)
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

  if (!summary_only) {
    # Return H5ClusteredVecSeq
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
    # Return H5ReducedClusteredVecSeq
    obj <- new("H5ReducedClusteredVecSeq",
               obj = file,
               scan_names = scan_names,
               mask = mask,
               clusters = clusters,
               scan_metadata = scan_metadata,
               cluster_metadata = out_clust_md,
               # Possibly fill cluster_names, cluster_ids, n_time, etc. if known
               cluster_names = character(length(cluster_ids)), # or real names
               cluster_ids = as.integer(cluster_ids),
               n_time = rep(nTimes, length(scan_names)),
               space = sp)
    return(obj)
  }
}

# Helper guess_h5_type
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
