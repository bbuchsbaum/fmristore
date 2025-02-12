#' @include all_class.R
NULL

#' Create an H5ClusteredVec object
#'
#' @description
#' Constructs a single-scan \code{H5ClusteredVec} object referencing a clustered
#' time-series dataset stored in an HDF5 file.
#'
#' @param obj An \code{\link[hdf5r]{H5File}} object referencing the open HDF5 file.
#' @param scan_name \code{character}, the scan identifier (subgroup under \code{/scans/}).
#' @param mask A \code{\link[neuroim2]{LogicalNeuroVol}} representing the 3D brain mask.
#' @param clusters A \code{\link[neuroim2]{ClusteredNeuroVol}} describing cluster assignments.
#'
#' @return A new \code{H5ClusteredVec} instance.
#' @export
H5ClusteredVec <- function(obj, scan_name, mask, clusters) {
  new("H5ClusteredVec",
      obj = obj,
      scan_name = scan_name,
      mask = mask,
      clusters = clusters,
      space = space(mask))
}

#' Create an H5ClusteredVecSeq object
#'
#' @description
#' Constructs a multi-scan \code{H5ClusteredVecSeq} object referencing a set of
#' clustered time-series datasets stored in an HDF5 file.
#'
#' @param obj An \code{\link[hdf5r]{H5File}} object (open for reading).
#' @param scan_names A \code{character} vector of scan IDs to include.
#' @param mask A \code{\link[neuroim2]{LogicalNeuroVol}} giving the 3D mask.
#' @param clusters A \code{\link[neuroim2]{ClusteredNeuroVol}} describing cluster assignments.
#' @param scan_metadata A \code{list} of metadata for each scan.
#' @param cluster_metadata A \code{data.frame} of cluster metadata.
#'
#' @return A new \code{H5ClusteredVecSeq} instance.
#' @export
H5ClusteredVecSeq <- function(obj, scan_names, mask, clusters,
                              scan_metadata, cluster_metadata) {
  new("H5ClusteredVecSeq",
      obj = obj,
      scan_names = scan_names,
      mask = mask,
      clusters = clusters,
      scan_metadata = scan_metadata,
      cluster_metadata = cluster_metadata,
      space = space(mask))
}

# Extract method for H5ClusteredVec (first definition)
setMethod("[", signature(x = "H5ClusteredVec"),
          function(x, i, j, k, l, ..., drop = TRUE) {
            # Get cluster IDs for requested voxels
            clus_ids <- x@clusters@clusters[i]

            # Initialize result matrix
            result <- matrix(0, nrow=length(i), ncol=length(l))

            # Read data by cluster
            for (cid in unique(clus_ids)) {
              idx <- which(clus_ids == cid)
              dset <- paste0("/scans/", x@scan_name, "/clusters/cluster_", cid)
              clus_data <- x@obj$read(dset)[idx, l, drop=FALSE]
              result[idx,] <- clus_data
            }

            return(result)
          })

# Get a single scan from H5ClusteredVecSeq
setMethod("[[", signature(x = "H5ClusteredVecSeq"),
          function(x, i) {
            scan_name <- x@scan_names[i]
            H5ClusteredVec(obj = x@obj,
                           scan_name = scan_name,
                           mask = x@mask,
                           clusters = x@clusters)
          })

# Get total length of H5ClusteredVecSeq
setMethod("length", signature(x = "H5ClusteredVecSeq"),
          function(x) {
            sum(sapply(x@scan_names, function(name) {
              first_clus <- x@clusters@clusters[1]
              dset <- paste0("/scans/", name, "/clusters/cluster_", first_clus)
              x@obj$get(dset)$dims()[2]  # second dim is time
            }))
          })

#' Subset an H5ClusteredVec (second definition, more detailed)
#'
#' @description
#' Allows subsetting a single clustered time-series scan with standard 4D or partial
#' 4D indexing. Internally, it maps (i, j, k) to voxel indices in the mask,
#' finds which clusters those voxels belong to, then partially reads each \code{cluster_<id>}
#' dataset from the HDF5 file for the requested time dimension (l).
#'
#' @param x An \code{H5ClusteredVec} object (one scan).
#' @param i Indices for dimension 1 (x-axis), or a mask-based index if j,k are missing.
#' @param j Optional indices for dimension 2 (y-axis).
#' @param k Optional indices for dimension 3 (z-axis).
#' @param l Optional indices for time dimension. If missing, all timepoints are taken.
#' @param drop Logical (default=TRUE), whether to drop dimensions of size 1.
#' @param ... Not used.
#'
#' @return
#' If i,j,k,l are all provided, returns an array \code{[length(i), length(j), length(k), length(l)]}.
#' If some spatial dims are omitted, the result is fewer dims. If i is treated as a
#' mask-based index, you get a 2D result \code{[length(i), length(l)]}.
#'
#' @export
setMethod(
  f = "[",
  signature = signature(x="H5ClusteredVec", i="numeric", j="numeric", drop="ANY"),
  definition = function(x, i, j, k, l, ..., drop=TRUE) {
    # If user provided k, do bounding box approach for (i,j,k).
    # If user didn't provide k (and j is not missing?), or if they'd prefer 1D indexing,
    # interpret i as mask-based, j as time, etc.

    # 1) Identify dims of mask
    dims_mask <- dim(x@mask)
    missing_jk <- (missing(j) && missing(k))

    # 2) Convert i,j,k -> mask-based indices if needed
    if (!missing_j) {
      if (missing(k)) {
        k <- seq_len(dims_mask[3]) # entire z
      }
      if (max(i) > dims_mask[1] || max(j) > dims_mask[2] || (missing(k)==FALSE && max(k)>dims_mask[3])) {
        stop("Subscript out of range for this H5ClusteredVec.")
      }
      coords_3d <- as.matrix(expand.grid(x=i, y=j, z=k, KEEP.OUT.ATTRS=FALSE))
      linear_idx <- coords_3d[,1] + (coords_3d[,2]-1)*dims_mask[1] + (coords_3d[,3]-1)*dims_mask[1]*dims_mask[2]
      mask_array <- as.logical(as.array(x@mask))
      valid_idx <- linear_idx[ mask_array[linear_idx] ]
      coords_3d <- coords_3d[mask_array[linear_idx],,drop=FALSE]
    } else {
      coords_3d <- NULL
      valid_idx <- as.integer(i)
    }

    if (length(valid_idx) < 1) {
      # Return empty
      if (!missing(l)) {
        out_dim <- c(length(valid_idx), length(l))
      } else {
        out_dim <- c(length(valid_idx), 1)
      }
      return(array(0, dim=out_dim))
    }

    # 3) Time dimension
    first_clus <- x@clusters@clusters[ valid_idx[1] ]
    full_time_length <- dim(x)[4]
    if (missing(l)) {
      time_idx <- seq_len(full_time_length)
    } else {
      time_idx <- as.integer(l)
      if (max(time_idx) > full_time_length) {
        stop("Time subscript out of range for H5ClusteredVec.")
      }
    }

    # 4) Allocate result
    n_spatial <- length(valid_idx)
    n_time <- length(time_idx)
    raw_result <- matrix(NA_real_, nrow=n_spatial, ncol=n_time)

    # 5) For each cluster => partial read
    all_clus_ids <- x@clusters@clusters[ valid_idx ]
    uniq_cids <- unique(all_clus_ids)
    scans_grp <- x@obj$get(paste0("/scans/", x@scan_name, "/clusters"))

    for (cid in uniq_cids) {
      row_idx <- which(all_clus_ids == cid)
      dset_name <- paste0("cluster_", cid)
      ds <- scans_grp$get(dset_name)

      cluster_vox_all <- which(x@clusters@clusters == cid)
      row_offsets <- match(valid_idx[row_idx], cluster_vox_all)

      for (m in seq_along(row_idx)) {
        rr <- row_offsets[m]
        row_data <- ds[rr, time_idx, drop=TRUE]
        raw_result[row_idx[m], ] <- row_data
      }
    }

    # 6) Reshape if necessary
    if (!missing_j) {
      li <- length(i)
      lj <- length(j)
      lk <- if (missing(k)) 1 else length(k)
      lt <- n_time
      arr_4d <- array(raw_result, dim=c(li, lj, lk, lt))
      if (drop) {
        arr_4d <- drop(arr_4d)
      }
      return(arr_4d)
    } else {
      if (drop && (length(i) == 1 || length(time_idx) == 1)) {
        return(drop(raw_result))
      }
      return(raw_result)
    }
  }
)

#' 4D Subsetting of H5ClusteredVecSeq
#'
#' @description
#' Interprets the combined time dimension of all scans in the sequence as one
#' concatenated axis. Then extracts (i,j,k) in space plus t in time. If t spans
#' multiple scans, partial subsets are read from each \code{H5ClusteredVec} object
#' and concatenated appropriately.
#'
#' @param x An \code{H5ClusteredVecSeq} object.
#' @param i,j,k Numeric vectors specifying 3D spatial indexing (if missing, full range).
#'   If only \code{i} is provided, it may be interpreted as mask-based indices if
#'   \code{j,k} are missing.
#' @param l Numeric vector for the global time dimension across the entire sequence.
#'   E.g., if the first scan has 200 timepoints, the second 180, then times 1..200
#'   map to scan 1, times 201..380 map to scan 2, etc.
#' @param drop Logical, if \code{TRUE}, drop single-length dimensions in the result.
#' @param ... Not used.
#'
#' @return A 4D array of shape \code{[length(i), length(j), length(k), length(l)]}
#'   if fully specified, or fewer dims if partial. If only a single i or a single l
#'   is requested, dimensions may be dropped if \code{drop=TRUE}.
#'
#' @importFrom methods callNextMethod
#' @export
setMethod(
  f = "[",
  signature = signature(x="H5ClusteredVecSeq", i="numeric", j="numeric", drop="ANY"),
  definition = function(x, i, j, k, l, ..., drop=TRUE) {
    # 1) Identify total time across all scans:
    scan_names <- x@scan_names
    n_scans <- length(scan_names)
    time_lengths <- numeric(n_scans)
    for (s in seq_len(n_scans)) {
      hv <- x[[s]]
      time_lengths[s] <- dim(hv)[4]
    }
    cumsums <- c(0, cumsum(time_lengths))
    total_time <- sum(time_lengths)

    if (missing(l)) {
      l <- seq_len(total_time)
    }

    # 2) Gather valid_idx from (i,j,k) as with H5ClusteredVec
    dims_mask <- dim(x@mask)
    if (missing(k)) {
      if (!missing(j)) {
        warning("Ambiguous subsetting: j provided but k missing; using full z-range.")
        k <- seq_len(dims_mask[3])
      }
    }
    getValidMaskIndices <- function(i, j, k) {
      if (!missing(j) && !missing(k)) {
        coords_3d <- as.matrix(expand.grid(i=i, j=j, k=k))
        linear_idx <- coords_3d[,1] + (coords_3d[,2]-1)*dims_mask[1] + (coords_3d[,3]-1)*dims_mask[1]*dims_mask[2]
        mask_arr <- as.logical(as.array(x@mask))
        valid_idx <- linear_idx[mask_arr[linear_idx]]
        list(valid_idx=valid_idx, coords=coords_3d[mask_arr[linear_idx],,drop=FALSE])
      } else {
        list(valid_idx=as.integer(i), coords=NULL)
      }
    }
    sp_ix <- getValidMaskIndices(i, j, k)
    valid_idx <- sp_ix$valid_idx
    n_spatial <- length(valid_idx)
    if (n_spatial < 1) {
      return(array(0, dim=c(0,0,0,length(l))))
    }

    # 3) We'll build a [n_spatial, length(l)] matrix
    n_t <- length(l)
    raw_mat <- matrix(NA_real_, nrow=n_spatial, ncol=n_t)

    findScanForTime <- function(tval) {
      s2 <- sum(tval > cumsums)
      return(s2)
    }
    for (col_i in seq_along(l)) {
      t_global <- l[col_i]
      sc_id <- findScanForTime(t_global)
      local_time <- t_global - cumsums[sc_id]
      hv <- x[[sc_id]]
      column_data <- linear_access(hv, valid_idx)[,local_time, drop=TRUE]
      raw_mat[, col_i] <- column_data
    }

    # 4) Reshape if (i,j,k) given
    if (!missing(j) && !missing(k)) {
      li <- length(i)
      lj <- length(j)
      lk <- length(k)
      out_arr <- array(raw_mat, dim=c(li, lj, lk, n_t))
      if (drop) {
        out_arr <- drop(out_arr)
      }
      return(out_arr)
    } else {
      if (drop && (length(valid_idx)==1 || n_t==1)) {
        return(drop(raw_mat))
      }
      return(raw_mat)
    }
  }
)

#' linear_access for H5ClusteredVecSeq
#'
#' @description
#' Returns the time series for a set of mask-based voxel indices across all
#' scans in the sequence, concatenated along the time dimension.
#'
#' If the sequence has scans S1, S2, ..., each with T1, T2, ... time points,
#' the total time dimension T = T1 + T2 + ... .
#'
#' @param x An \code{H5ClusteredVecSeq} object.
#' @param i A numeric vector of mask-based voxel indices \code{(1..sum(x@mask))}.
#'
#' @return A numeric matrix of size \code{[length(i), sum_of_all_scan_timepoints]}.
#'
#' @importFrom assertthat assert_that
#' @export
setMethod(
  f = "linear_access",
  signature = signature(x="H5ClusteredVecSeq", i="numeric"),
  definition = function(x, i) {
    nVox <- sum(x@mask)
    if (any(i < 1 | i > nVox)) {
      stop("Some voxel index in 'i' is out of range [1..sum(mask)].")
    }

    scan_names <- x@scan_names
    n_scans <- length(scan_names)
    time_lengths <- numeric(n_scans)
    for (s in seq_len(n_scans)) {
      hv <- x[[s]]
      time_lengths[s] <- dim(hv)[4]
    }
    total_time <- sum(time_lengths)

    n_spat <- length(i)
    result_mat <- matrix(NA_real_, nrow=n_spat, ncol=total_time)

    cumsums <- c(0, cumsum(time_lengths))
    for (s in seq_len(n_scans)) {
      hv <- x[[s]]
      scan_nt <- time_lengths[s]
      start_col <- cumsums[s] + 1
      end_col   <- cumsums[s] + scan_nt

      # partial read => linear_access(hv, i) => shape [length(i), scan_nt]
      sub_mat <- linear_access(hv, i)
      result_mat[, start_col:end_col] <- sub_mat
    }
    result_mat
  }
)



#' Pretty Printer for H5ClusteredVec
#'
#' @description
#' Displays a concise, nicely formatted summary of an \code{H5ClusteredVec},
#' including the scan name, the cluster assignments, and the underlying
#' HDF5 file information.
#'
#' @importFrom crayon bold blue silver yellow green magenta red
#' @importFrom methods show
#' @export
setMethod(
  f = "show",
  signature = "H5ClusteredVec",
  definition = function(object) {
    # Heading
    cat("\n", crayon::bold(crayon::blue("H5ClusteredVec")), "\n", sep = "")

    # Basic line
    cat(crayon::silver("────────────────────────────────────────\n"))
    cat(crayon::bold(crayon::yellow("Basic Info")), "\n")

    # 1) Scan Name
    cat(crayon::silver(" • "), crayon::green("Scan Name:"), object@scan_name, "\n")

    # 2) Active voxels in mask
    n_vox <- sum(object@mask)
    cat(crayon::silver(" • "), crayon::green("Active voxels in mask:"), n_vox, "\n")

    # 3) Number of clusters
    cluster_ids <- unique(object@clusters@clusters)
    n_clusters  <- length(cluster_ids)
    cat(crayon::silver(" • "), crayon::green("Number of clusters:"), n_clusters, "\n")

    # 4) Attempt to read a time dimension
    #    We'll assume each cluster dataset is shaped (nVoxInCluster, nTime)
    #    so we fetch the second dimension from the first cluster (if we can).
    n_time <- NA  # if file is closed or no clusters, remains NA
    if (object@obj$is_valid) {
      if (n_clusters > 0) {
        grp_path <- file.path("/scans", object@scan_name, "clusters")
        first_cid <- cluster_ids[1]  # pick the first cluster
        ds_name   <- paste0("cluster_", first_cid)

        # Safely attempt to read dims
        ds <- object@obj$get(file.path(grp_path, ds_name))
        if (!is.null(ds)) {
          dims_ds <- ds$dims()
          # dims_ds is e.g. [nVoxInThisCluster, nTime]
          if (length(dims_ds) == 2) {
            n_time <- dims_ds[2]
          }
          ds$close()
        }
      }
    }

    # Print # time points (if known)
    if (!is.na(n_time)) {
      cat(crayon::silver(" • "), crayon::green("Time points:"), n_time, "\n")
    } else {
      cat(crayon::silver(" • "), crayon::green("Time points:"), "Unknown (file closed or no clusters)\n")
    }

    # Storage Info
    cat(crayon::bold("\nStorage:"), "\n")
    if (object@obj$is_valid) {
      cat(crayon::silver(" • "), "HDF5 file: ",
          crayon::magenta(object@obj$get_filename()), "\n", sep="")
    } else {
      cat(crayon::silver(" • "), "HDF5 file is ",
          crayon::red("CLOSED"), "\n", sep="")
    }

    cat("\n")
  }
)


#' Pretty Printer for H5ClusteredVecSeq
#'
#' @description
#' Displays a concise, nicely formatted summary of an \code{H5ClusteredVecSeq},
#' including the number of scans, the cluster assignments, and the underlying
#' HDF5 file information.
#'
#' @importFrom crayon bold blue silver yellow green magenta red
#' @importFrom methods show
#' @export
setMethod(
  f = "show",
  signature = "H5ClusteredVecSeq",
  definition = function(object) {
    # Header line
    cat("\n", crayon::bold(crayon::blue("H5ClusteredVecSeq")), "\n", sep="")

    # Basic Info section
    cat(crayon::silver("────────────────────────────────────────\n"))
    cat(crayon::bold(crayon::yellow("Basic Info")), "\n")

    # Number of scans
    n_scans <- length(object@scan_names)
    cat(crayon::silver(" • "), crayon::green("Number of scans:"), n_scans, "\n")

    # Number of active voxels from the mask
    n_vox <- sum(object@mask)
    cat(crayon::silver(" • "), crayon::green("Active voxels in mask:"), n_vox, "\n")

    # Number of clusters (based on the cluster assignment vector)
    # If object@clusters is a ClusteredNeuroVol, it has @clusters slot
    # storing the per-voxel cluster assignment:
    n_clusters <- length(unique(object@clusters@clusters))
    cat(crayon::silver(" • "), crayon::green("Number of clusters:"), n_clusters, "\n")

    # Possibly show partial scan names
    cat(crayon::bold("\nScan Names:"), "\n")
    if (n_scans <= 5) {
      cat("  ", paste(object@scan_names, collapse=", "), "\n")
    } else {
      preview <- paste(utils::head(object@scan_names, 5), collapse=", ")
      remainder <- n_scans - 5
      cat("  ", preview, crayon::silver(paste0(" ... (", remainder, " more)")), "\n")
    }

    # Show snippet of cluster metadata columns
    cat(crayon::bold("\nCluster Metadata:"), "\n")
    if (nrow(object@cluster_metadata) == 0) {
      cat(crayon::silver("  No cluster_metadata table.\n"))
    } else {
      md_cols <- names(object@cluster_metadata)
      cat(crayon::silver("  Columns:"), paste(md_cols, collapse=", "), "\n")
      cat(crayon::silver("  #rows="), nrow(object@cluster_metadata), "\n")
    }

    # Show HDF5 storage info
    cat(crayon::bold("\nStorage:"), "\n")
    if (object@obj$is_valid) {
      cat(crayon::silver(" • "), "HDF5 file: ",
          crayon::magenta(object@obj$get_filename()), "\n", sep="")
    } else {
      cat(crayon::silver(" • "), "HDF5 file is ",
          crayon::red("CLOSED"), "\n", sep="")
    }

    cat("\n")
  }
)
