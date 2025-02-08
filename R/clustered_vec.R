#' @include all_class.R
NULL



# Constructor for H5ClusteredVec
H5ClusteredVec <- function(obj, scan_name, mask, clusters) {
  new("H5ClusteredVec",
      obj = obj,
      scan_name = scan_name,
      mask = mask,
      clusters = clusters,
      space = space(mask))
}

# Constructor for H5ClusteredVecSeq
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

# Extract method for H5ClusteredVec
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
      x@obj$get(dset)$dims()[2]
    }))
})



#' Subset an H5ClusteredVec
#'
#' @description
#' Allows subsetting a single clustered time-series scan with standard 4D or partial
#' 4D indexing. Internally, it maps (i, j, k) to a set of voxel indices in the mask,
#' finds which clusters those voxels belong to, then partially reads each cluster_<>
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
#' @return An array of dimension [length(i), length(j), length(k), length(l)] if
#'   all four dims are provided. If some dims are omitted, the result might have fewer dims.
#'   If all of i,j,k are collapsed into a single "mask-based" index, you'll get a 2D array
#'   [length(i), length(l)].
#'
#' @export
setMethod(
  f = "[",
  signature = signature(x="H5ClusteredVec", i="numeric", j="numeric", drop="ANY"),
  definition = function(x, i, j, k, l, ..., drop=TRUE) {
    # If user provided k, do bounding box approach for (i,j,k).
    # If user didn't provide k (and j is not missing?), or if they'd prefer 1D indexing:
    #   - We interpret i as mask-based indices and j as time, etc.
    #   For clarity, we'll handle the full 4D logic here.

    # 1) Identify dimension sizes in the mask
    dims_mask <- dim(x@mask)
    # check if user gave i,j,k or just i
    missing_jk <- (missing(j) && missing(k))

    # 2) Convert (i,j,k) -> 1D mask indices if j,k are provided
    #    or interpret i as 1D mask-based if j,k are missing
    if (!missing_j) {
      # If user gave j but not k => partial subsetting? We'll handle the logic
      if (missing(k)) {
        k <- seq_len(dims_mask[3]) # entire z
      }
      # convert i,j,k => arrayInd -> then map to mask-based index
      # Must check i,j,k are within [1..dims_mask[1]], [1..dims_mask[2]], etc.
      # ...
      # For brevity, let's do minimal bounds checks:
      if (max(i) > dims_mask[1] || max(j) > dims_mask[2] || (missing(k)==FALSE && max(k)>dims_mask[3])) {
        stop("Subscript out of range for this H5ClusteredVec.")
      }
      # Expand 3D grid
      coords_3d <- as.matrix(expand.grid(x=i, y=j, z=k, KEEP.OUT.ATTRS=FALSE))
      # Now find which of those are in the mask
      # Convert (x,y,z) to 1D => idx
      # Because x@mask is a 3D array, we do:
      linear_idx <- coords_3d[,1] + (coords_3d[,2]-1)*dims_mask[1] + (coords_3d[,3]-1)*dims_mask[1]*dims_mask[2]
      # Filter to only those that are in the mask => mask_array[linear_idx] == TRUE
      mask_array <- as.logical(as.array(x@mask))
      valid_idx <- linear_idx[ mask_array[linear_idx] ]
      # So valid_idx is a subset of row positions in the flattened volume
      # Then the final set of coords is coords_3d for those indices
      coords_3d <- coords_3d[mask_array[linear_idx],,drop=FALSE]
    } else {
      # if j,k missing => interpret i as mask-based indices
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

    # 3) Time dimension: if user gave l, use that; else all time
    # cluster data sets have shape [#voxInCluster, #timepoints]
    first_clus <- x@clusters@clusters[ valid_idx[1] ]  # just to see #time
    # read dims from one cluster dset or from space:
    # Actually we know from dim(x) => last dimension is the time length
    full_time_length <- dim(x)[4]
    if (missing(l)) {
      time_idx <- seq_len(full_time_length)
    } else {
      time_idx <- as.integer(l)
      if (max(time_idx)>full_time_length) {
        stop("Time subscript out of range for H5ClusteredVec.")
      }
    }

    # 4) We'll build a result array
    # If user gave i,j,k => shape => [length(i), length(j), length(k), length(time_idx)]
    # if they gave only i => shape => [length(i), length(time_idx)]
    # We'll store in a vector, fill in correct column-major ordering after the loop
    # For simplicity: we'll do a direct approach => a final array in R.
    # Then we rearrange it if i,j,k approach is used.

    n_spatial <- length(valid_idx)
    n_time    <- length(time_idx)
    raw_result <- matrix(NA_real_, nrow=n_spatial, ncol=n_time)

    # 5) For each cluster => partial read from "/scans/<x@scan_name>/clusters/cluster_<cid>"
    # group the valid_idx by cluster id
    all_clus_ids <- x@clusters@clusters[ valid_idx ]
    # unique clusters
    uniq_cids <- unique(all_clus_ids)

    # read partial row from each cluster for rowIndices = which(all_clus_ids== cid)
    scans_grp <- x@obj$get(paste0("/scans/", x@scan_name, "/clusters"))

    for (cid in uniq_cids) {
      row_idx <- which(all_clus_ids == cid)  # positions in valid_idx
      dset_name <- paste0("cluster_", cid)
      ds <- scans_grp$get(dset_name)

      # read partial => ds[rowIndices, time_idx]
      # but rowIndices here are the *row offsets in the cluster dataset*, i.e. 1..nVoxInCluster?
      # Problem: we do not have a direct row offset for each voxel in that cluster
      # unless we do the same approach as your existing method => the row i is "the i-th voxel in mask"
      # but that only works if you wrote cluster data in the same order as cluster_map?
      # We'll assume each cluster dataset is written in the order of mask-based indices that belong to that cluster.
      # So the row offset is the rank of each voxel among that cluster in ascending order:
      # let's create a local map for cluster cid => row index
      # we can do that once, but we might do it for each cid.
      # We'll do a small local approach:
      cluster_vox_all <- which(x@clusters@clusters == cid)
      # we want row_i for each voxel in row_idx =>
      # find match( valid_idx[row_idx], cluster_vox_all)
      # that gives partial row offsets
      row_offsets <- match(valid_idx[row_idx], cluster_vox_all)
      # Now we can do partial read for these row_offsets?
      # hdf5r doesn't allow a non-contiguous bracket subscript easily in one call.
      # We might do multiple small reads or a fancy "selection".
      # For simplicity, do multiple small reads in a loop:

      for (m in seq_along(row_idx)) {
        rr <- row_offsets[m]
        # read one row => ds[rr, time_idx]
        # or we can read all row_offsets at once if we define a sub-block. But let's keep it simple:
        row_data <- ds[rr, time_idx, drop=TRUE]
        raw_result[ row_idx[m], ] <- row_data
      }
    }

    # 6) Now raw_result is [n_spatial, n_time],
    # we rearrange it if user gave i,j,k => shape => [ length(i), length(j), length(k), length(l) ]
    if (!missing_j) {
      li <- length(i)
      lj <- length(j)
      lk <- if (missing(k)) 1 else length(k)
      lt <- n_time

      # raw_result was filled in the order valid_idx. The order of valid_idx in coords_3d is column-major if we used expand.grid(...).
      # so we can reshape accordingly:
      # We'll do an array( raw_result, c(n_spatial, n_time ) ), then we'd reorder if needed.
      # But we want final shape => c(li, lj, lk, lt).
      # We'll define an array -> dimension (li*lj*lk, lt). Then array it.
      arr_4d <- array(raw_result, dim=c(li, lj, lk, lt))
      if (drop) {
        arr_4d <- drop(arr_4d)
      }
      return(arr_4d)
    } else {
      # user only gave i => shape => [length(i), length(time_idx)]
      if (drop && (length(i) == 1 || length(time_idx)==1)) {
        return(drop(raw_result))
      }
      return(raw_result)
    }
  }
)


#' linear_access for H5ClusteredVec
#'
#' @param x An H5ClusteredVec
#' @param i A numeric vector of mask-based voxel indices (i.e. 1..sum(mask))
#' @return A numeric matrix `[length(i), nTime]` containing time series for those voxels
#'
setMethod(
  f = "linear_access",
  signature = signature(x="H5ClusteredVec", i="numeric"),
  definition = function(x, i) {
    # 1) cluster IDs for each i
    clus_ids <- x@clusters@clusters[i]

    # 2) Prepare result => [length(i), nTime]
    n_spatial <- length(i)
    # we can get #time from dim(x)[4] or read from the first cluster dataset:
    n_time <- dim(x)[4]
    result <- matrix(NA_real_, nrow=n_spatial, ncol=n_time)

    # 3) Group i by cluster => partial read from /scans/<scan_name>/clusters/cluster_<cid>
    scan_clust_grp <- x@obj$get(paste0("/scans/", x@scan_name, "/clusters"))
    for (cid in unique(clus_ids)) {
      rowidx <- which(clus_ids == cid)
      # read partial rows from cluster_<cid>
      dset_name <- paste0("cluster_", cid)
      dset <- scan_clust_grp$get(dset_name)
      # need to find row offsets => rank of each voxel in that cluster
      # e.g. full cluster membership => which(x@clusters@clusters==cid)
      clus_vox <- which(x@clusters@clusters==cid)  # all mask-based indices for that cluster
      # row_offsets => match(i[rowidx], clus_vox)
      row_offsets <- match(i[rowidx], clus_vox)

      # read row by row or attempt a selection
      # For simplicity:
      for (m in seq_along(rowidx)) {
        rr <- row_offsets[m]
        row_data <- dset[rr, , drop=TRUE]  # all time
        result[rowidx[m], ] <- row_data
      }
    }

    return(result)
  }
)


# Additional methods for H5ClusteredVec

# Get dimensions
setMethod("dim", signature(x = "H5ClusteredVec"),
  function(x) {
    first_clus <- x@clusters@clusters[1]
    dset <- paste0("/scans/", x@scan_name, "/clusters/cluster_", first_clus)
    c(dim(x@mask), x@obj$get(dset)$dims()[2])
})



# Get number of scans
setMethod("n_scans", signature(x = "H5ClusteredVecSeq"),
  function(x) {
    length(x@scan_names)
})

# Get scan names
setMethod("scan_names", signature(x = "H5ClusteredVecSeq"),
  function(x) {
    x@scan_names
})

# Get scan metadata
setMethod("scan_metadata", signature(x = "H5ClusteredVecSeq"),
  function(x) {
    x@scan_metadata
})

# Get cluster metadata
setMethod("cluster_metadata", signature(x = "H5ClusteredVecSeq"),
  function(x) {
    x@cluster_metadata
})


#' Subset an H5ClusteredVec
#'
#' @description
#' Allows subsetting a single clustered time-series scan with standard 4D or partial
#' 4D indexing. Internally, it maps (i, j, k) to a set of voxel indices in the mask,
#' finds which clusters those voxels belong to, then partially reads each cluster_<>
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
#' @return An array of dimension [length(i), length(j), length(k), length(l)] if
#'   all four dims are provided. If some dims are omitted, the result might have fewer dims.
#'   If all of i,j,k are collapsed into a single "mask-based" index, you'll get a 2D array
#'   [length(i), length(l)].
#'
#' @export
setMethod(
  f = "[",
  signature = signature(x="H5ClusteredVec", i="numeric", j="numeric", drop="ANY"),
  definition = function(x, i, j, k, l, ..., drop=TRUE) {
    # If user provided k, do bounding box approach for (i,j,k).
    # If user didn't provide k (and j is not missing?), or if they'd prefer 1D indexing:
    #   - We interpret i as mask-based indices and j as time, etc.
    #   For clarity, we'll handle the full 4D logic here.

    # 1) Identify dimension sizes in the mask
    dims_mask <- dim(x@mask)
    # check if user gave i,j,k or just i
    missing_jk <- (missing(j) && missing(k))

    # 2) Convert (i,j,k) -> 1D mask indices if j,k are provided
    #    or interpret i as 1D mask-based if j,k are missing
    if (!missing_j) {
      # If user gave j but not k => partial subsetting? We'll handle the logic
      if (missing(k)) {
        k <- seq_len(dims_mask[3]) # entire z
      }
      # convert i,j,k => arrayInd -> then map to mask-based index
      # Must check i,j,k are within [1..dims_mask[1]], [1..dims_mask[2]], etc.
      # ...
      # For brevity, let's do minimal bounds checks:
      if (max(i) > dims_mask[1] || max(j) > dims_mask[2] || (missing(k)==FALSE && max(k)>dims_mask[3])) {
        stop("Subscript out of range for this H5ClusteredVec.")
      }
      # Expand 3D grid
      coords_3d <- as.matrix(expand.grid(x=i, y=j, z=k, KEEP.OUT.ATTRS=FALSE))
      # Now find which of those are in the mask
      # Convert (x,y,z) to 1D => idx
      # Because x@mask is a 3D array, we do:
      linear_idx <- coords_3d[,1] + (coords_3d[,2]-1)*dims_mask[1] + (coords_3d[,3]-1)*dims_mask[1]*dims_mask[2]
      # Filter to only those that are in the mask => mask_array[linear_idx] == TRUE
      mask_array <- as.logical(as.array(x@mask))
      valid_idx <- linear_idx[ mask_array[linear_idx] ]
      # So valid_idx is a subset of row positions in the flattened volume
      # Then the final set of coords is coords_3d for those indices
      coords_3d <- coords_3d[mask_array[linear_idx],,drop=FALSE]
    } else {
      # if j,k missing => interpret i as mask-based indices
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

    # 3) Time dimension: if user gave l, use that; else all time
    # cluster data sets have shape [#voxInCluster, #timepoints]
    first_clus <- x@clusters@clusters[ valid_idx[1] ]  # just to see #time
    # read dims from one cluster dset or from space:
    # Actually we know from dim(x) => last dimension is the time length
    full_time_length <- dim(x)[4]
    if (missing(l)) {
      time_idx <- seq_len(full_time_length)
    } else {
      time_idx <- as.integer(l)
      if (max(time_idx)>full_time_length) {
        stop("Time subscript out of range for H5ClusteredVec.")
      }
    }

    # 4) We'll build a result array
    # If user gave i,j,k => shape => [length(i), length(j), length(k), length(time_idx)]
    # if they gave only i => shape => [length(i), length(time_idx)]
    # We'll store in a vector, fill in correct column-major ordering after the loop
    # For simplicity: we'll do a direct approach => a final array in R.
    # Then we rearrange it if i,j,k approach is used.

    n_spatial <- length(valid_idx)
    n_time    <- length(time_idx)
    raw_result <- matrix(NA_real_, nrow=n_spatial, ncol=n_time)

    # 5) For each cluster => partial read from "/scans/<x@scan_name>/clusters/cluster_<cid>"
    # group the valid_idx by cluster id
    all_clus_ids <- x@clusters@clusters[ valid_idx ]
    # unique clusters
    uniq_cids <- unique(all_clus_ids)

    # read partial row from each cluster for rowIndices = which(all_clus_ids== cid)
    scans_grp <- x@obj$get(paste0("/scans/", x@scan_name, "/clusters"))

    for (cid in uniq_cids) {
      row_idx <- which(all_clus_ids == cid)  # positions in valid_idx
      dset_name <- paste0("cluster_", cid)
      ds <- scans_grp$get(dset_name)

      # read partial => ds[rowIndices, time_idx]
      # but rowIndices here are the *row offsets in the cluster dataset*, i.e. 1..nVoxInCluster?
      # Problem: we do not have a direct row offset for each voxel in that cluster
      # unless we do the same approach as your existing method => the row i is "the i-th voxel in mask"
      # but that only works if you wrote cluster data in the same order as cluster_map?
      # We'll assume each cluster dataset is written in the order of mask-based indices that belong to that cluster.
      # So the row offset is the rank of each voxel among that cluster in ascending order:
      # let's create a local map for cluster cid => row index
      # we can do that once, but we might do it for each cid.
      # We'll do a small local approach:
      cluster_vox_all <- which(x@clusters@clusters == cid)
      # we want row_i for each voxel in row_idx =>
      # find match( valid_idx[row_idx], cluster_vox_all)
      # that gives partial row offsets
      row_offsets <- match(valid_idx[row_idx], cluster_vox_all)
      # Now we can do partial read for these row_offsets?
      # hdf5r doesn't allow a non-contiguous bracket subscript easily in one call.
      # We might do multiple small reads or a fancy "selection".
      # For simplicity, do multiple small reads in a loop:

      for (m in seq_along(row_idx)) {
        rr <- row_offsets[m]
        # read one row => ds[rr, time_idx]
        # or we can read all row_offsets at once if we define a sub-block. But let's keep it simple:
        row_data <- ds[rr, time_idx, drop=TRUE]
        raw_result[ row_idx[m], ] <- row_data
      }
    }

    # 6) Now raw_result is [n_spatial, n_time],
    # we rearrange it if user gave i,j,k => shape => [ length(i), length(j), length(k), length(l) ]
    if (!missing_j) {
      li <- length(i)
      lj <- length(j)
      lk <- if (missing(k)) 1 else length(k)
      lt <- n_time

      # raw_result was filled in the order valid_idx. The order of valid_idx in coords_3d is column-major if we used expand.grid(...).
      # so we can reshape accordingly:
      # We'll do an array( raw_result, c(n_spatial, n_time ) ), then we'd reorder if needed.
      # But we want final shape => c(li, lj, lk, lt).
      # We'll define an array -> dimension (li*lj*lk, lt). Then array it.
      arr_4d <- array(raw_result, dim=c(li, lj, lk, lt))
      if (drop) {
        arr_4d <- drop(arr_4d)
      }
      return(arr_4d)
    } else {
      # user only gave i => shape => [length(i), length(time_idx)]
      if (drop && (length(i) == 1 || length(time_idx)==1)) {
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
#' multiple scans, partial subsets are read from each relevant \code{H5ClusteredVec}
#' object and assembled.
#'
#' @param x An \code{H5ClusteredVecSeq} object.
#' @param i,j,k Numeric vectors specifying 3D spatial indexing (if missing, full range).
#'   If only \code{i} is provided, it may be interpreted as mask-based indices if
#'   \code{j,k} are missing.
#' @param l Numeric vector for the global time dimension across the entire sequence.
#'   For instance, if the first scan has 200 timepoints, the second has 180,
#'   then times 1..200 map to scan1, times 201..380 map to scan2, etc.
#' @param drop Logical, if \code{TRUE}, drop single-length dimensions in the result.
#' @param ... Not used.
#'
#' @return A 4D array of shape [length(i), length(j), length(k), length(l)] if
#'   fully specified, or fewer dims if partial. If only a single i or a single l
#'   dimension is requested, dimensions may be dropped if \code{drop=TRUE}.
#'
#' @importFrom methods callNextMethod
#' @export
setMethod(
  f = "[",
  signature = signature(x="H5ClusteredVecSeq", i="numeric", j="numeric", drop="ANY"),
  definition = function(x, i, j, k, l, ..., drop=TRUE) {
    # 1) If user didn't provide k, or partial, interpret in the same style as
    #    [ for H5ClusteredVec. We'll unify the approach for a 4D bounding box.
    #    Then the tricky part is the time dimension across multiple scans.

    # 2) Identify total time across all scans:
    # Each scan (H5ClusteredVec) has dim(...) => [X,Y,Z, nTime_s].
    # We'll store a vector of cumulative time offsets, so we know which scan
    # the requested global time indices map to.

    scan_names <- x@scan_names
    n_scans <- length(scan_names)

    # For each scan, get nTime_s:
    # We can do dim(x[[s]])[4] or read from x@clusters@...
    # But let's do a small function:
    time_lengths <- numeric(n_scans)
    for (s in seq_len(n_scans)) {
      # Use x[[s]] => returns an H5ClusteredVec
      # Then dim(...) => c(X,Y,Z,T)
      hv <- x[[s]]  # an H5ClusteredVec
      time_lengths[s] <- dim(hv)[4]
    }
    # cumulative offsets
    cumsums <- c(0, cumsum(time_lengths))  # c(0, 200, 380, ...)

    # 3) If user didn't specify l => full range
    total_time <- sum(time_lengths)
    if (missing(l)) {
      l <- seq_len(total_time)
    }

    # 4) We will build a final array. If user gave (i,j,k,l) => shape => [li, lj, lk, ll].
    # If user gave partial, we adapt. We'll do a 4D approach for clarity,
    # then drop if requested.

    # figure out how many spatial indices => we'll do the same approach as [ for H5ClusteredVec].
    # But we won't replicate all code; let's define a helper that given (i,j,k)
    # returns the list of valid mask-based indices in the 3D volume.
    # for brevity, we inline a small approach:

    dims_mask <- dim(x@mask)
    if (missing(k)) {
      # interpret i as mask-based if j is also missing
      if (!missing(j)) {
        warning("Ambiguous subsetting: j provided but k missing. Full dimension for k used.")
        k <- seq_len(dims_mask[3])
      } else if (length(i) == 1 && i[1] < 0) {
        stop("No negative or out-of-range index allowed.")
      }
    }
    # We'll unify as we did in H5ClusteredVec.
    # For brevity here, let's define a small function:
    getValidMaskIndices <- function(i, j, k) {
      if (!missing(j) && !missing(k)) {
        # bounding box approach
        coords_3d <- as.matrix(expand.grid(i=i, j=j, k=k))
        linear_idx <- coords_3d[,1] + (coords_3d[,2]-1)*dims_mask[1] + (coords_3d[,3]-1)*dims_mask[1]*dims_mask[2]
        mask_arr <- as.logical(as.array(x@mask))
        valid_idx <- linear_idx[mask_arr[linear_idx]]
        list(valid_idx=valid_idx, coords=coords_3d[mask_arr[linear_idx],,drop=FALSE])
      } else {
        # user just gave i => interpret as mask-based
        list(valid_idx=as.integer(i), coords=NULL)
      }
    }

    sp_ix <- getValidMaskIndices(i, j, k)
    valid_idx <- sp_ix$valid_idx
    n_spatial <- length(valid_idx)
    if (n_spatial < 1) {
      # return empty
      # dimension might be c(length(i), length(j), length(k), length(l)) => all zero => empty array
      return(array(0, dim=c(0,0,0,length(l))))
    }

    # 5) We'll gather results in a (n_spatial, length(l)) matrix, then reshape
    #    to [li, lj, lk, ll].
    n_t <- length(l)
    raw_mat <- matrix(NA_real_, nrow=n_spatial, ncol=n_t)

    # 6) Partition l by scan => e.g. if l = [150..250], then partial from scan1, partial from scan2, etc.
    # For each requested time index, find which scan it belongs to:
    #    if cumsums[s] < t <= cumsums[s+1], => belongs to scan s
    # We'll define a helper:
    findScanForTime <- function(tval) {
      # find s where cumsums[s] < tval <= cumsums[s+1]
      # we can do a simple binary search or which().
      s <- which(tval > cumsums & tval <= cumsums[-1] + 1, arr.ind=TRUE)
      # or simpler approach:
      s2 <- sum(tval > cumsums)
      return(s2)
    }

    # We'll do it in a loop or vectorized:
    # We might do a chunk approach for efficiency, but let's do a direct approach:
    for (col_i in seq_along(l)) {
      t_global <- l[col_i]
      sc_id <- findScanForTime(t_global) # which scan
      # local_time => t_global - cumsums[sc_id]
      local_time <- t_global - cumsums[sc_id]
      # Now partial read from x[[sc_id]] => an H5ClusteredVec
      # We'll do "linear_access" or bracket subsetting for space
      hv <- x[[sc_id]]  # an H5ClusteredVec

      # We want space= valid_idx => if user wants bounding box => we can do hv[ coords_3d, local_time] ???
      # More direct approach => linear_access with the mask-based indices 'valid_idx'
      column_data <- linear_access(hv, valid_idx)[,local_time, drop=TRUE]
      # ^ but linear_access(hv, valid_idx) => [ length(valid_idx), nTime_sc_id ]
      # We then pick local_time as a column from that
      # so:
      raw_mat[, col_i] <- column_data
    }

    # 7) Reshape if i,j,k approach => [ length(i), length(j), length(k), length(l ) ]
    # if we only have i => [ length(i), length(l) ]
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
      # just i => [ length(i), n_t ]
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
#' Returns the time series for a set of mask-based voxel indices across
#' all scans in the sequence, concatenated along the time dimension.
#'
#' Specifically, if the sequence has scans S1, S2, ..., each with T1, T2, ... time points,
#' the total time dimension is T = T1 + T2 + ... .
#' The returned matrix is [length(i), T], where T is the sum of times across scans.
#'
#' @param x An \code{H5ClusteredVecSeq} object.
#' @param i A numeric vector of mask-based voxel indices (1..sum(\code{x@mask})).
#'
#' @return A numeric matrix of size [length(i), sum_of_all_scan_timepoints].
#'
#' @importFrom assertthat assert_that
#' @export
setMethod(
  f = "linear_access",
  signature = signature(x="H5ClusteredVecSeq", i="numeric"),
  definition = function(x, i) {
    # 1) Validate indices 'i' are within [1..sum(mask)]
    nVox <- sum(x@mask)
    if (any(i < 1 | i > nVox)) {
      stop("Some voxel index in 'i' is out of range [1..sum(mask)].")
    }

    # 2) We gather total time from all scans.
    # We'll do a small helper to find the time dimension for each scan.
    scan_names <- x@scan_names
    n_scans <- length(scan_names)

    time_lengths <- numeric(n_scans)
    for (s in seq_len(n_scans)) {
      # each sub-scan is H5ClusteredVec => x[[s]]
      hv <- x[[s]]  # an H5ClusteredVec
      # 'dim(hv)' => c(X, Y, Z, T_s)
      time_lengths[s] <- dim(hv)[4]
    }
    total_time <- sum(time_lengths)

    # 3) Prepare result => [ length(i), total_time ]
    n_spat <- length(i)
    result_mat <- matrix(NA_real_, nrow=n_spat, ncol=total_time)

    # 4) Loop over scans, read partial time series from each =>
    #    linear_access(H5ClusteredVec, i). That yields [length(i), T_s].
    #    Then we place it in the columns of result_mat.

    cumsums <- c(0, cumsum(time_lengths))
    for (s in seq_len(n_scans)) {
      hv <- x[[s]]  # H5ClusteredVec
      scan_nt <- time_lengths[s]

      # columns for scan s => (cumsums[s]+1) : cumsums[s+1]
      start_col <- cumsums[s] + 1
      end_col   <- cumsums[s] + scan_nt

      # partial read => linear_access(hv, i) => shape [length(i), scan_nt]
      sub_mat <- linear_access(hv, i)

      # place sub_mat in result_mat
      result_mat[, start_col:end_col] <- sub_mat
    }

    return(result_mat)
  }
)




