# test_fmristore_extras.R

library(testthat)
library(neuroim2)
#library(fmristore)  # or whatever your package is called
library(Matrix)

test_that("H5NeuroVol handles empty or single-slice subsetting", {
  # Prepare a small 3D array
  arr <- array(rnorm(2*3*4), dim=c(2,3,4))
  spc <- NeuroSpace(dim=c(2,3,4))
  vol <- NeuroVol(arr, spc)

  # Write to HDF5
  tmpfile <- tempfile(fileext=".h5")
  on.exit(unlink(tmpfile))
  h5vol <- to_nih5_vol(vol, tmpfile, data_type="DOUBLE", chunk_dim=c(2,3,1))

  # Single slice
  slice1 <- h5vol[1,,]
  expect_equal(dim(slice1), c(3,4))
  expect_equal(slice1, arr[1,,])

  # Subset that is empty
  empty_sub <- h5vol[ integer(0), 2, ]
  expect_equal(length(empty_sub), 0)
  expect_equal(dim(empty_sub), c(0,4))

  # Single voxel
  single_vox <- h5vol[1,2,3]
  expect_length(single_vox, 1)
  expect_equal(single_vox, arr[1,2,3])

  # Make sure file is still valid
  expect_true(h5vol@h5obj$is_valid)
})

test_that("H5NeuroVol error handling for out-of-range indices", {
  dims <- c(4,4,4)
  spc  <- NeuroSpace(dims)
  arr  <- array(runif(prod(dims)), dim=dims)
  vol  <- NeuroVol(arr, spc)

  tmp <- tempfile(fileext=".h5")
  on.exit(unlink(tmp))
  h5vol <- to_nih5_vol(vol, tmp, data_type="FLOAT")

  # Out-of-range
  expect_error(h5vol[5,1,1], regexp="Subscript")
  expect_error(h5vol[1:4,1:4,6], regexp="Subscript")
  # Negative index
  expect_error(h5vol[-1, , ], regexp="Subscript")

  # linear_access out of range
  expect_error(linear_access(h5vol, 0), regexp="out of range")
  expect_error(linear_access(h5vol, 200), regexp="out of range")
})

test_that("H5NeuroVec partial dimension subsetting and zero-size slices", {
  # 4D data: 3 x 3 x 2 x 3
  arr <- array(rnorm(3*3*2*3), dim=c(3,3,2,3))
  spc <- NeuroSpace(dim=c(3,3,2,3))
  vec <- NeuroVec(arr, spc)

  tmp <- tempfile(fileext=".h5")
  on.exit(unlink(tmp))

  # Force no chunking to test unchunked usage
  h5vec <- to_nih5_vec(vec, file_name=tmp, chunk_dim=c(3,3,2,3), compression=0)
  expect_true(inherits(h5vec, "H5NeuroVec"))

  # Zero-slice in one dimension
  esub <- h5vec[ , , , integer(0)]
  expect_equal(dim(esub), c(3,3,2,0))

  # Single time slice
  sub1 <- h5vec[1:2, , 1, 1]
  expect_equal(dim(sub1), c(2,3,1,1))
  expect_equal(sub1, arr[1:2, , 1, 1,drop=FALSE], tolerance=1e-6)

  # i dimension only
  only_i <- h5vec[1, , , ,drop=FALSE]
  expect_equal(dim(only_i), c(1,3,2,3))

  # negative or zero index => error
  expect_error(h5vec[-1,1,1,1], "out of range")
  expect_error(h5vec[0,1,1,1],   "out of range")
})

test_that("H5NeuroVec linear access corner cases", {
  arr <- array(rnorm(2*2*2*2), dim=c(2,2,2,2))
  spc <- NeuroSpace(dim=c(2,2,2,2))
  vec <- NeuroVec(arr, spc)

  tmp <- tempfile(fileext=".h5")
  on.exit(unlink(tmp))

  h5vec <- to_nih5_vec(vec, file_name=tmp)

  # linear_access with all possible indices
  tot <- prod(dim(h5vec))
  lv  <- linear_access(h5vec, 1:tot)
  expect_equal(lv, as.vector(arr), tolerance=1e-7)

  # random sample
  idx_samp <- c(1, 2, tot-1, tot)
  lv_samp  <- linear_access(h5vec, idx_samp)
  expect_equal(lv_samp, arr[idx_samp], tolerance=1e-7)
})

test_that("LatentNeuroVec partial subsetting and out-of-mask voxels", {
  # Suppose a 2x2x2 volume with 3 timepoints => shape c(2,2,2,3)
  n_time <- 3
  n_basis <- 2
  # basis => [3 x 2]
  basis <- Matrix(matrix(rnorm(n_time * n_basis), nrow=n_time, ncol=n_basis))

  # mask => 2x2x2 but let's mask only half
  mask_arr <- array(c(TRUE, FALSE,
                      TRUE, FALSE,
                      FALSE, TRUE,
                      FALSE, TRUE), dim=c(2,2,2))
  mask_vol <- LogicalNeuroVol(mask_arr, NeuroSpace(c(2,2,2)))
  # => sum(mask) = 4

  # loadings => [4 x 2]
  loadings <- matrix(rnorm(4*2), nrow=4, ncol=2)
  offset   <- rnorm(4)

  # Build LatentNeuroVec => space => 2,2,2,3
  spc <- NeuroSpace(dim=c(2,2,2,n_time))
  lat <- LatentNeuroVec(basis, loadings, spc, mask=mask_vol, offset=offset)

  # Subset => i=1:2, j=1:2, k=1:2 => entire space
  # We test that the "masked-out" voxels are zero
  full_4d <- lat[1:2, 1:2, 1:2, 1:3]
  # dimension => c(2,2,2,3)
  # For each time slice, out-of-mask is zero
  for (t in 1:3) {
    # Recompute the known partialVol = basis[t, ] %*% t(loadings) + offset
    slice <- as.numeric(basis[t,,drop=FALSE] %*% t(loadings) + offset)
    # Then reorder them into a 2x2x2 but masked voxels appear => rest 0
    arrt <- array(0, c(2,2,2))
    # mask row => which(mask_arr) => e.g. c(1,3,5,7) in linear indexing
    mv <- which(mask_arr)
    arrt[mv] <- slice
    # compare to full_4d[..,t]
    expect_equal(arrt, full_4d[,,,t], tolerance=1e-7)
  }

  idx_1d <- 1 + (2-1)*2 + (1-1)*2*2  # = 3
  s1 <- series(lat, 1, 2, 1)
  expect_length(s1, 3)

  rowid <- match(idx_1d, which(mask_arr))  # rowid in loadings; should now be 2
  manval <- numeric(3)
  for (t in seq_len(3)) {
    manval[t] <- sum(basis[t, ] * loadings[rowid, ]) + offset[rowid]
  }
  expect_equal(s1, manval, tolerance = 1e-8)


  # A voxel that is out of mask => e.g. i=2, j=1, k=1 => check mask
  idx_1d2 <- 2 + (1-1)*2 + (1-1)*2*2  # =2 => mask[2]=FALSE
  s2 <- series(lat, 2,1,1)
  expect_equal(s2, c(0,0,0))  # entire time series zero
})

test_that("LatentNeuroVec error handling", {
  # dimension mismatch in basis vs space
  expect_error({
    basis <- Matrix(rnorm(20), nrow=5)
    loadings <- Matrix(rnorm(50), nrow=10) # mismatch in #cols
    spc <- NeuroSpace(c(2,2,2,5))
    mask_arr <- array(TRUE, c(2,2,2))
    mask_vol <- LogicalNeuroVol(mask_arr, NeuroSpace(c(2,2,2)))
    LatentNeuroVec(basis, loadings, spc, mask_vol)
  }, regexp="must have")

  # offset length mismatch
  expect_error({
    basis <- Matrix(rnorm(10), nrow=5, ncol=2)
    loadings <- Matrix(rnorm(2*4), nrow=4, ncol=2)
    spc <- NeuroSpace(c(2,2,1,5))
    mask_arr <- array(c(TRUE,TRUE,TRUE,TRUE), dim=c(2,2,1))
    mask_vol <- LogicalNeuroVol(mask_arr, NeuroSpace(c(2,2,1)))
    LatentNeuroVec(basis, loadings, spc, mask_vol, offset=rnorm(10))
  }, "must match")
})

test_that("LabeledVolumeSet partial usage + memoise=TRUE", {
  # create a small 4D data: (4,4,2,3) => 3 volumes
  arr <- array(rnorm(4*4*2*3), dim=c(4,4,2,3))
  spc <- NeuroSpace(dim=c(4,4,2,3))
  vec <- NeuroVec(arr, spc)

  # mask => pick random subset
  mask_arr <- array(sample(c(TRUE,FALSE), 4*4*2, replace=TRUE), dim=c(4,4,2))
  mask_vol <- LogicalNeuroVol(mask_arr, NeuroSpace(c(4,4,2)))

  labels <- c("volA", "volB", "volC")

  tmpfile <- tempfile(fileext=".h5")
  on.exit(unlink(tmpfile))

  # write as labeled
  write_labeled_vec(vec=vec, mask=mask_vol, labels=labels, file=tmpfile, compression=1)

  # read it with memoise=TRUE
  lvs <- read_labeled_vec(tmpfile, memoise=TRUE)
  expect_s4_class(lvs, "LabeledVolumeSet")

  # check subsetting => [i,j,k,l]
  sub_4d <- lvs[1:2,2:3,1:2, 1:2]
  expect_equal(dim(sub_4d), c(2,2,2,2))
  # compare manually
  gold <- arr[1:2,2:3,1:2,1:2] * rep(mask_arr[1:2,2:3,1:2], 2)
  expect_equal(sub_4d, gold, tolerance=1e-7)

  # check linear_access
  lin_idx <- c(1,4,10,16,20,30)
  lv_vals <- linear_access(lvs, lin_idx)

  # construct "expected_full"
  # each vol => arr[,,,vol] * mask
  arr_masked <- arr
  for (v in seq_len(3)) {
    arr_masked[,,,v] <- arr_masked[,,,v] * mask_arr
  }
  gold_lin <- arr_masked[ lin_idx ]
  expect_equal(lv_vals, gold_lin, tolerance=1e-7)

  # check single volume => [[
  vC <- lvs[[3]]
  expect_s4_class(vC, "DenseNeuroVol")
  # verify it matches arr[,,,3]*mask
  expect_equal(as.array(vC@.Data), arr[,,,3]*mask_arr, tolerance=1e-7)
})

test_that("Clustered dataset read/write with summary_only=TRUE", {
  # Create a simple 3D mask
  dims_3d <- c(4, 4, 4)
  mask_arr <- array(FALSE, dims_3d)
  set.seed(123)
  mask_arr[sample(length(mask_arr), length(mask_arr)/5)] <- TRUE
  mvol <- LogicalNeuroVol(mask_arr, NeuroSpace(dims_3d))

  # Create a cluster map:
  # Let's say we intended to sample from 1:10, but only the active voxels count.
  cluster_ids <- 1:10
  cl_map <- integer(sum(mask_arr))
  for (j in seq_along(cl_map)) {
    cl_map[j] <- sample(cluster_ids, 1)
  }
  cVol <- ClusteredNeuroVol(mvol, cl_map)

  # Print cVol to see how many clusters are active (e.g., 7)
  print(cVol)

  # For summary data: each scan is a [nTime, nClusters_actual] matrix.
  nTime <- 5
  actual_clusters <- sort(unique(cVol@clusters))
  nClusters <- length(actual_clusters)  # This is the number of active clusters

  # Now build two summary matrices with nClusters columns (not 10)
  sc1 <- matrix(rnorm(nTime * nClusters), nrow = nTime, ncol = nClusters)
  sc2 <- matrix(rnorm(nTime * nClusters), nrow = nTime, ncol = nClusters)
  scans_list <- list(sc1, sc2)

  # Metadata
  sc_md <- list(
    list(TR = 2.0, subject = "subjA"),
    list(TR = 2.0, subject = "subjA")
  )
  cl_md <- data.frame(cluster_id = actual_clusters,
                      description = paste("Cluster", actual_clusters))

  # Write the clustered dataset
  tmpf <- tempfile(fileext = ".h5")
  on.exit(unlink(tmpf))

  wobj <- write_clustered_dataset(
    file             = tmpf,
    vecs            = scans_list,
    scan_names      = c("runA", "runB"),
    mask            = mvol,
    clusters        = cVol,
    scan_metadata   = sc_md,
    cluster_metadata = cl_md,
    summary_only    = TRUE,  # store only cluster-level summary
    compression     = 4,
    chunk_size      = 1024
  )
  expect_s4_class(wobj, "H5ReducedClusteredVecSeq")

  # read
  rseq <- read_clustered_dataset(tmpf)
  expect_s4_class(rseq, "H5ReducedClusteredVecSeq")

  # check as.matrix => list or bind
  df_list <- as.matrix(rseq, bind_scans=FALSE)
  expect_length(df_list, 2)
  expect_equal(df_list[[1]], sc1, tolerance=1e-7)
  expect_equal(df_list[[2]], sc2, tolerance=1e-7)
})

test_that("Clustered dataset full-voxel read/write", {
  # make a small 3D mask, cluster_map
  dims_3d <- c(3,3,2)
  mask_arr <- array(TRUE, dims_3d)
  mask <- LogicalNeuroVol(mask_arr, NeuroSpace(dims_3d))

  # define clusters => each voxel belongs to 1..4
  cMap <- integer(prod(dims_3d))
  cMap[] <- sample(1:4, length(cMap), replace=TRUE)
  cVol <- ClusteredNeuroVol(mask, cMap)

  # build e.g. 2 scans => each is NeuroVec
  # say each has 2 timepoints
  arrA <- array(rnorm(prod(dims_3d)*2), dim=c(dims_3d,2))
  arrB <- array(rnorm(prod(dims_3d)*2), dim=c(dims_3d,2))
  spA  <- NeuroSpace(c(dims_3d,2))
  spB  <- NeuroSpace(c(dims_3d,2))

  vecA <- NeuroVec(arrA, spA)
  vecB <- NeuroVec(arrB, spB)
  scans_list <- list(vecA, vecB)
  scan_names <- c("scanX","scanY")
  # metadata
  smd_list <- list(
    list(TR=1.5, subject="Bob"),
    list(TR=2.0, subject="Sue")
  )

  tmpf <- tempfile(fileext=".h5")
  on.exit(unlink(tmpf))

  wobj <- write_clustered_dataset(
    file=file.path(tmpf),
    vecs=scans_list,
    scan_names=scan_names,
    mask=mask,
    clusters=cVol,
    scan_metadata=smd_list,
    summary_only=FALSE
  )
  expect_s4_class(wobj, "H5ClusteredVecSeq")

  # read back => prefer_reduced=FALSE
  rseq <- read_clustered_dataset(tmpf, prefer_reduced=FALSE)
  expect_s4_class(rseq, "H5ClusteredVecSeq")

  # check times
  # for each scan => dimension => [ #vox in cluster, #time ] => or partial reading
  # We'll do linear_access or [i=..]

  # check single scan => rseq[[1]] => H5ClusteredVec
  hv1 <- rseq[[1]]
  expect_s4_class(hv1, "H5ClusteredVec")
  # partial sub => i=1:3, j=2:3, k=2, l=1 => compare arrA
  sub_mem <- arrA[1:3,2:3,2,1, drop=FALSE]
  sub_h5  <- hv1[1:3, 2:3, 2, 1]
  # sub_h5 => shape c(3,2,1,1) if drop=FALSE, or c(3,2) if drop=TRUE
  # depending on how you wrote your subset method.
  dim_check <- c(3,2)
  expect_equal(dim(sub_h5), dim_check)
  sub_mem2d <- drop(sub_mem)  # => c(3,2)
  expect_equal(sub_h5, sub_mem2d, tolerance=1e-7)

  # linear access => pick random
  totA <- prod(dims_3d)*2
  idx_samp <- c(1,2,9,10, totA)
  hvvals   <- hv1[idx_samp]
  memvalsA <- arrA[idx_samp]
  expect_equal(hvvals, memvalsA, tolerance=1e-7)
})
