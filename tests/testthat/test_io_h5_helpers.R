library(testthat)
library(hdf5r)
library(neuroim2)

# Test read_h5_mask_to_LogicalNeuroVol

test_that("read_h5_mask_to_LogicalNeuroVol reads and validates masks", {
  skip_if_not_installed("hdf5r")
  sp <- NeuroSpace(c(2,2,2))
  mask_arr <- array(c(TRUE,FALSE,FALSE,TRUE, TRUE,FALSE,TRUE,FALSE), dim=c(2,2,2))

  tmp <- tempfile(fileext=".h5")
  on.exit(unlink(tmp), add=TRUE)
  h5 <- H5File$new(tmp, mode="w")
  h5$create_dataset("mask", robj=as.integer(mask_arr))
  h5$create_dataset("mask_bad", robj=array(0L, dim=c(2,2,3)))
  h5$create_dataset("mask_2d", robj=array(0L, dim=c(2,2)))

  m <- read_h5_mask_to_LogicalNeuroVol(h5, "mask", sp)
  expect_s4_class(m, "LogicalNeuroVol")
  expect_equal(as.logical(as.array(m)), as.logical(mask_arr))

  expect_error(read_h5_mask_to_LogicalNeuroVol(h5, "mask_bad", sp),
               "Mask dimensions in HDF5")
  expect_error(read_h5_mask_to_LogicalNeuroVol(h5, "mask_2d", sp),
               "expected 3 dimensions")
  h5$close_all()
})

# Test read_h5_clusters_to_ClusteredNeuroVol

test_that("read_h5_clusters_to_ClusteredNeuroVol enforces length match", {
  skip_if_not_installed("hdf5r")
  sp <- NeuroSpace(c(2,2,1))
  mask_arr <- array(c(TRUE, TRUE, FALSE, TRUE), dim=c(2,2,1))
  mask <- LogicalNeuroVol(mask_arr, sp)

  tmp <- tempfile(fileext=".h5")
  on.exit(unlink(tmp), add=TRUE)
  h5 <- H5File$new(tmp, mode="w")
  h5$create_dataset("cluster_map", robj=1:sum(mask))
  h5$create_dataset("cluster_bad", robj=1:2)

  clus <- read_h5_clusters_to_ClusteredNeuroVol(h5, mask, "cluster_map")
  expect_s4_class(clus, "ClusteredNeuroVol")
  expect_equal(clus@clusters, 1:sum(mask))

  expect_error(read_h5_clusters_to_ClusteredNeuroVol(h5, mask, "cluster_bad"),
               "Length of")
  h5$close_all()
})
