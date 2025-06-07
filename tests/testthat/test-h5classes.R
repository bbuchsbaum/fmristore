
library(neuroim2)
library(fmristore)

test_that("H5NeuroVol construction and subsetting works", {
  # Create test data
  dims <- c(10, 10, 10)
  space <- NeuroSpace(dims, spacing=c(2,2,2), origin=c(1,1,1))
  data <- array(rnorm(prod(dims)), dim=dims)
  vol <- NeuroVol(data, space)

  # Create temporary HDF5 file
  tmp <- tempfile(fileext=".h5")
  on.exit(unlink(tmp))

  # Convert to H5NeuroVol
  h5vol <- as_h5(vol, tmp)

  # Test basic properties
  expect_true(inherits(h5vol, "H5NeuroVol"))
  expect_equal(dim(h5vol), dims)
  expect_equal(spacing(space(h5vol)), spacing(space(vol)))
  expect_equal(origin(space(h5vol)), origin(space(vol)))

  # Test subsetting
  subset <- h5vol[1:5, 1:5, 1:5]
  expect_equal(dim(subset), c(5,5,5))
  expect_equal(as.array(subset), data[1:5, 1:5, 1:5], tolerance=1e-6)
})

test_that("coercion via as() uses as_h5 for DenseNeuroVol", {
  dims <- c(4, 4, 4)
  sp   <- NeuroSpace(dims)
  arr  <- array(rnorm(prod(dims)), dim = dims)
  dvol <- DenseNeuroVol(arr, sp)

  h5vol <- as(dvol, "H5NeuroVol")

  expect_true(inherits(h5vol, "H5NeuroVol"))
  expect_equal(dim(h5vol), dims)
  expect_equal(h5vol[2,2,2], arr[2,2,2], tolerance = 1e-6)

  filepath <- h5vol@h5obj$get_filename()
  close(h5vol)
  unlink(filepath)
})

test_that("H5NeuroVec construction and series extraction works", {
  # Create test data
  dims <- c(10, 10, 10)
  nvols <- 5
  space <- NeuroSpace(c(dims, nvols), spacing=c(2,2,2), origin=c(1,1,1))
  data <- array(rnorm(prod(dims) * nvols), dim=c(dims, nvols))
  vec <- NeuroVec(data, space)

  # Create temporary HDF5 file
  tmp <- tempfile(fileext=".h5")
  on.exit(unlink(tmp))

  # Convert to H5NeuroVec
  h5vec <- as_h5(vec, tmp)

  # Test basic properties
  expect_true(inherits(h5vec, "H5NeuroVec"))
  expect_equal(dim(h5vec), c(dims, nvols))
  expect_equal(spacing(space(h5vec)), spacing(space(vec)))
  expect_equal(origin(space(h5vec)), origin(space(vec)))
  expect_equal(trans(space(h5vec)), trans(space(vec)))

  # Test series extraction
  series1 <- series(h5vec, 5, 5, 5)
  expect_equal(length(series1), nvols)
  expect_equal(series1, data[5,5,5,], tolerance=1e-6)

  # Test subsetting
  subset <- h5vec[1:5, 1:5, 1:5, 1:3]
  expect_equal(dim(subset), c(5,5,5,3))
  expect_equal(as.array(subset), data[1:5, 1:5, 1:5, 1:3], tolerance=1e-6)
})

test_that("LatentNeuroVec single timepoint reconstruction", {
  n_basis       <- 5
  n_voxels      <- 1000
  n_timepoints  <- 100

  # basis => (n_timepoints x n_basis) => (100 x 5)
  basis <- Matrix(rnorm(n_timepoints * n_basis),
                  nrow = n_timepoints,
                  ncol = n_basis)

  # loadings => (p x k) => (1000 x 5)
  loadings <- Matrix(rnorm(n_voxels * n_basis),
                     nrow = n_voxels,
                     ncol = n_basis,
                     sparse = TRUE)

  offset <- rnorm(n_voxels)  # length=1000

  # Create a NeuroSpace with dims (10,10,10,100)
  space <- NeuroSpace(c(10,10,10,n_timepoints), spacing=c(2,2,2), origin=c(1,1,1))
  mask_vol <- LogicalNeuroVol(array(TRUE, dim=c(10,10,10)), drop_dim(space))

  latent_vec <- LatentNeuroVec(
    basis    = basis,      # (100 x 5)
    loadings = loadings,   # (1000 x 5)
    space    = space,
    mask     = mask_vol,
    offset   = offset
  )

  # Test basic properties
  expect_true(inherits(latent_vec, "LatentNeuroVec"))
  expect_equal(dim(latent_vec), c(10,10,10,n_timepoints))

  # Fix: Multiply as (1×n_basis)×(n_basis×n_voxels) => (1×n_voxels)
  # "time_point_1" is the full volume at time=1, i.e. basis[1, ] times loadings plus offset
  # Instead of basis[1,,drop=FALSE] %*% loadings:
  time_point_1 <- as.vector( basis[1,,drop=FALSE] %*% t(loadings) + offset )
  #time_point_1 <- as.vector(basis[1,,drop=FALSE] %*% loadings + offset)
  # Reconstruct time=1 from the latent vector
  reconstructed_1 <- as.vector(latent_vec[,,,1])

  # They should match
  expect_equal(time_point_1, reconstructed_1, tolerance=1e-8)

  # Test series extraction
  # e.g. single voxel's time series
  series1 <- series(latent_vec, 5, 5, 5)
  expect_equal(length(series1), n_timepoints)
})


test_that("H5NeuroVol indexing matches in-memory NeuroVol", {
  # 1) Create a random 3D volume in memory
  set.seed(42)
  arr <- array(rnorm(1000), dim = c(10, 10, 10))
  sp  <- NeuroSpace(dim=c(10, 10, 10))
  mem_vol <- NeuroVol(arr, sp)  # DenseNeuroVol

  # 2) Convert it to H5NeuroVol on the fly
  tmpfile <- tempfile(fileext = ".h5")
  h5vol   <- as_h5(mem_vol, file=tmpfile, data_type="FLOAT")

  # 3) Check basic dimension
  expect_equal(dim(h5vol), c(10,10,10))

  # 4) Test single voxel access
  expect_equal(h5vol[1,1,1], arr[1,1,1], tolerance=1e-7)
  expect_equal(h5vol[10,10,10], arr[10,10,10], tolerance=1e-7)

  # 5) Test slices
  slice_2 <- h5vol[1:5, 3:4, 2]
  expect_equal(slice_2, arr[1:5, 3:4, 2], tolerance=1e-6)

  # 6) Test linear indexing
  lin_idx <- c(1, 10, 100, 500, 999)
  h5_vals <- h5vol[lin_idx]
  mem_vals <- arr[lin_idx]
  expect_equal(h5_vals, mem_vals, tolerance=1e-6)

  # 7) Test matrix-based subsetting (arrayInd style)
  coords <- matrix(c(
    1,1,1,
    2,3,4,
    10,10,10
  ), ncol=3, byrow=TRUE)
  expect_equal(h5vol[coords], arr[coords], tolerance=1e-6)

  # 8) Random index checks
  i <- sample(1:10, 5)
  j <- sample(1:10, 4)
  k <- sample(1:10, 3)
  expect_equal(h5vol[i,j,k], arr[i,j,k], tolerance=1e-6)

  # 9) Check partial missing indices (like h5vol[i, , ])
  i2 <- 2:5
  sub_arr <- h5vol[i2,,]
  expect_equal(sub_arr, arr[i2,,], tolerance=1e-6)

  # 10) Cleanup
  hdf5r::h5close(h5vol@h5obj)
  file.remove(tmpfile)
})


test_that("H5NeuroVec indexing matches in-memory DenseNeuroVec", {

  ## 1) Construct a 4D in-memory DenseNeuroVec
  set.seed(42)
  # e.g. 10x10x5 volume with 8 timepoints => total dims c(10,10,5,8)
  arr <- array(rnorm(10*10*5*8), dim=c(10,10,5,8))
  sp  <- NeuroSpace(dim=c(10,10,5,8))
  mem_vec <- DenseNeuroVec(data=arr, space=sp)

  ## 2) Convert it to an H5NeuroVec file
  tmpfile <- tempfile(fileext=".h5")
  h5_vec <- as_h5(mem_vec, file=tmpfile, data_type="FLOAT",
                        chunk_dim=c(4,4,4,8),
                        compression=6)

  ## Basic checks
  expect_s4_class(h5_vec, "H5NeuroVec")
  expect_equal(dim(h5_vec), c(10,10,5,8))

  ## 3) Test straightforward subsetting
  # e.g. a small sub-block
  sub_arr <- mem_vec[1:5, 1:4, 2:3, 1:2]
  h5sub   <- h5_vec[1:5, 1:4, 2:3, 1:2]
  expect_equal(h5sub, sub_arr, tolerance=1e-7)

  # Another sub-block
  sub_arr2 <- mem_vec[2:3, 5:10, 1:5, 8]
  h5sub2   <- h5_vec[2:3, 5:10, 1:5, 8]
  expect_equal(h5sub2, sub_arr2, tolerance=1e-7)

  ## 4) Test partial indexing (missing one dimension => full range):
  # example: x[i, j, k] => implies full range in 4th dimension
  sub3d_mem <- mem_vec[2:3, 1:2, 4:5,]  # i.e. ignoring time => all volumes
  sub3d_h5  <- h5_vec[2:3, 1:2, 4:5,]
  expect_equal(sub3d_mem, sub3d_h5, tolerance=1e-7)

  ## 5) Test linear indexing
  # pick some scattered linear indices, compare
  lin_idx <- c(1, 50, 100, 234, 999, 1000, 2000)
  lin_mem <- mem_vec[ lin_idx ]
  lin_h5  <- h5_vec[ lin_idx ]
  expect_equal(lin_h5, lin_mem, tolerance=1e-7)

  ## 6) Test more random scattered subsetting across dimensions
  # We'll define random i, j, k, l
  i_scat <- sample(1:10, 4)
  j_scat <- sample(1:10, 5)
  k_scat <- sample(1:5,  3)
  l_scat <- sample(1:8,  2)

  scat_mem <- mem_vec[i_scat, j_scat, k_scat, l_scat]
  scat_h5  <- h5_vec[i_scat, j_scat, k_scat, l_scat]
  expect_equal(scat_h5, scat_mem, tolerance=1e-7)

  ## 7) Test random linear indexing for 4D
  lin_scat <- sample(prod(dim(mem_vec)), 20) # 20 random voxels
  lin_mem2 <- mem_vec[ lin_scat ]
  lin_h5_2 <- h5_vec[ lin_scat ]
  expect_equal(lin_h5_2, lin_mem2, tolerance=1e-7)

  ## Cleanup
  h5_vec@obj$close()
  file.remove(tmpfile)
})
