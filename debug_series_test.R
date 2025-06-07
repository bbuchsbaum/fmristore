library(fmristore)
library(neuroim2)
library(Matrix)

# Create simple test case
dims <- c(4, 4, 2, 5)
sp <- NeuroSpace(dims)

# Create mask - let's be very explicit
mask_arr <- array(FALSE, dim = dims[1:3])
mask_arr[2:3, 2:3, 1:2] <- TRUE  # Small central region
mask_vol <- LogicalNeuroVol(mask_arr, drop_dim(sp))
nVox_mask <- sum(mask_vol)

print(paste("Number of masked voxels:", nVox_mask))
print("Mask indices:")
print(which(mask_arr, arr.ind = TRUE))

# Create simple data
set.seed(123)
original_data <- array(rnorm(prod(dims), mean=10, sd=2), dims)

# Create LatentNeuroVec
# For testing, let's use a simple full-rank decomposition
masked_indices <- which(mask_arr, arr.ind = TRUE)
masked_data <- matrix(0, nrow = nVox_mask, ncol = dims[4])

for (i in seq_len(nVox_mask)) {
  x <- masked_indices[i, 1]
  y <- masked_indices[i, 2] 
  z <- masked_indices[i, 3]
  masked_data[i, ] <- original_data[x, y, z, ]
}

print("Masked data dimensions:")
print(dim(masked_data))
print("First few values of masked data:")
print(masked_data[1:min(3, nrow(masked_data)), ])

# Use identity decomposition for exact match
k <- min(dims[4], nVox_mask)
basis_mat <- Matrix(diag(k), sparse = FALSE)  # Identity matrix
loadings_mat <- Matrix(masked_data[1:k, ], sparse = FALSE)  # Take first k voxels
offset_vec <- rep(0, k)

# Adjust if we have fewer voxels than time points
if (nVox_mask < dims[4]) {
  basis_mat <- Matrix(diag(nVox_mask), nrow = dims[4], ncol = nVox_mask, sparse = FALSE)
  loadings_mat <- Matrix(masked_data, sparse = FALSE)
  offset_vec <- rep(0, nVox_mask)
}

print("Basis dimensions:")
print(dim(basis_mat))
print("Loadings dimensions:")
print(dim(loadings_mat))

# Create LatentNeuroVec
lvec <- LatentNeuroVec(basis = basis_mat,
                      loadings = loadings_mat,
                      space = sp,
                      mask = mask_vol,
                      offset = offset_vec,
                      label = "debug_test")

# Create DenseNeuroVec
dvec <- DenseNeuroVec(original_data, sp)

print("\n=== Testing series extraction ===")

# Test a known masked voxel
test_coords <- c(2L, 2L, 1L)
print(paste("Testing coordinates:", paste(test_coords, collapse=", ")))
print(paste("Is in mask:", mask_vol[test_coords[1], test_coords[2], test_coords[3]]))

latent_series <- series(lvec, test_coords[1], test_coords[2], test_coords[3])
dense_series <- series(dvec, test_coords[1], test_coords[2], test_coords[3])

print("LatentNeuroVec series:")
print(latent_series)
print("DenseNeuroVec series:")  
print(dense_series)
print("Original data at this voxel:")
print(original_data[test_coords[1], test_coords[2], test_coords[3], ])

# Check if the issue is with reconstruction
print("\n=== Testing manual reconstruction ===")
# Manual reconstruction: loadings[voxel_idx, ] %*% t(basis) + offset[voxel_idx]
# Find which masked voxel index corresponds to coordinates (2,2,1)
voxel_idx <- which(apply(masked_indices, 1, function(row) all(row == test_coords)))
print(paste("Voxel index in masked data:", voxel_idx))

if (length(voxel_idx) > 0) {
  manual_reconstruction <- as.vector(loadings_mat[voxel_idx, ] %*% t(basis_mat)) + offset_vec[voxel_idx]
  print("Manual reconstruction:")
  print(manual_reconstruction)
}

# Check the as.array method
print("\n=== Testing as.array reconstruction ===")
reconstructed_array <- as.array(lvec)
print("Reconstructed array at test coordinates:")
print(reconstructed_array[test_coords[1], test_coords[2], test_coords[3], ]) 