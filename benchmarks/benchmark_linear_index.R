# Benchmark linear index calculations for LatentNeuroVec subsetting
#
# This script compares the previous expand.grid approach with the
# new loop-based method found in R/latent_vec.R.  It is intended to be
# run manually and requires the microbenchmark package.

library(microbenchmark)

# Old implementation using expand.grid
old_index_calc <- function(i, j, k, dims) {
  coords <- expand.grid(i, j, k)
  coords[,1] + (coords[,2]-1) * dims[1] + (coords[,3]-1) * dims[1] * dims[2]
}

# New implementation mirroring latent_vec.R
new_index_calc <- function(i, j, k, dims) {
  n <- length(i) * length(j) * length(k)
  lin <- integer(n)
  ij_grid <- outer(i, j, function(ii, jj) ii + (jj - 1L) * dims[1])
  step <- length(i) * length(j)
  for (kk in seq_along(k)) {
    start <- (kk - 1L) * step + 1L
    end <- kk * step
    lin[start:end] <- ij_grid + (k[kk] - 1L) * dims[1] * dims[2]
  }
  lin
}

# Example large ranges
I <- seq_len(200)
J <- seq_len(200)
K <- seq_len(150)
DIMS <- c(256, 256, 180)

bm <- microbenchmark(
  old = old_index_calc(I, J, K, DIMS),
  new = new_index_calc(I, J, K, DIMS),
  times = 10
)
print(bm)

