# Test suite for scale and performance critical operations
# This addresses gaps in testing large-scale operations and performance characteristics

library(fmristore)

test_that("H5ClusterExperiment handles large numbers of scans efficiently", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("neuroim2")
  skip_on_cran()  # These tests may be too intensive for CRAN
  
  # Test 1: Create experiment with many scans
  temp_file <- tempfile(fileext = ".h5")
  on.exit(unlink(temp_file), add = TRUE)
  
  # Parameters for scale test
  n_scans <- 50  # Moderate number to test scaling
  mask_dims <- c(10L, 10L, 5L)
  n_time_per_scan <- 20L
  n_clusters <- 5L
  
  # Create base objects
  mask <- fmristore:::create_minimal_LogicalNeuroVol(mask_dims)
  clusters <- fmristore:::create_minimal_ClusteredNeuroVol(mask, num_clusters = n_clusters)
  
  # Generate scan data
  scan_list <- list()
  scan_names <- character(n_scans)
  scan_metadata <- list()
  
  # Time the creation of many scans
  create_time <- system.time({
    for (i in seq_len(n_scans)) {
      scan_name <- sprintf("scan_%03d", i)
      scan_names[i] <- scan_name
      
      # Create scan data - alternating between full and summary
      if (i %% 2 == 1) {
        # Full scan
        scan_list[[scan_name]] <- list(
          type = "full",
          data = lapply(seq_len(n_clusters), function(cid) {
            n_vox <- sum(clusters@clusters == cid)
            matrix(rnorm(n_vox * n_time_per_scan), nrow = n_vox, ncol = n_time_per_scan)
          })
        )
      } else {
        # Summary scan
        scan_list[[scan_name]] <- list(
          type = "summary",
          data = matrix(rnorm(n_time_per_scan * n_clusters), 
                       nrow = n_time_per_scan, ncol = n_clusters)
        )
      }
      
      scan_metadata[[scan_name]] <- list(
        subject_id = sprintf("sub_%02d", (i-1) %/% 5 + 1),
        session = sprintf("ses_%02d", (i-1) %% 5 + 1),
        task = sample(c("rest", "motor", "visual"), 1),
        run = i,
        TR = 2.0
      )
    }
  })
  
  expect_true(create_time["elapsed"] < 60, 
              info = sprintf("Creating %d scans took %.1f seconds", n_scans, create_time["elapsed"]))
  
  # Write to HDF5
  write_time <- system.time({
    write_clustered_experiment_h5(
      file = temp_file,
      runs_data = scan_list,
      scan_names = scan_names,
      mask = mask,
      clusters = clusters,
      scan_metadata = scan_metadata,
      compression = 4  # Moderate compression for speed
    )
  })
  
  expect_true(write_time["elapsed"] < 120,
              info = sprintf("Writing %d scans took %.1f seconds", n_scans, write_time["elapsed"]))
  
  # Test 2: Efficient access to specific scans
  exp_obj <- H5ClusterExperiment(temp_file)
  
  # Time accessing individual scans
  access_times <- numeric(10)
  for (i in 1:10) {
    scan_idx <- sample(n_scans, 1)
    scan_name <- scan_names[scan_idx]
    
    access_times[i] <- system.time({
      scan_data <- exp_obj[[scan_name]]
      # Perform basic operation
      if (inherits(scan_data, "H5ClusterRun")) {
        dim_check <- dim(scan_data)
      } else {
        mat_check <- as.matrix(scan_data)
      }
    })["elapsed"]
  }
  
  expect_true(mean(access_times) < 1.0,
              info = sprintf("Average scan access time: %.3f seconds", mean(access_times)))
  
  # Test 3: Memory efficiency with subset operations
  gc()  # Clear memory
  mem_before <- gc()[2, 2]  # Used memory in MB
  
  # Access multiple scans without storing all data
  subset_results <- lapply(1:20, function(i) {
    scan_name <- scan_names[i]
    scan_obj <- exp_obj[[scan_name]]
    
    if (inherits(scan_obj, "H5ClusterRun")) {
      # Extract just a few time series
      series(scan_obj, i = 1:5)
    } else {
      # Get summary for one cluster
      as.matrix(scan_obj)[, 1]
    }
  })
  
  gc()
  mem_after <- gc()[2, 2]
  mem_increase <- mem_after - mem_before
  
  expect_true(mem_increase < 100,
              info = sprintf("Memory increase for 20 scan operations: %.1f MB", mem_increase))
  
  close(exp_obj)
  
  # Test 4: File size efficiency
  file_size_mb <- file.info(temp_file)$size / 1024^2
  expected_size_mb <- n_scans * n_time_per_scan * mean(c(sum(mask@.Data), n_clusters)) * 8 / 1024^2
  
  expect_true(file_size_mb < expected_size_mb,
              info = sprintf("Compressed size (%.1f MB) should be less than uncompressed (%.1f MB)",
                           file_size_mb, expected_size_mb))
})

test_that("Chunking strategies affect performance appropriately", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("neuroim2")
  skip_on_cran()
  
  # Create test data
  dims <- c(50L, 50L, 30L, 100L)  # Larger dimensions for chunking tests
  vec_data <- array(rnorm(prod(dims[1:3]) * 10), dim = c(dims[1:3], 10L))  # Subset for speed
  space_obj <- neuroim2::NeuroSpace(dims[1:3])
  vec <- neuroim2::NeuroVec(vec_data, space_obj)
  
  chunk_configs <- list(
    # Different chunking strategies
    time_optimized = c(10, 10, 10, dims[4]),     # Full time dimension
    space_optimized = c(dims[1], dims[2], dims[3], 1),  # Single time point
    balanced = c(25, 25, 15, 20),                # Balanced chunks
    small = c(5, 5, 5, 10),                      # Small chunks
    large = c(50, 50, 30, 50)                    # Large chunks
  )
  
  results <- list()
  
  for (chunk_name in names(chunk_configs)) {
    temp_file <- tempfile(fileext = ".h5")
    chunk_dim <- chunk_configs[[chunk_name]]
    
    # Time the write operation
    write_time <- system.time({
      h5_vec <- as_h5(vec, file = temp_file, 
                      chunk_dim = chunk_dim[1:4],
                      compression = 4)
      close(h5_vec)
    })["elapsed"]
    
    # Time different access patterns
    h5_vec <- H5NeuroVec(temp_file)
    
    # 1. Time series access (should favor time_optimized)
    ts_time <- system.time({
      for (i in 1:20) {
        idx <- sample(prod(dims[1:3]), 1)
        ts <- series(h5_vec, idx)
      }
    })["elapsed"]
    
    # 2. Spatial slice access (should favor space_optimized)
    slice_time <- system.time({
      for (t in 1:5) {
        slice <- h5_vec[, , , t]
      }
    })["elapsed"]
    
    # 3. Random access
    random_time <- system.time({
      for (i in 1:50) {
        x <- sample(dims[1], 1)
        y <- sample(dims[2], 1) 
        z <- sample(dims[3], 1)
        t <- sample(10, 1)  # Using subset
        val <- h5_vec[x, y, z, t]
      }
    })["elapsed"]
    
    close(h5_vec)
    
    results[[chunk_name]] <- list(
      write_time = write_time,
      ts_time = ts_time,
      slice_time = slice_time,
      random_time = random_time,
      file_size = file.info(temp_file)$size / 1024^2
    )
    
    unlink(temp_file)
  }
  
  # Verify performance characteristics
  # Time-optimized should be fastest for time series
  ts_times <- sapply(results, function(x) x$ts_time)
  expect_true(
    ts_times["time_optimized"] <= min(ts_times) * 1.5,
    info = "Time-optimized chunking should be among fastest for time series access"
  )
  
  # Space-optimized should be fastest for spatial slices
  slice_times <- sapply(results, function(x) x$slice_time)
  expect_true(
    slice_times["space_optimized"] <= min(slice_times) * 1.5,
    info = "Space-optimized chunking should be among fastest for spatial slice access"
  )
  
  # Small chunks should have larger file size due to overhead
  file_sizes <- sapply(results, function(x) x$file_size)
  expect_true(
    file_sizes["small"] > file_sizes["balanced"],
    info = "Small chunks should result in larger file size due to metadata overhead"
  )
})

test_that("Large LatentNeuroVec operations scale appropriately", {
  skip_if_not_installed("Matrix")
  skip_if_not_installed("neuroim2")
  skip_on_cran()
  
  # Test scaling with increasing data sizes
  test_configs <- list(
    small = list(dims = c(20L, 20L, 10L), n_time = 50L, n_comp = 10L),
    medium = list(dims = c(40L, 40L, 20L), n_time = 100L, n_comp = 20L),
    large = list(dims = c(60L, 60L, 30L), n_time = 200L, n_comp = 40L)
  )
  
  scaling_results <- list()
  
  for (size_name in names(test_configs)) {
    config <- test_configs[[size_name]]
    
    # Create mask covering ~30% of volume
    mask_array <- array(FALSE, dim = config$dims)
    n_voxels <- prod(config$dims)
    mask_indices <- sample(n_voxels, ceiling(n_voxels * 0.3))
    mask_array[mask_indices] <- TRUE
    mask <- neuroim2::LogicalNeuroVol(mask_array, neuroim2::NeuroSpace(config$dims))
    n_mask_voxels <- sum(mask@.Data)
    
    # Create basis and loadings
    create_time <- system.time({
      basis <- matrix(rnorm(config$n_time * config$n_comp), 
                     nrow = config$n_time, ncol = config$n_comp)
      
      # Create sparse loadings with ~10% density
      n_nonzero <- ceiling(n_mask_voxels * config$n_comp * 0.1)
      loadings <- Matrix::sparseMatrix(
        i = sample(n_mask_voxels, n_nonzero, replace = TRUE),
        j = sample(config$n_comp, n_nonzero, replace = TRUE),
        x = rnorm(n_nonzero),
        dims = c(n_mask_voxels, config$n_comp)
      )
      
      lnv <- LatentNeuroVec(basis = basis, loadings = loadings, mask = mask)
    })["elapsed"]
    
    # Test reconstruction time
    reconstruct_time <- system.time({
      # Reconstruct subset of voxels
      subset_indices <- sample(which(mask@.Data), min(100, n_mask_voxels))
      series_data <- series(lnv, subset_indices)
    })["elapsed"]
    
    # Test write time
    temp_file <- tempfile(fileext = ".h5")
    write_time <- system.time({
      write_vec(lnv, temp_file, compression = 4)
    })["elapsed"]
    
    file_size_mb <- file.info(temp_file)$size / 1024^2
    unlink(temp_file)
    
    # Test component extraction time
    component_time <- system.time({
      # Extract first 5 components
      comps <- components(lnv)
      comp_subset <- comps[1:min(5, length(comps))]
    })["elapsed"]
    
    scaling_results[[size_name]] <- list(
      n_voxels = n_mask_voxels,
      n_elements = n_mask_voxels * config$n_time,
      create_time = create_time,
      reconstruct_time = reconstruct_time,
      write_time = write_time,
      component_time = component_time,
      file_size_mb = file_size_mb
    )
  }
  
  # Check that operations scale reasonably (not worse than O(n^2))
  sizes <- sapply(scaling_results, function(x) x$n_elements)
  create_times <- sapply(scaling_results, function(x) x$create_time)
  reconstruct_times <- sapply(scaling_results, function(x) x$reconstruct_time)
  
  # Fit linear model to log-log plot to estimate scaling exponent
  if (length(sizes) >= 3) {
    create_model <- lm(log(create_times) ~ log(sizes))
    create_exponent <- coef(create_model)[2]
    
    reconstruct_model <- lm(log(reconstruct_times) ~ log(sizes))
    reconstruct_exponent <- coef(reconstruct_model)[2]
    
    # Exponent should be close to 1 for linear scaling, definitely less than 2
    expect_true(create_exponent < 1.5,
                info = sprintf("Creation scaling exponent: %.2f (should be < 1.5)", create_exponent))
    
    expect_true(reconstruct_exponent < 1.5,
                info = sprintf("Reconstruction scaling exponent: %.2f (should be < 1.5)", 
                             reconstruct_exponent))
  }
  
  # Memory efficiency check
  # File size should scale linearly with data size due to compression and sparsity
  file_sizes <- sapply(scaling_results, function(x) x$file_size_mb)
  size_model <- lm(file_sizes ~ sizes)
  size_r_squared <- summary(size_model)$r.squared
  
  expect_true(size_r_squared > 0.9,
              info = sprintf("File size should scale linearly with data size (R²=%.3f)", size_r_squared))
})