# Regression tests for parcellated (clustered) HDF5 writers, readers and
# subsetting. Each block targets a specific bug that previously went
# unnoticed.

library(testthat)
library(neuroim2)

# Fixture: 4x3x2 volume, 9 mask voxels, non-contiguous cluster IDs 3/7/10.
# Voxel data value = mask_index * 100 + t, so every value identifies its voxel.
make_parcel_fixture <- function(origin = c(-10, -20, -30), n_time = 5L,
                                trans = NULL) {
  dims <- c(4L, 3L, 2L)
  sp <- if (is.null(trans)) {
    NeuroSpace(dims, spacing = c(2, 2, 2), origin = origin)
  } else {
    NeuroSpace(dims, trans = trans)
  }
  arr <- array(FALSE, dims)
  lin <- c(1, 2, 4, 6, 7, 9, 13, 18, 24)
  arr[lin] <- TRUE
  mask <- LogicalNeuroVol(arr, sp)
  cl <- c(7L, 3L, 7L, 3L, 3L, 7L, 7L, 3L, 10L)
  cvol <- ClusteredNeuroVol(mask, cl)
  full <- outer(seq_along(cl) * 100, seq_len(n_time), "+")
  ref <- array(0, c(dims, n_time))
  for (t in seq_len(n_time)) {
    v <- array(0, dims)
    v[lin] <- full[, t]
    ref[, , , t] <- v
  }
  list(
    mask = mask, cvol = cvol, cl = cl, lin = lin, full = full,
    n_time = n_time, dims = dims, sp = sp, ref = ref
  )
}

parcel_cluster_list <- function(fx, ids = sort(unique(fx$cl))) {
  dl <- lapply(ids, function(id) fx$full[which(fx$cl == id), , drop = FALSE])
  names(dl) <- paste0("cluster_", ids)
  dl
}

parcel_full_run <- function(fx, scan_name = "run1", meta = list(),
                            data = parcel_cluster_list(fx)) {
  list(scan_name = scan_name, type = "full", data = data, metadata = meta)
}

write_parcel_fixture <- function(fx, path, runs = list(parcel_full_run(fx))) {
  write_parcellated_experiment_h5(path, fx$mask, fx$cvol, runs,
    overwrite = TRUE, verbose = FALSE
  )
  invisible(path)
}

open_parcel_exp <- function(path, ...) {
  suppressMessages(H5ParcellatedMultiScan(path, ...))
}

# --- write_parcellated_scan_h5: full data matched by name --------------------

test_that("write_parcellated_scan_h5 assigns full data to clusters by name", {
  fx <- make_parcel_fixture()
  tf <- tempfile(fileext = ".h5")
  on.exit(unlink(tf), add = TRUE)

  dl <- parcel_cluster_list(fx)[c("cluster_10", "cluster_7", "cluster_3")]
  run <- write_parcellated_scan_h5(tf, fx$mask, fx$cvol, dl,
    scan_name = "s", data_type = "full"
  )
  on.exit(close(run), add = TRUE, after = FALSE)

  expect_equal(unname(series(run, seq_len(9))), t(fx$full))
  h <- hdf5r::H5File$new(tf, "r")
  on.exit(h$close_all(), add = TRUE, after = FALSE)
  for (id in c(3L, 7L, 10L)) {
    dset_dims <- h[[sprintf("scans/s/clusters/cluster_%d", id)]]$dims
    expect_equal(dset_dims, c(sum(fx$cl == id), fx$n_time))
  }
})

test_that("write_parcellated_scan_h5 rejects missing or mis-sized clusters", {
  fx <- make_parcel_fixture()
  tf <- tempfile(fileext = ".h5")
  on.exit(unlink(tf), add = TRUE)
  dl <- parcel_cluster_list(fx)
  write_full <- function(d) {
    write_parcellated_scan_h5(tf, fx$mask, fx$cvol, d, data_type = "full")
  }

  renamed <- dl
  names(renamed)[3] <- "cluster_11"
  expect_error(write_full(renamed), "cluster_10")

  short <- dl
  short$cluster_3 <- short$cluster_3[1:2, , drop = FALSE]
  expect_error(write_full(short), "cluster_3.*2 rows.*4 voxels")

  wide <- dl
  wide$cluster_7 <- cbind(wide$cluster_7, 0)
  expect_error(write_full(wide), "same number of columns")

  expect_error(write_full(unname(dl)), "cluster_<id>")
})
