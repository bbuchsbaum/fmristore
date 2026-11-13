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

# --- `[` on H5ParcellatedScan with j missing -------------------------------

test_that("[ on H5ParcellatedScan keeps k and l when j is missing", {
  fx <- make_parcel_fixture()
  tf <- tempfile(fileext = ".h5")
  on.exit(unlink(tf), add = TRUE)
  write_parcel_fixture(fx, tf)
  exp <- open_parcel_exp(tf)
  on.exit(close(exp), add = TRUE, after = FALSE)
  run <- exp@runs[[1]]

  expect_equal(run[1:3, , , 2], fx$ref[1:3, , , 2])
  expect_equal(run[, , , 1], fx$ref[, , , 1])
  expect_equal(run[2, , 2, ], fx$ref[2, , 2, ])
  expect_equal(run[, , 1:2, 3:4], fx$ref[, , 1:2, 3:4])
  expect_equal(run[4, , ], fx$ref[4, , , ])

  # A single subscript still means mask indices.
  expect_equal(run[c(1, 2)], fx$full[c(1, 2), ])
  expect_equal(run[c(9, 1), drop = FALSE], fx$full[c(9, 1), , drop = FALSE])
  # Fully specified subscripts are unchanged.
  expect_equal(run[1:4, 1:3, 1:2, 2], fx$ref[, , , 2])
})

# --- overwrite = TRUE must not destroy the file on invalid input ------------

test_that("overwrite = TRUE keeps the existing file when validation fails", {
  fx <- make_parcel_fixture()
  tf <- tempfile(fileext = ".h5")
  on.exit(unlink(tf), add = TRUE)
  write_parcel_fixture(fx, tf)
  before <- unname(tools::md5sum(tf))
  overwrite_with <- function(runs, ...) {
    write_parcellated_experiment_h5(tf, fx$mask, fx$cvol, runs,
      overwrite = TRUE, verbose = FALSE, ...
    )
  }
  unchanged <- function() {
    file.exists(tf) && identical(unname(tools::md5sum(tf)), before)
  }

  bad_names <- list(list(
    scan_name = "run1", type = "full",
    data = list(bad_name = matrix(1, 1, 1))
  ))
  expect_error(overwrite_with(bad_names), "invalid cluster names")
  expect_true(unchanged())

  expect_error(
    overwrite_with(list(parcel_full_run(fx)),
      cluster_metadata = data.frame(cluster_id = 1:2)
    ),
    "Cluster IDs in metadata"
  )
  expect_true(unchanged())

  dl <- parcel_cluster_list(fx)
  dl$cluster_3 <- dl$cluster_3[1:2, , drop = FALSE]
  expect_error(overwrite_with(list(parcel_full_run(fx, data = dl))), "cluster_3")
  expect_true(unchanged())

  bad_summary <- list(list(
    scan_name = "s", type = "summary", data = matrix(0, 5, 2)
  ))
  expect_error(overwrite_with(bad_summary), "columns")
  expect_true(unchanged())

  dup <- list(parcel_full_run(fx, "a"), parcel_full_run(fx, "a"))
  expect_error(overwrite_with(dup), "unique")
  expect_true(unchanged())
})

# --- Full-data writer validation ---------------------------------------------

test_that("full-data writer validates matrices against clusters and metadata", {
  fx <- make_parcel_fixture()
  tf <- tempfile(fileext = ".h5")
  on.exit(unlink(tf), add = TRUE)
  write_run <- function(rd) {
    write_parcellated_experiment_h5(tf, fx$mask, fx$cvol, list(rd),
      overwrite = TRUE, verbose = FALSE
    )
  }

  expect_error(
    write_run(parcel_full_run(fx, meta = list(n_time = 3L))),
    "n_time \\(3\\).*5 columns"
  )

  dl <- parcel_cluster_list(fx)
  dl$cluster_10 <- NULL
  expect_error(write_run(parcel_full_run(fx, data = dl)), "missing.*cluster_10")

  dl <- parcel_cluster_list(fx)
  dl$cluster_3 <- dl$cluster_3[1:2, , drop = FALSE]
  expect_error(
    write_run(parcel_full_run(fx, data = dl)),
    "cluster_3.*2 rows.*4 voxels"
  )

  dl <- parcel_cluster_list(fx)
  dl$cluster_7 <- cbind(dl$cluster_7, 0)
  expect_error(
    write_run(parcel_full_run(fx, data = dl)),
    "same number of columns"
  )

  dl <- parcel_cluster_list(fx)
  dl$cluster_99 <- matrix(0, 1, fx$n_time)
  expect_error(write_run(parcel_full_run(fx, data = dl)), "cluster_99")

  dl <- parcel_cluster_list(fx)
  dl$cluster_7 <- array(as.character(dl$cluster_7), dim(dl$cluster_7))
  expect_error(write_run(parcel_full_run(fx, data = dl)), "numeric matrix")

  # Valid data, in any order and with matching metadata, still round-trips.
  dl <- rev(parcel_cluster_list(fx))
  expect_no_error(
    write_run(parcel_full_run(fx, meta = list(n_time = 5L), data = dl))
  )
  exp <- open_parcel_exp(tf)
  on.exit(close(exp), add = TRUE, after = FALSE)
  expect_equal(exp@runs[[1]]@n_time, 5L)
  expect_equal(unname(series(exp@runs[[1]], seq_len(9))), t(fx$full))
})

# --- Spatial metadata and header round trip ---------------------------------

test_that("parcellated writers preserve origin, spacing and orientation", {
  ang <- pi / 6
  rot <- matrix(c(cos(ang), sin(ang), 0, -sin(ang), cos(ang), 0, 0, 0, 1), 3)
  oblique <- diag(4)
  oblique[1:3, 1:3] <- rot %*% diag(c(-2, 2.5, 3)) # rotated, x flipped
  oblique[1:3, 4] <- c(12.5, -7, 3)

  fixtures <- list(
    shifted = make_parcel_fixture(origin = c(-10, -20, -30)),
    oblique = make_parcel_fixture(trans = oblique)
  )
  for (nm in names(fixtures)) {
    fx <- fixtures[[nm]]
    tf_exp <- tempfile(fileext = ".h5")
    tf_scan <- tempfile(fileext = ".h5")
    write_parcel_fixture(fx, tf_exp)
    run <- write_parcellated_scan_h5(tf_scan, fx$mask, fx$cvol,
      parcel_cluster_list(fx),
      data_type = "full"
    )
    close(run)

    for (tf in c(tf_exp, tf_scan)) {
      exp <- open_parcel_exp(tf)
      sp_read <- space(mask(exp))
      close(exp)
      expect_equal(trans(sp_read), trans(fx$sp), tolerance = 1e-6, info = nm)
      expect_equal(origin(sp_read), origin(fx$sp), tolerance = 1e-6, info = nm)
      expect_equal(spacing(sp_read), spacing(fx$sp), tolerance = 1e-6, info = nm)

      h <- hdf5r::H5File$new(tf, "r")
      hdr_dim <- h[["header/dim"]][]
      h$close_all()
      expect_equal(hdr_dim[5], fx$n_time, info = nm)
    }
    unlink(c(tf_exp, tf_scan))
  }
})

# --- `$` accessor --------------------------------------------------------------

test_that("$ on H5ParcellatedMultiScan returns shared run properties", {
  fx <- make_parcel_fixture()
  tf <- tempfile(fileext = ".h5")
  on.exit(unlink(tf), add = TRUE)
  write_parcel_fixture(fx, tf)
  exp <- open_parcel_exp(tf)
  on.exit(close(exp), add = TRUE, after = FALSE)

  expect_identical(exp$mask, mask(exp))
  expect_identical(exp$clusters, clusters(exp))
  expect_true(inherits(exp$obj, "H5File"))
  expect_identical(exp$h5file, exp$obj)
  expect_true(exp$h5file$is_valid)
  expect_equal(exp$n_voxels, 9L)
  expect_identical(exp$runs, exp@runs)
  expect_null(exp$not_a_property)
})
