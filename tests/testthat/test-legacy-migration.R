library(neuroim2)

.migration_fixture <- function() {
  dims <- c(4L, 3L, 2L)
  trans <- diag(c(2, 2, 3, 1))
  trans[1:3, 4] <- c(-10, -20, -30)
  sp <- NeuroSpace(dims, trans = trans)
  arr <- array(FALSE, dims)
  arr[c(1, 2, 4, 6, 7, 9, 13, 18, 24)] <- TRUE
  mask <- LogicalNeuroVol(arr, sp)
  cl <- c(7L, 3L, 7L, 3L, 3L, 7L, 7L, 3L, 10L)
  list(mask = mask, cvol = ClusteredNeuroVol(mask, cl), cl = cl, ids = sort(unique(cl)), trans = trans)
}

.full_run <- function(fx, name, n_time, offset = 0) {
  full <- outer(seq_along(fx$cl) * 100 + offset, seq_len(n_time), "+")
  data <- lapply(fx$ids, function(id) full[fx$cl == id, , drop = FALSE])
  names(data) <- paste0("cluster_", fx$ids)
  list(scan_name = name, type = "full", data = data, full = full)
}

.summary_run <- function(fx, name, n_time, offset = 0) {
  data <- matrix(offset + seq_len(n_time * length(fx$ids)), n_time, length(fx$ids))
  list(scan_name = name, type = "summary", data = data)
}

test_that("summary scans become a parcel frame with runs in written order", {
  fx <- .migration_fixture()
  legacy <- tempfile(fileext = ".h5")
  frame_path <- tempfile(fileext = ".h5")
  on.exit(unlink(c(legacy, frame_path)), add = TRUE)
  runs <- list(.summary_run(fx, "run-2", 4L), .summary_run(fx, "run-10", 3L, offset = 1000))
  write_parcellated_experiment_h5(legacy, fx$mask, fx$cvol, runs, overwrite = TRUE, verbose = FALSE)

  migrate_parcellated_h5(legacy, frame_path)
  frame <- open_frame_h5(frame_path)
  sp <- fmridataset::space(frame)
  expect_s3_class(sp, "parcel_space")
  expect_identical(as.integer(Matrix::colSums(fmridataset::parcel_membership(sp))), as.integer(table(fx$cl)))
  expect_identical(unique(frame$observations$data$run_id), c("run-2", "run-10"))

  values <- fmridataset::collect_assay(frame, "summary")
  expect_equal(values[1:4, ], runs[[1]]$data, ignore_attr = TRUE, tolerance = 1e-6)
  expect_equal(values[5:7, ], runs[[2]]$data, ignore_attr = TRUE, tolerance = 1e-6)
  expect_equal(sp$parent$affine, fx$trans, ignore_attr = TRUE)
  expect_identical(sp$parent$support, which(fx$mask@.Data))
})

test_that("full scans become a voxel frame with a cluster_id feature column", {
  fx <- .migration_fixture()
  legacy <- tempfile(fileext = ".h5")
  frame_path <- tempfile(fileext = ".h5")
  on.exit(unlink(c(legacy, frame_path)), add = TRUE)
  run_a <- .full_run(fx, "a", 5L)
  run_b <- .full_run(fx, "b", 2L, offset = 5000)
  write_parcellated_experiment_h5(legacy, fx$mask, fx$cvol, list(run_a[1:4], run_b[1:4]),
    overwrite = TRUE, verbose = FALSE
  )

  migrate_parcellated_h5(legacy, frame_path)
  frame <- open_frame_h5(frame_path)
  expect_identical(frame$features$data$cluster_id, fx$cl)
  values <- fmridataset::collect_assay(frame, "signal")
  expect_equal(values, rbind(t(run_a$full), t(run_b$full)), ignore_attr = TRUE, tolerance = 1e-6)
})

test_that("mixed files are converted one kind at a time", {
  fx <- .migration_fixture()
  legacy <- tempfile(fileext = ".h5")
  frame_path <- tempfile(fileext = ".h5")
  on.exit(unlink(c(legacy, frame_path)), add = TRUE)
  write_parcellated_experiment_h5(legacy, fx$mask, fx$cvol,
    list(.full_run(fx, "full", 3L)[1:4], .summary_run(fx, "summ", 3L)),
    overwrite = TRUE, verbose = FALSE
  )
  expect_error(migrate_parcellated_h5(legacy, frame_path), "mixes full and summary")
  migrate_parcellated_h5(legacy, frame_path, kind = "summary")
  expect_identical(unique(open_frame_h5(frame_path)$observations$data$run_id), "summ")
  expect_error(migrate_parcellated_h5(legacy, frame_path, kind = "full"), "already exists")
  migrate_parcellated_h5(legacy, frame_path, kind = "full", overwrite = TRUE)
  expect_identical(unique(open_frame_h5(frame_path)$observations$data$run_id), "full")
})

test_that("labeled volume sets become a label x voxel frame", {
  fx <- .migration_fixture()
  dims <- dim(fx$mask)
  data <- array(seq_len(prod(dims) * 3), c(dims, 3))
  vec <- DenseNeuroVec(data, NeuroSpace(c(dims, 3), trans = fx$trans))
  legacy <- tempfile(fileext = ".h5")
  frame_path <- tempfile(fileext = ".h5")
  on.exit(unlink(c(legacy, frame_path)), add = TRUE)
  write_labeled_vec(vec, fx$mask, c("Faces", "Words", "Motion"), file = legacy)$close_all()

  migrate_labeled_h5(legacy, frame_path)
  frame <- open_frame_h5(frame_path)
  expect_identical(fmridataset::observation_ids(frame), c("Faces", "Words", "Motion"))
  support <- which(fx$mask@.Data)
  values <- fmridataset::collect_assay(frame, "value")
  expect_equal(values[2, ], as.numeric(data[, , , 2])[support])
  expect_identical(fmridataset::space(frame)$support, support)
})

test_that("the example parcellated fixture converts one kind at a time", {
  path <- suppressWarnings(create_minimal_h5_for_H5ParcellatedMultiScan())
  on.exit(unlink(path), add = TRUE)
  expect_error(
    suppressWarnings(migrate_parcellated_h5(path, tempfile(fileext = ".h5"))),
    "mixes full and summary"
  )
  for (kind in c("full", "summary")) {
    dest <- tempfile(fileext = ".h5")
    suppressWarnings(migrate_parcellated_h5(path, dest, kind = kind))
    values <- fmridataset::collect_assay(open_frame_h5(dest))
    expect_true(all(is.finite(values)))
    expect_gt(nrow(values), 0)
    unlink(dest)
  }
})
