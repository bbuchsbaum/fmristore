# test_fmristore_extras.R

library(testthat)
library(neuroim2)
library(fmristore)
library(Matrix)

test_that("LabeledVolumeSet partial usage + memoise=TRUE", {
  # create a small 4D data: (4,4,2,3) => 3 volumes
  arr <- array(rnorm(4 * 4 * 2 * 3), dim = c(4, 4, 2, 3))
  spc <- NeuroSpace(dim = c(4, 4, 2, 3))
  vec <- NeuroVec(arr, spc)

  # mask => pick random subset
  mask_arr <- array(sample(c(TRUE, FALSE), 4 * 4 * 2, replace = TRUE), dim = c(4, 4, 2))
  mask_vol <- LogicalNeuroVol(mask_arr, NeuroSpace(c(4, 4, 2)))

  labels <- c("volA", "volB", "volC")

  tmpfile <- tempfile(fileext = ".h5")
  on.exit(unlink(tmpfile))

  # write as labeled
  write_labeled_vec(vec = vec, mask = mask_vol, labels = labels, file = tmpfile, compression = 1)

  # read it with memoise=TRUE
  lvs <- read_labeled_vec(tmpfile)
  expect_s4_class(lvs, "LabeledVolumeSet")

  # check subsetting => [i,j,k,l]
  sub_4d <- lvs[1:2, 2:3, 1:2, 1:2]
  expect_equal(dim(sub_4d), c(2, 2, 2, 2))
  # compare manually
  gold <- arr[1:2, 2:3, 1:2, 1:2] * rep(mask_arr[1:2, 2:3, 1:2], 2)
  expect_equal(sub_4d, gold, tolerance = 1e-7)

  # check linear_access
  lin_idx <- c(1, 4, 10, 16, 20, 30)
  lv_vals <- linear_access(lvs, lin_idx)

  # construct "expected_full"
  # each vol => arr[,,,vol] * mask
  arr_masked <- arr
  for (v in seq_len(3)) {
    arr_masked[, , , v] <- arr_masked[, , , v] * mask_arr
  }
  gold_lin <- arr_masked[lin_idx]
  expect_equal(lv_vals, gold_lin, tolerance = 1e-7)

  # check single volume => [[
  vC <- lvs[[3]]
  expect_s4_class(vC, "DenseNeuroVol")
  # verify it matches arr[,,,3]*mask
  expect_equal(as.array(vC@.Data), arr[, , , 3] * mask_arr, tolerance = 1e-7)
})
