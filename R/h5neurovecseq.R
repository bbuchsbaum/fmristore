#' H5NeuroVecSeq constructor
#'
#' @description
#' Load a NeuroVecSeq stored in HDF5 format. Each scan is stored
#' under `/scans/<name>/data` and will be returned as an `H5NeuroVec`.
#'
#' @param file Path to the HDF5 file created by `neurovecseq_to_h5()`.
#' @return An object of class `H5NeuroVecSeq`.
#' @export
H5NeuroVecSeq <- function(file) {
  assertthat::assert_that(file.exists(file))
  h5obj <- hdf5r::H5File$new(file, mode = "r")
  rtype <- try(hdf5r::h5attr(h5obj, "rtype"), silent = TRUE)
  if (!is.character(rtype) || rtype != "NeuroVecSeq") {
    stop("Invalid NeuroVecSeq HDF5 file: ", file)
  }
  scans_grp <- h5obj[["scans"]]
  snames <- names(scans_grp)
  vecs <- lapply(snames, function(nm) {
    H5NeuroVec(file, dataset_name = file.path("scans", nm, "data"))
  })
  names(vecs) <- snames
  new("H5NeuroVecSeq", obj = h5obj, vecs = vecs)
}

#' @rdname n_scans-methods
#' @export
setMethod("n_scans", "H5NeuroVecSeq", function(x) length(x@vecs))

#' @rdname scan_names-methods
#' @export
setMethod("scan_names", "H5NeuroVecSeq", function(x) names(x@vecs))

#' @export
setMethod("[[", signature(x = "H5NeuroVecSeq", i = "ANY"),
  function(x, i, j, ...) {
    x@vecs[[i]]
  })

#' @rdname h5file-methods
#' @export
setMethod("h5file", "H5NeuroVecSeq", function(x) {
  if (is.null(x@obj)) {
    stop("H5File object is NULL")
  }
  x@obj$filename
})

#' @rdname close
#' @export
setMethod("close", "H5NeuroVecSeq", function(con, ...) {
  if (!is.null(con@obj)) safe_h5_close(con@obj)
  invisible(NULL)
})
