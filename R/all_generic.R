#' Get the number of scans
#'
#' This generic returns the number of scans in an object (e.g. a sequence of fMRI scans).
#'
#' @param x The object from which to retrieve the number of scans
#' @return An integer representing the number of scans
#' @export
setGeneric("n_scans", function(x) standardGeneric("n_scans"))

#' Get the scan names
#'
#' This generic returns a character vector of scan names or labels.
#'
#' @param x The object from which to retrieve the scan names
#' @return A character vector of scan names
#' @export
setGeneric("scan_names", function(x) standardGeneric("scan_names"))

#' Get scan metadata
#'
#' This generic returns any available metadata associated with each scan.
#'
#' @param x The object from which to retrieve the scan metadata
#' @return A list (or other structure) containing metadata for each scan
#' @export
setGeneric("scan_metadata", function(x) standardGeneric("scan_metadata"))

#' Get cluster metadata
#'
#' This generic returns any available metadata associated with clusters in an object.
#'
#' @param x The object from which to retrieve the cluster metadata
#' @return A data frame or other structure containing metadata for each cluster
#' @export
setGeneric("cluster_metadata", function(x) standardGeneric("cluster_metadata"))
