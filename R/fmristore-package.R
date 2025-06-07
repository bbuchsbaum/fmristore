#' fmristore: Efficient Storage of fMRI Data
#'
#' The fmristore package provides efficient storage and retrieval of functional
#' magnetic resonance imaging (fMRI) data using HDF5 format. It offers S4 classes
#' and methods for working with dense and sparse representations of 4D neuroimaging
#' data, latent space decompositions, and spatially clustered voxel time series.
#'
#' @section Main Features:
#' \itemize{
#'   \item HDF5-backed storage for memory-efficient access to large neuroimaging datasets
#'   \item Support for both dense and sparse 4D fMRI data representations
#'   \item Latent space decompositions using basis functions and loadings
#'   \item Spatial clustering approaches for organizing voxel time series
#'   \item Integration with the neuroim2 package for neuroimaging data structures
#' }
#'
#' @section Key Classes:
#' \itemize{
#'   \item \code{\link{H5NeuroVol}} and \code{\link{H5NeuroVec}}: HDF5-backed 3D/4D volumes
#'   \item \code{\link{LatentNeuroVec}}: Latent representations with spatial basis and temporal loadings
#'   \item \code{\link{H5ClusterExperiment}}: Container for multiple runs with clustered organization
#'   \item \code{LabeledVolumeSet}: Collection of labeled brain regions
#' }
#'
#' @section Getting Started:
#' For an introduction to the package, see the vignettes:
#' \itemize{
#'   \item \code{vignette("H5Neuro", package = "fmristore")}: Working with HDF5-backed neuroimaging data
#'   \item \code{vignette("LabeledVolumeSet", package = "fmristore")}: Managing labeled brain regions
#' }
#'
#' @keywords internal
#' @aliases fmristore-package fmristore
#' @author Bradley Buchsbaum \email{brad.buchsbaum@@gmail.com} [aut, cre]
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
