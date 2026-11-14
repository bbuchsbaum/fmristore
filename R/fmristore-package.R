#' fmristore: Efficient Storage of fMRI Data
#'
#' fmristore persists \code{fmridataset::fmri_frame} objects in HDF5: atomic
#' single-file snapshots (\code{\link{write_frame_h5}},
#' \code{\link{open_frame_h5}}), append-only sharded stores
#' (\code{\link{write_sharded_frame_h5}}), recovery and migration
#' (\code{\link{recover_frame_store}}, \code{\link{migrate_frame_h5}}), and
#' lazy HDF5 array sources (\code{\link{h5_array_source}}).
#'
#' @section Legacy classes:
#' \code{H5ParcellatedMultiScan}, \code{H5ParcellatedScan},
#' \code{H5ParcellatedScanSummary} and \code{LabeledVolumeSet} remain only
#' for packages that have not yet moved to frames. They are frozen and will be
#' removed; see \code{design/deficiency-remediation-plan.md} in the source
#' repository.
#'
#' @keywords internal
#' @aliases fmristore-package fmristore
#' @author Bradley Buchsbaum \email{brad.buchsbaum@@gmail.com} [aut, cre]
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
