% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cluster_experiment.R
\name{H5ClusterExperiment}
\alias{H5ClusterExperiment}
\title{Constructor for H5ClusterExperiment Objects}
\usage{
H5ClusterExperiment(
  file,
  scan_names = NULL,
  mask = NULL,
  clusters = NULL,
  scan_metadata = NULL,
  cluster_metadata = NULL,
  summary_preference = NULL,
  keep_handle_open = TRUE
)
}
\arguments{
\item{file}{Either a character string path to the HDF5 file or an
open \code{H5File} object.}

\item{scan_names}{(Optional) A character vector specifying which scans under `/scans/`
to include in the experiment. If `NULL` (default), the constructor attempts
to discover all available scan groups under `/scans/`.}

\item{mask}{(Optional) A `LogicalNeuroVol` object for the brain mask. If `NULL`,
the constructor attempts to load it from `/mask` in the HDF5 file.}

\item{clusters}{(Optional) A `ClusteredNeuroVol` object for cluster assignments.
If `NULL`, the constructor attempts to load it from `/cluster_map` and
potentially `/voxel_coords` in the HDF5 file.}

\item{scan_metadata}{(Optional) A list to override or supplement metadata read
from the HDF5 file. If provided, its length should match the number of scans.}

\item{cluster_metadata}{(Optional) A data.frame to override or supplement
cluster metadata read from `/clusters/cluster_meta` in the HDF5 file.}

\item{summary_preference}{(Optional) Character string controlling which run type to load.
If \code{NULL} (default), the constructor reads the \code{/scans@summary_only}
attribute to choose a default: \code{"require"} if \code{summary_only=TRUE},
\code{"ignore"} if \code{summary_only=FALSE}, otherwise \code{"prefer"}.
Explicit values can be "require" (only load summary runs, error if missing),
"prefer" (load summary if available, else full), or "ignore" (load full runs only).
This influences whether \code{make_run_summary} or \code{make_run_full} is called.
*Note: This parameter requires careful implementation based on HDF5 content checks.*}

\item{keep_handle_open}{(Logical) Only relevant if \code{file_source} is a path.
If \code{TRUE} (default), the HDF5 file handle is kept open within the returned
object. If \code{FALSE}, the handle is closed after reading metadata.
*Note:* For most operations, the handle needs to remain open.}
}
\value{
A new \code{H5ClusterExperiment} object.
}
\description{
Creates a new `H5ClusterExperiment` object, representing a collection of
clustered neuroimaging runs sharing a common HDF5 file, mask, and cluster map.

This function handles opening the HDF5 file (if a path is provided),
identifying available scans, and creating the appropriate run objects
(`H5ClusterRun` or `H5ClusterRunSummary`) for each scan based on
the available data within the HDF5 file structure (following the
ClusteredTimeSeriesSpec).
}
\examples{
\dontrun{
# Create temporary HDF5 file with minimal experiment structure
temp_file <- tempfile(fileext = ".h5")
exp_file <- fmristore:::create_minimal_h5_for_H5ClusterExperiment(file_path = temp_file)

# Create H5ClusterExperiment object
experiment <- H5ClusterExperiment(exp_file)

# Access scan names
print(scan_names(experiment))

# Get number of scans
print(n_scans(experiment))

# Access runs
print(names(experiment@runs))

# Extract data for specific voxels (first 5 mask indices)
voxel_data <- series_concat(experiment, mask_idx = 1:5)
print(dim(voxel_data))

# Clean up
close(experiment)
unlink(temp_file)
}

}
\seealso{
Other H5Cluster: 
\code{\link{$,H5ClusterExperiment-method}},
\code{\link{H5ClusterExperiment-class}},
\code{\link{H5ClusterRun-class}},
\code{\link{H5ClusterRunSummary-class}},
\code{\link{H5ClusteredArray-class}},
\code{\link{[,H5ClusterRun,ANY,ANY,ANY-method}},
\code{\link{[,H5ClusterRun,ANY,missing,ANY-method}},
\code{\link{as.data.frame}()},
\code{\link{as.matrix}()},
\code{\link{close}()},
\code{\link{cluster_metadata,H5ClusterExperiment-method}},
\code{\link{clusters}()},
\code{\link{dim}()},
\code{\link{h5file,H5ClusteredArray-method}},
\code{\link{linear_access-methods}},
\code{\link{make_run_full}()},
\code{\link{make_run_summary}()},
\code{\link{mask}()},
\code{\link{matrix_concat}()},
\code{\link{n_scans,H5ClusterExperiment-method}},
\code{\link{scan_metadata,H5ClusterExperiment-method}},
\code{\link{scan_names,H5ClusterExperiment-method}},
\code{\link{series}()},
\code{\link{series_concat}()},
\code{\link{show,H5ClusterRun-method}}
}
\concept{H5Cluster}
