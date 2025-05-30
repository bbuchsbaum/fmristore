% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/h5neurovec.R, R/h5neurovol.R, R/labeled_vec.R, R/cluster_experiment.R
\name{close,H5NeuroVec-method}
\alias{close,H5NeuroVec-method}
\alias{close,H5NeuroVol-method}
\alias{close,LabeledVolumeSet-method}
\alias{close,H5ClusterExperiment-method}
\title{Close the HDF5 file associated with an H5NeuroVec}
\usage{
\S4method{close}{H5NeuroVec}(con, ...)

\S4method{close}{H5NeuroVol}(con, ...)

\S4method{close}{LabeledVolumeSet}(con, ...)
\S4method{close}{H5ClusterExperiment}(con, ...)
}
\arguments{
\item{con}{A \\code{LabeledVolumeSet} object.}

\item{...}{Additional arguments (ignored).}
}
\value{
Invisibly returns \code{NULL}.

Invisibly returns \code{NULL}.

Invisibly returns \\code{NULL}.
Invisibly returns \code{NULL}.
}
\description{
This method manually closes the HDF5 file handle stored within the
H5NeuroVec object. It uses the \code{safe_h5_close} helper.

This method manually closes the HDF5 file handle stored within the
H5NeuroVol object. It uses the \code{safe_h5_close} helper.

This method manually closes the HDF5 file handle stored within the
LabeledVolumeSet object. It uses the \\code{safe_h5_close} helper to
ensure the handle is valid before attempting to close. After closing,
the internal handle reference is nulled to prevent accidental reuse.
This method manually closes the shared HDF5 file handle stored within the H5ClusterExperiment object. It uses the \code{safe_h5_close} helper.
}
\details{
**Important:** If this \code{LabeledVolumeSet} object was created from
a file path using \code{\link{read_labeled_vec}}, the user is responsible
for calling this \code{close} method when finished with the object to release
the file handle. Failure to do so will leave the file open until the R session ends.
If the object was created using an existing \code{H5File} handle, closing
remains the responsibility of the code that originally opened the handle.
}
