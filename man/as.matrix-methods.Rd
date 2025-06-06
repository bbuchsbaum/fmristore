% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_generic.R, R/cluster_array.R
\name{as.matrix}
\alias{as.matrix}
\alias{as.matrix,H5ClusterRunSummary-method}
\title{Convert to Matrix}
\usage{
as.matrix(x, ...)

\S4method{as.matrix}{H5ClusterRunSummary}(x)
}
\arguments{
\item{x}{An object for which a `as.matrix` method is defined.}

\item{...}{Additional arguments passed to specific as.matrix methods.}
}
\description{
These methods convert `fmristore` specific objects to matrices,
implementing the S4 generic function \code{\link[base]{as.matrix}} from the `base` package.
}
\seealso{
Other H5Cluster: 
\code{\link{$,H5ClusterExperiment-method}},
\code{\link{H5ClusterExperiment}},
\code{\link{H5ClusterExperiment-class}},
\code{\link{H5ClusterRun-class}},
\code{\link{H5ClusterRunSummary-class}},
\code{\link{H5ClusteredArray-class}},
\code{\link{[,H5ClusterRun,ANY,ANY,ANY-method}},
\code{\link{[,H5ClusterRun,ANY,missing,ANY-method}},
\code{\link{as.data.frame}()},
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
