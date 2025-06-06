% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_generic.R, R/cluster_array.R,
%   R/cluster_experiment.R
\name{clusters}
\alias{clusters}
\alias{clusters,H5ClusteredArray-method}
\alias{clusters,H5ClusterExperiment-method}
\title{Get the cluster map object}
\usage{
clusters(x, ...)

\S4method{clusters}{H5ClusteredArray}(x)

\S4method{clusters}{H5ClusterExperiment}(x)
}
\arguments{
\item{x}{H5ClusterExperiment object}

\item{...}{Additional arguments}
}
\value{
The clusters object (e.g., a ClusteredNeuroVol)

The clusters object from the first run.
}
\description{
Get the cluster map object

Get the clusters object via generic
}
\examples{
\dontrun{
# For H5ClusterExperiment:
if (!is.null(fmristore:::create_minimal_h5_for_H5ClusterExperiment)) {
  temp_exp_file <- NULL
  exp_obj <- NULL
  tryCatch({
    temp_exp_file <- fmristore:::create_minimal_h5_for_H5ClusterExperiment(
      master_mask_dims = c(4L,4L,3L),
      num_master_clusters = 2L
    )
    exp_obj <- fmristore::H5ClusterExperiment(file_path = temp_exp_file)

    # Get the master cluster map from the experiment
    cluster_vol <- clusters(exp_obj)
    print(cluster_vol)
    # if (requireNamespace("neuroim2", quietly=TRUE)) print(is(cluster_vol, "ClusteredNeuroVol"))

    # Individual runs also have cluster information, potentially accessible via their own methods
    # run1 <- runs(exp_obj)[["Run1_Full"]]
    # run1_clusters <- clusters(run1) # Assuming a method for H5ClusterRun
    # print(run1_clusters)

  }, error = function(e) {
    message("clusters example for H5ClusterExperiment failed: ", e$message)
  }, finally = {
    if (!is.null(exp_obj)) try(close(exp_obj), silent = TRUE)
    if (!is.null(temp_exp_file) && file.exists(temp_exp_file)) {
      unlink(temp_exp_file)
    }
  })
} else {
  message("Skipping clusters example for H5ClusterExperiment: helper not available.")
}
}

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
\code{\link{as.matrix}()},
\code{\link{close}()},
\code{\link{cluster_metadata,H5ClusterExperiment-method}},
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
