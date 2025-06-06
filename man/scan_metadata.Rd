% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_generic.R
\name{scan_metadata}
\alias{scan_metadata}
\title{Get scan metadata}
\usage{
scan_metadata(x)
}
\arguments{
\item{x}{The object from which to retrieve the scan metadata}
}
\value{
A list (or other structure) containing metadata for each scan
}
\description{
This generic returns any available metadata associated with each scan.
}
\examples{
\dontrun{
if (requireNamespace("neuroim2", quietly = TRUE) &&
    requireNamespace("hdf5r", quietly = TRUE) &&
    exists("H5ClusterExperiment", where = "package:fmristore") &&
    exists("scan_metadata", where = "package:fmristore") &&
    !is.null(fmristore:::create_minimal_h5_for_H5ClusterExperiment)) {

  temp_exp_file <- NULL
  exp_obj <- NULL
  tryCatch({
    # Create a minimal H5ClusterExperiment
    temp_exp_file <- fmristore:::create_minimal_h5_for_H5ClusterExperiment()
    exp_obj <- fmristore::H5ClusterExperiment(file_path = temp_exp_file)

    # Get the scan metadata
    s_meta <- scan_metadata(exp_obj)
    print(s_meta)
    # The helper currently doesn't add rich scan_metadata,
    # so this might be an empty list or list of NULLs by default.
    # length(s_meta) == n_scans(exp_obj) # This should hold TRUE

  }, error = function(e) {
    message("scan_metadata example failed: ", e$message)
  }, finally = {
    if (!is.null(exp_obj)) try(close(exp_obj), silent = TRUE)
    if (!is.null(temp_exp_file) && file.exists(temp_exp_file)) {
      unlink(temp_exp_file)
    }
  })
} else {
  message("Skipping scan_metadata example: dependencies or helper not available.")
}
}

}
