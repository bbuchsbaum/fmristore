% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_generic.R
\name{scan_names}
\alias{scan_names}
\title{Get the scan names}
\usage{
scan_names(x)
}
\arguments{
\item{x}{The object from which to retrieve the scan names}
}
\value{
A character vector of scan names
}
\description{
This generic returns a character vector of scan names or labels.
}
\examples{
\dontrun{
if (requireNamespace("neuroim2", quietly = TRUE) &&
    requireNamespace("hdf5r", quietly = TRUE) &&
    exists("H5ClusterExperiment", where = "package:fmristore") &&
    exists("scan_names", where = "package:fmristore") &&
    !is.null(fmristore:::create_minimal_h5_for_H5ClusterExperiment)) {

  temp_exp_file <- NULL
  exp_obj <- NULL
  tryCatch({
    # Create a minimal H5ClusterExperiment
    temp_exp_file <- fmristore:::create_minimal_h5_for_H5ClusterExperiment()
    exp_obj <- fmristore::H5ClusterExperiment(file_path = temp_exp_file)

    # Get the scan names
    s_names <- scan_names(exp_obj)
    print(s_names) # Should be c("Run1_Full", "Run2_Summary") or similar

  }, error = function(e) {
    message("scan_names example failed: ", e$message)
  }, finally = {
    if (!is.null(exp_obj)) try(close(exp_obj), silent = TRUE)
    if (!is.null(temp_exp_file) && file.exists(temp_exp_file)) {
      unlink(temp_exp_file)
    }
  })
} else {
  message("Skipping scan_names example: dependencies or helper not available.")
}
}

}
