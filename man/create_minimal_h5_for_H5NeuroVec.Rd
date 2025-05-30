% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zzz_example_helpers.R
\name{create_minimal_h5_for_H5NeuroVec}
\alias{create_minimal_h5_for_H5NeuroVec}
\title{Create a minimal HDF5 file suitable for H5NeuroVec examples}
\usage{
create_minimal_h5_for_H5NeuroVec(dims = c(3L, 3L, 2L, 5L), file_path = NULL)
}
\arguments{
\item{dims}{4D dimensions, e.g., c(3L, 3L, 2L, 5L).}

\item{file_path}{Optional: path to HDF5 file. If \code{NULL}, a temp file is created.}
}
\value{
Path to the created HDF5 file.
}
\description{
This function creates a temporary HDF5 file with the minimal structure
expected by the \code{H5NeuroVec} constructor (rtype, /space, /data).
}
\examples{
# Not typically called directly by users, used in other examples via fmristore:::
if (requireNamespace("neuroim2", quietly = TRUE) && 
    requireNamespace("hdf5r", quietly = TRUE)) {
  temp_file <- fmristore:::create_minimal_h5_for_H5NeuroVec()
  print(temp_file)
  # Example usage:
  # h5vec <- H5NeuroVec(temp_file)
  # print(h5vec)
  # close(h5vec)
  if (file.exists(temp_file)) unlink(temp_file)
}
}
\keyword{internal}
