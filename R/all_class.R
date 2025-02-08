#' @import methods
#' @import Matrix
#' @importFrom hdf5r H5File
#' @import neuroim2
NULL

# S4 needs to know 'H5File' is an S3 class from hdf5r
setOldClass("H5File")

#' H5NeuroVol Class
#'
#' @description
#' A class representing a three-dimensional brain image backed by an HDF5 dataset.
#' \code{H5NeuroVol} objects provide efficient storage and access for large-scale
#' brain imaging data in the HDF5 format.
#'
#' @slot h5obj An instance of class \code{H5File} from the \pkg{hdf5r} package,
#'   representing the underlying HDF5 file containing the brain image data.
#'
#' @details
#' The \code{H5NeuroVol} class inherits a \code{space} slot from
#' \code{\link[neuroim2]{NeuroVol-class}}, which specifies the spatial domain
#' (dimensions, orientation, etc.). Data I/O is performed by reading/writing
#' subsets of the HDF5 dataset, allowing efficient handling of large 3D volumes.
#'
#' @section Inheritance:
#' \code{H5NeuroVol} inherits from:
#' \itemize{
#'   \item \code{\link[neuroim2]{NeuroVol-class}}: Base class for 3D brain volumes
#'   \item \code{\link[neuroim2]{ArrayLike3D-class}}: Interface for 3D array-like operations
#' }
#'
#' @seealso
#' \code{\link[neuroim2]{NeuroVol-class}} for the base 3D brain volume class.
#' \code{\link[hdf5r]{H5File}} for details on HDF5 file objects.
#'
#' @examples
#' \dontrun{
#' library(hdf5r)
#' library(neuroim2)
#'
#' # Create or open an HDF5 file
#' h5_file <- H5File$new("brain_volume.h5", mode = "w")
#'
#' # Create a dataset in the HDF5 file
#' h5_file$create_dataset("brain_data",
#'                        dims = c(64, 64, 64),
#'                        dtype = h5types$H5T_NATIVE_FLOAT)
#'
#' # Create an H5NeuroVol object (requires a NeuroSpace)
#' h5_vol <- H5NeuroVol(h5obj = h5_file, space = NeuroSpace(dim = c(64, 64, 64)))
#'
#' # Access a subset of the data
#' subset <- h5_vol[1:10, 1:10, 1:10]
#'
#' # Close the HDF5 file when done
#' h5_file$close()
#' }
#'
#' @export
#' @rdname H5NeuroVol-class
#' @importClassesFrom neuroim2 ArrayLike3D
#' @importMethodsFrom neuroim2 [
#' @importMethodsFrom neuroim2 linear_access
setClass("H5NeuroVol",
         slots = c(
           h5obj = "H5File"  # underlying HDF5 file handle
         ),
         contains = c("NeuroVol", "ArrayLike3D"))

#' H5NeuroVec Class
#'
#' @description
#' A class representing a four-dimensional brain image backed by an HDF5 file.
#' \code{H5NeuroVec} objects provide efficient storage and access for large-scale
#' 4D neuroimaging data using the HDF5 file format.
#'
#' @slot obj An instance of class \code{H5File} from the \pkg{hdf5r} package,
#'   representing the underlying HDF5 file containing the 4D brain image data.
#'
#' @details
#' \code{H5NeuroVec} inherits a \code{space} slot from
#' \code{\link[neuroim2]{NeuroVec-class}}, defining the 4D dimensions
#' (e.g., \code{x}, \code{y}, \code{z}, \code{time}). Data are stored in an HDF5
#' dataset and accessed on demand.
#'
#' @section Inheritance:
#' \code{H5NeuroVec} inherits from:
#' \itemize{
#'   \item \code{\link[neuroim2]{NeuroVec-class}}: Base class for 4D brain images
#'   \item \code{\link[neuroim2]{ArrayLike4D-class}}: Interface for 4D array-like operations
#' }
#'
#' @seealso
#' \code{\link[neuroim2]{NeuroVec-class}} for the base 4D brain image class.
#' \code{\link{H5NeuroVol-class}} for the 3D counterpart.
#' \code{\link[hdf5r]{H5File}} for details on HDF5 file objects.
#'
#' @examples
#' \dontrun{
#' library(hdf5r)
#' library(neuroim2)
#'
#' # Create or open an HDF5 file
#' h5_file <- H5File$new("brain_timeseries.h5", mode = "w")
#'
#' # Create a dataset in the HDF5 file
#' h5_file$create_dataset("fmri_data",
#'                        dims = c(64, 64, 32, 100),
#'                        dtype = h5types$H5T_NATIVE_FLOAT)
#'
#' # Create an H5NeuroVec object
#' h5_vec <- H5NeuroVec(obj = h5_file,
#'                      space = NeuroSpace(dim = c(64, 64, 32, 100)))
#'
#' # Access a subset of the data
#' subset <- h5_vec[1:10, 1:10, 1:10, 1:5]
#'
#' # Close the HDF5 file when done
#' h5_file$close()
#' }
#'
#' @export
#' @rdname H5NeuroVec-class
setClass("H5NeuroVec",
         slots = c(
           obj = "H5File"  # underlying HDF5 file handle
         ),
         contains = c("NeuroVec", "ArrayLike4D"))

#' H5Format Class
#'
#' @description
#' This class represents the HDF5 (Hierarchical Data Format version 5) file format.
#' It extends \code{FileFormat-class} with HDF5-specific attributes.
#'
#' @seealso \code{\link[neuroim2]{FileFormat-class}}
#'
#' @keywords internal
#' @noRd
setClass("H5Format",
         contains = c("FileFormat"))

#' H5NeuroVecSource Class
#'
#' @description
#' A class used internally to produce an \code{\linkS4class{H5NeuroVec}} from an HDF5 file.
#'
#' @slot file_name A \code{character} string specifying the name of the HDF5 file.
#'
#' @seealso \code{\link{H5NeuroVec-class}}
#'
#' @keywords internal
#' @noRd
setClass("H5NeuroVecSource",
         representation(file_name = "character"))

#' LatentNeuroVec Class
#'
#' @description
#' A class that represents a 4-dimensional neuroimaging array using a latent space
#' decomposition. It stores the data as a set of basis functions (dictionary) and
#' a corresponding set of loadings (coefficients), enabling efficient representation
#' and manipulation of high-dimensional data.
#'
#' @slot basis A \code{Matrix} object where each column represents a basis vector
#'   in the latent space.
#' @slot loadings A \code{Matrix} object (often sparse) containing the coefficients
#'   for each basis vector across the spatial dimensions.
#' @slot offset A \code{numeric} vector representing a constant offset term for
#'   each voxel or spatial location.
#'
#' @details
#' \code{LatentNeuroVec} inherits from \code{\link[neuroim2]{NeuroVec-class}}
#' and \code{\link[neuroim2]{AbstractSparseNeuroVec-class}}. The original 4D data
#' can be reconstructed as:
#' \deqn{data[x,y,z,t] = \sum_i \bigl(basis[t,i] \times loadings[i,v]\bigr) + offset[v]}.
#'
#' This approach is especially useful for large datasets where storing the full
#' 4D array is expensive.
#'
#' @section Inheritance:
#' \code{LatentNeuroVec} inherits from:
#' \itemize{
#'   \item \code{\link[neuroim2]{NeuroVec-class}}
#'   \item \code{\link[neuroim2]{AbstractSparseNeuroVec-class}}
#' }
#'
#' @seealso
#' \code{\link[neuroim2]{NeuroVec-class}},
#' \code{\link[neuroim2]{AbstractSparseNeuroVec-class}}.
#'
#' @examples
#' \dontrun{
#' library(Matrix)
#' library(neuroim2)
#'
#' # Create a LatentNeuroVec object
#' basis <- Matrix(rnorm(100 * 10), nrow = 100, ncol = 10)
#' loadings <- Matrix(rnorm(10 * 1000), nrow = 10, ncol = 1000, sparse = TRUE)
#' offset <- rnorm(1000)
#'
#' latent_vec <- new("LatentNeuroVec",
#'                   basis = basis,
#'                   loadings = loadings,
#'                   offset = offset,
#'                   space = NeuroSpace(dim = c(10, 10, 10, 100)))
#'
#' # Access the basis
#' basis_functions <- latent_vec@basis
#'
#' # Reconstruct a single time point
#' time_point_1 <- latent_vec@loadings %*% latent_vec@basis[1,] + latent_vec@offset
#' }
#'
#' @export
#' @rdname LatentNeuroVec-class
setClass("LatentNeuroVec",
         slots = c(
           basis = "Matrix",
           loadings = "Matrix",
           offset = "numeric"
         ),
         contains = c("NeuroVec", "AbstractSparseNeuroVec"))

#' LatentNeuroVecSource Class
#'
#' @description
#' A class used internally to produce a \code{\linkS4class{LatentNeuroVec}} instance.
#'
#' @slot file_name A \code{character} string specifying the file name.
#'
#' @seealso \code{\link{LatentNeuroVec-class}}
#'
#' @keywords internal
#' @noRd
setClass("LatentNeuroVecSource",
         representation(file_name="character"))

#' H5ClusteredVec Class
#'
#' @description
#' A class representing a single clustered time series scan backed by HDF5.
#' Provides a NeuroVec-like interface while maintaining HDF5 backing.
#'
#' @slot obj The \code{H5File} object referencing the open file.
#' @slot scan_name A \code{character} label identifying this scan.
#' @slot mask A \code{LogicalNeuroVol} defining the brain mask.
#' @slot clusters A \code{ClusteredNeuroVol} containing cluster assignments.
#'
#' @section Inheritance:
#' Inherits from \code{H5NeuroVec}, which is in turn a \code{NeuroVec}.
#'
#' @importFrom hdf5r H5File
#' @export
setClass("H5ClusteredVec",
         slots = c(
           obj = "H5File",
           scan_name = "character",
           mask = "LogicalNeuroVol",
           clusters = "ClusteredNeuroVol"
         ),
         contains = c("H5NeuroVec"))

#' H5ClusteredVecSeq Class
#'
#' @description
#' A class representing a sequence of clustered time series scans backed by an HDF5 file.
#' It provides a \code{NeuroVecSeq}-like interface while maintaining HDF5 backing.
#'
#' @slot obj The \code{H5File} object referencing the open file.
#' @slot scan_names A \code{character} vector of scan IDs.
#' @slot mask A \code{LogicalNeuroVol} defining the brain mask.
#' @slot clusters A \code{ClusteredNeuroVol} containing cluster assignments.
#' @slot scan_metadata A \code{list} of metadata, one entry per scan.
#' @slot cluster_metadata A \code{data.frame} containing metadata for each cluster.
#'
#' @importFrom hdf5r H5File
#' @export
setClass("H5ClusteredVecSeq",
         slots = c(
           obj = "H5File",
           scan_names = "character",
           mask = "LogicalNeuroVol",
           clusters = "ClusteredNeuroVol",
           scan_metadata = "list",
           cluster_metadata = "data.frame"
         ),
         contains = c("NeuroVecSeq"))

#' H5ReducedClusteredVec Class
#'
#' @description
#' Represents a single scan's *summaries*—i.e., one time-series per cluster
#' in a 2D dataset [nTime, nClusters]. Each column is a cluster, each row
#' is a time point.
#'
#' @details
#' Typically stored in
#' \code{/scans/<scan_name>/clusters_summary/summary_data}, with shape
#' [nTime, nClusters]. A \code{cluster_names} vector can label each column.
#'
#' @slot obj The \code{H5File} object referencing the open file on disk.
#' @slot scan_name The scan identifier (e.g., a subgroup name).
#' @slot mask A \code{LogicalNeuroVol} with 3D shape (optional, for reference).
#' @slot clusters A \code{ClusteredNeuroVol} describing the cluster partition.
#' @slot cluster_names A \code{character} vector naming each column.
#' @slot cluster_ids An \code{integer} vector (same length) if numeric IDs are tracked.
#' @slot n_time An \code{numeric} specifying how many time points.
#'
#' @section Inheritance:
#' Inherits from \code{H5NeuroVec}, so it is considered a 4D object with the 4th dimension
#' effectively replaced by cluster columns. Typically used for summarized timeseries.
#'
#' @seealso \code{H5ReducedClusteredVecSeq}
#'
#' @importFrom hdf5r H5File
#' @export
setClass("H5ReducedClusteredVec",
         slots = c(
           obj = "H5File",
           scan_name = "character",
           mask = "LogicalNeuroVol",
           clusters = "ClusteredNeuroVol",
           cluster_names = "character",
           cluster_ids = "integer",
           n_time = "numeric"
         ),
         contains = c("H5NeuroVec"))



#' H5ReducedClusteredVecSeq Class
#'
#' @description
#' Represents a multi-scan collection of reduced (summary) timeseries in an HDF5 file.
#' Each scan has [nTime, nClusters] data at /scans/<scan_name>/clusters_summary/summary_data.
#'
#' @slot obj The \code{H5File} object referencing the open file.
#' @slot scan_names A \code{character} vector of scan IDs.
#' @slot mask A \code{LogicalNeuroVol} for 3D geometry reference.
#' @slot clusters A \code{ClusteredNeuroVol} describing the global cluster partition.
#' @slot scan_metadata A \code{list} of metadata for each scan.
#' @slot cluster_metadata A \code{data.frame} or \code{list} describing each cluster.
#' @slot cluster_names A \code{character} vector labeling columns (common across scans).
#' @slot cluster_ids An \code{integer} vector of cluster IDs, if needed.
#' @slot n_time A \code{numeric} vector specifying time lengths per scan.
#'
#' @section Inheritance:
#' Inherits from \code{NeuroVecSeq}, since it is effectively multiple (reduced) timeseries.
#'
#' @seealso \code{H5ReducedClusteredVec}
#'
#' @importFrom hdf5r H5File
#' @export
setClass("H5ReducedClusteredVecSeq",
         slots = c(
           obj = "H5File",
           scan_names = "character",
           mask = "LogicalNeuroVol",
           clusters = "ClusteredNeuroVol",
           scan_metadata = "list",
           cluster_metadata = "data.frame",
           cluster_names = "character",
           cluster_ids = "integer",
           n_time = "numeric"
         ),
         contains = c("NeuroVecSeq"))

#' LabeledVolumeSet Class
#'
#' @description
#' A class representing a multi-volume dataset stored in HDF5, where each volume
#' is labeled (similar to a named 4th dimension). We store:
#' \itemize{
#'   \item a 3D mask
#'   \item a set of label strings
#'   \item for each label, data only at the mask's nonzero entries
#' }
#'
#' This extends \code{NeuroVec}, so it is logically a 4D object with dimension
#' \code{[X, Y, Z, #labels]}.
#'
#' @slot obj An \code{H5File} reference (the file handle).
#' @slot mask A \code{LogicalNeuroVol} of shape [X, Y, Z].
#' @slot labels A \code{character} vector for the volume labels.
#' @slot load_env An \code{environment} storing references for lazy loading.
#'
#' @seealso \code{\link[neuroim2]{NeuroVec-class}}
#' @export
setClass("LabeledVolumeSet",
         slots = c(
           obj      = "H5File",          # pointer to open file
           mask     = "LogicalNeuroVol", # 3D mask
           labels   = "character",       # volume (sub-volume) names
           load_env = "environment"      # environment for lazy loading
         ),
         contains = c("NeuroVec"))  # extends NeuroVec
