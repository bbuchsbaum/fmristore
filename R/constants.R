#' HDF5 Path Constants
#'
#' @description
#' Central definition of HDF5 file layout constants used throughout the package.
#' This provides a single source of truth for the HDF5 file structure.
#'
#' @format A named list of character strings containing HDF5 path templates.
#' @return A named list with HDF5 path constants.
#' @export
#' @examples
#' # Access constants
#' H5_PATHS$SCANS_GROUP # "/scans"
#' sprintf(H5_PATHS$CLUSTER_DSET_TPL, "run1", 5) # "/scans/run1/clusters/cluster_5"
H5_PATHS <- list(
  # Top-level groups
  SCANS_GROUP = "/scans",
  BASIS_GROUP = "/basis",
  CLUSTER_MAP = "/cluster_map",
  MASK = "/mask",
  HEADER_GROUP = "/header",
  SPACE_GROUP = "/space",
  VOXEL_COORDS = "/voxel_coords",

  # Scan-level path templates (use with sprintf)
  SCAN_GROUP_TPL = "/scans/%s",
  CLUSTERS_GROUP_TPL = "/scans/%s/clusters",
  SUMMARY_GROUP_TPL = "/scans/%s/clusters_summary",
  EMBEDDING_TPL = "/scans/%s/embedding",
  METADATA_GROUP_TPL = "/scans/%s/metadata",

  # Dataset templates (use with sprintf)
  CLUSTER_DSET_TPL = "/scans/%s/clusters/cluster_%d",
  SUMMARY_DSET_TPL = "/scans/%s/clusters_summary/%s",
  N_TIME_ATTR = "n_time",

  # Basis and latent structure
  BASIS_MATRIX = "/basis/basis_matrix",
  BASIS_OFFSET = "/basis/offset",
  BASIS_LABELS = "/basis/labels",

  # Metadata paths
  CLUSTER_METADATA = "/cluster_metadata",
  SCAN_METADATA_TPL = "/scans/%s/metadata",
  N_TIME_METADATA_TPL = "/scans/%s/metadata/n_time"
)

#' HDF5 Attribute Constants
#'
#' @description
#' Named constants for HDF5 attribute names used throughout the package.
#'
#' @format A named list of character strings.
#' @return A named list with HDF5 attribute name constants.
#' @export
#' @examples
#' H5_ATTRS$N_TIME
H5_ATTRS <- list(
  N_TIME = "n_time",
  COMPRESS = "compress",
  CLUSTER_NAMES = "cluster_names",
  CLUSTER_IDS = "cluster_ids",
  SUMMARY_ONLY = "summary_only",
  N_VOXELS = "n_voxels",
  N_CLUSTERS = "n_clusters"
)

#' HDF5 Dataset Name Constants
#'
#' @description
#' Named constants for HDF5 dataset names used throughout the package.
#'
#' @format A named list of character strings.
#' @return A named list with HDF5 dataset name constants.
#' @export
#' @examples
#' H5_DSETS$DATA
H5_DSETS <- list(
  SUMMARY_DATA = "summary_data",
  ELEMENTS = "elements",
  DATA = "data"
)
