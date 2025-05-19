#' @rdname as_h5-methods
#' @export
setMethod(
  f = "as_h5",
  signature = signature(object = "NeuroVec"),
  definition = function(object, file = NULL, data_type = "FLOAT",
                         chunk_dim = c(4, 4, 4, dim(object)[4]),
                         compression = 6) {
    to_nih5_vec(object, file_name = file, data_type = data_type,
                chunk_dim = chunk_dim, compression = compression)
  }
)

#' @rdname as_h5-methods
#' @export
setMethod(
  f = "as_h5",
  signature = signature(object = "LatentNeuroVec"),
  definition = function(object, file = NULL, data_type = "FLOAT",
                         compression = 6) {
    to_h5_latentvec(object, file_name = file, data_type = data_type,
                    compression = compression)
  }
)

#' @rdname as_h5-methods
#' @export
setMethod(
  f = "as_h5",
  signature = signature(object = "LabeledVolumeSet"),
  definition = function(object, file, mask, labels, compression = 4,
                         dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE,
                         chunk_size = 1024, header_values = list()) {
    write_labeled_vec(vec = object, mask = mask, labels = labels, file = file,
                      compression = compression, dtype = dtype,
                      chunk_size = chunk_size, header_values = header_values)
  }
)

#' @rdname as_h5-methods
#' @export
setMethod(
  f = "as_h5",
  signature = signature(object = "list"),
  definition = function(object, file, scan_names, mask, clusters,
                         scan_metadata, cluster_metadata = NULL,
                         summary_only = FALSE, compression = 4,
                         chunk_size = 1024) {
    write_clustered_dataset(file = file, vecs = object, scan_names = scan_names,
                           mask = mask, clusters = clusters,
                           scan_metadata = scan_metadata,
                           cluster_metadata = cluster_metadata,
                           summary_only = summary_only,
                           compression = compression,
                           chunk_size = chunk_size)
  }
) 