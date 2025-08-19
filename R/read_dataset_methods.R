#' @include all_generic.R
#' @include all_class.R
#' @include h5_utils.R
#' @importFrom hdf5r H5File h5attr h5attr_names
#' @importFrom neuroim2 LogicalNeuroVol NeuroSpace
#' @importFrom utils packageVersion
NULL

# Helper Functions ----

#' Check if a path exists in an HDF5 file
#'
#' Safely checks if a path exists in an HDF5 file
#'
#' @param h5 An H5File object
#' @param path Character string path to check
#' @return Logical indicating if path exists
#' @keywords internal
exists_in_h5 <- function(h5, path) {
  tryCatch(
    {
      h5$exists(path)
    },
    error = function(e) {
      FALSE
    }
  )
}

#' Check if an attribute exists in an HDF5 object
#'
#' Safely checks if an attribute exists
#'
#' @param h5 An H5File or H5Group object
#' @param attr_name Character string attribute name
#' @return Logical indicating if attribute exists
#' @keywords internal
h5_attr_exists <- function(h5, attr_name) {
  tryCatch(
    {
      attr_name %in% h5attr_names(h5)
    },
    error = function(e) {
      FALSE
    }
  )
}

#' Get HDF5 attribute safely
#'
#' Gets an attribute value or returns NULL if not found
#'
#' @param h5 An H5File or H5Group object
#' @param attr_name Character string attribute name
#' @return Attribute value or NULL
#' @keywords internal
h5_attr_safe <- function(h5, attr_name) {
  if (h5_attr_exists(h5, attr_name)) {
    h5attr(h5, attr_name)
  } else {
    NULL
  }
}

# Type Detection ----

#' Detect the type of data in an HDF5 file
#'
#' @description
#' Examines the structure and attributes of an HDF5 file to determine
#' what type of fmristore object it contains.
#'
#' @param h5obj An H5File object or file path to an HDF5 file
#' @return Character string indicating the type, one of:
#'   - "H5NeuroVol" - 3D brain volume
#'   - "H5NeuroVec" - 4D brain time series
#'   - "H5NeuroVecSeq" - Sequence of 4D scans
#'   - "H5ParcellatedMultiScan" - Multi-run parcellated experiment
#'   - "H5ParcellatedScan" - Single parcellated run (full data)
#'   - "H5ParcellatedScanSummary" - Single parcellated run (summary)
#'   - "LatentNeuroVec" - Latent representation
#'   - "LabeledVolumeSet" - Labeled brain regions
#'   - "unknown" - Type could not be determined
#'
#' @export
#' @examples
#' \dontrun{
#' # Detect type from file path
#' type <- detect_h5_type("data.h5")
#'
#' # Detect type from open H5File
#' h5 <- H5File$new("data.h5", "r")
#' type <- detect_h5_type(h5)
#' h5$close_all()
#' }
#' @export
detect_h5_type <- function(h5obj) {
  # Handle both file paths and H5File objects
  if (is.character(h5obj)) {
    if (!file.exists(h5obj)) {
      stop("File does not exist: ", h5obj)
    }
    h5 <- H5File$new(h5obj, "r")
    on.exit(h5$close_all())
  } else if (!is(h5obj, "H5File")) {
    stop("h5obj must be a file path or H5File object")
  } else {
    h5 <- h5obj
  }

  # First check for explicit class attribute (newest format)
  if (h5_attr_exists(h5, "fmristore_class")) {
    return(h5attr(h5, "fmristore_class"))
  }

  # Check for rtype attribute (older format)
  if (h5_attr_exists(h5, "rtype")) {
    rtype <- h5attr(h5, "rtype")

    # Map rtype to our class names
    type_map <- list(
      "DenseNeuroVol" = "H5NeuroVol",
      "DenseNeuroVec" = "H5NeuroVec",
      "SparseNeuroVol" = "H5NeuroVol",
      "SparseNeuroVec" = "H5NeuroVec"
    )

    if (rtype %in% names(type_map)) {
      return(type_map[[rtype]])
    }
  }

  # Heuristic detection based on HDF5 structure

  # Check for H5ParcellatedMultiScan structure
  if (exists_in_h5(h5, "/scans") && exists_in_h5(h5, "/cluster_map")) {
    return("H5ParcellatedMultiScan")
  }

  # Check for H5NeuroVecSeq structure (has scans but no cluster_map)
  if (exists_in_h5(h5, "/scans") &&
    exists_in_h5(h5, "/header") &&
    !exists_in_h5(h5, "/cluster_map")) {
    return("H5NeuroVecSeq")
  }

  # Check for LatentNeuroVec structure
  if (exists_in_h5(h5, "/basis/basis_matrix") ||
    exists_in_h5(h5, "/basis/basis_data")) {
    return("LatentNeuroVec")
  }

  # Check for LabeledVolumeSet structure
  if (exists_in_h5(h5, "/labels")) {
    # Check if it has the expected structure
    labels_exist <- tryCatch(
      {
        labels <- h5_read(h5, "/labels")
        length(labels) > 0
      },
      error = function(e) FALSE
    )

    if (labels_exist) {
      return("LabeledVolumeSet")
    }
  }

  # Check for basic NeuroVol/NeuroVec by data dimensions
  if (exists_in_h5(h5, "/data")) {
    # Get data dimensions
    data_dims <- tryCatch(
      {
        data_dset <- h5[["data"]]
        dims <- data_dset$dims
        data_dset$close()
        dims
      },
      error = function(e) NULL
    )

    if (!is.null(data_dims)) {
      if (length(data_dims) == 3) {
        return("H5NeuroVol")
      } else if (length(data_dims) == 4) {
        return("H5NeuroVec")
      }
    }
  }

  # Check space group for dimension clues
  if (exists_in_h5(h5, "/space/dim")) {
    dims <- tryCatch(
      {
        h5_read(h5, "/space/dim")
      },
      error = function(e) NULL
    )

    if (!is.null(dims)) {
      if (length(dims) == 3) {
        return("H5NeuroVol")
      } else if (length(dims) == 4) {
        return("H5NeuroVec")
      }
    }
  }

  # Check header for clues
  if (exists_in_h5(h5, "/header/dim")) {
    dims <- tryCatch(
      {
        h5_read(h5, "/header/dim")
      },
      error = function(e) NULL
    )

    if (!is.null(dims) && length(dims) >= 3) {
      # Check the first element which indicates number of dimensions
      if (dims[1] == 3 || (length(dims) >= 4 && dims[5] == 1)) {
        return("H5NeuroVol")
      } else if (dims[1] == 4 || (length(dims) >= 5 && dims[5] > 1)) {
        return("H5NeuroVec")
      }
    }
  }

  return("unknown")
}

# Read Dataset Methods ----

#' Read Dataset Method Implementation
#'
#' Internal function to read based on known type
#'
#' @param file Path to HDF5 file
#' @param type Character string indicating type
#' @param ... Additional arguments passed to constructors
#' @return The appropriate fmristore object
#' @keywords internal
read_dataset_typed <- function(file, type, ...) {
  switch(type,
    # H5-backed neuroimaging volumes
    "H5NeuroVol" = H5NeuroVol(file),

    # H5-backed neuroimaging vectors
    "H5NeuroVec" = H5NeuroVec(file, ...),

    # H5-backed sequence of vectors
    "H5NeuroVecSeq" = H5NeuroVecSeq(file),

    # Parcellated multi-scan experiments
    "H5ParcellatedMultiScan" = H5ParcellatedMultiScan(file, ...),

    # Individual parcellated scans (usually not read directly)
    "H5ParcellatedScan" = {
      warning("H5ParcellatedScan is typically accessed through H5ParcellatedMultiScan")
      H5ParcellatedScan(file, ...)
    },
    "H5ParcellatedScanSummary" = {
      warning("H5ParcellatedScanSummary is typically accessed through H5ParcellatedMultiScan")
      H5ParcellatedScanSummary(file, ...)
    },

    # Latent representations
    "LatentNeuroVec" = {
      # LatentNeuroVec is not H5-backed, need to read and construct
      h5 <- H5File$new(file, "r")
      on.exit(h5$close_all())

      # Read mask
      mask <- if (exists_in_h5(h5, "/mask")) {
        mask_data <- h5_read(h5, "/mask")
        # Need space information to create LogicalNeuroVol
        space_dims <- h5_read(h5, "/header/dim")[2:4]
        mask_array <- array(as.logical(mask_data), dim = space_dims)
        # Create LogicalNeuroVol (assuming neuroim2 is available)
        neuroim2::LogicalNeuroVol(
          mask_array,
          neuroim2::NeuroSpace(dim = space_dims)
        )
      } else {
        stop("LatentNeuroVec file missing mask")
      }

      # Read basis
      basis <- if (exists_in_h5(h5, "/basis/basis_matrix")) {
        t(h5_read(h5, "/basis/basis_matrix")) # Transpose to get time x components
      } else if (exists_in_h5(h5, "/basis/basis_data")) {
        t(h5_read(h5, "/basis/basis_data"))
      } else {
        stop("LatentNeuroVec file missing basis data")
      }

      # Read loadings from first scan
      loadings <- if (exists_in_h5(h5, "/scans")) {
        scan_names <- h5[["scans"]]$ls()$name
        if (length(scan_names) > 0) {
          first_scan <- scan_names[1]
          embedding_path <- paste0("/scans/", first_scan, "/embedding")
          if (exists_in_h5(h5, embedding_path)) {
            t(h5_read(h5, embedding_path)) # Transpose to get voxels x components
          } else {
            stop("LatentNeuroVec file missing embedding data")
          }
        } else {
          stop("LatentNeuroVec file has no scans")
        }
      } else {
        stop("LatentNeuroVec file missing scans")
      }

      # Read offset if present
      offset <- if (exists_in_h5(h5, "/offset")) {
        h5_read(h5, "/offset")
      } else {
        NULL
      }

      # Create space
      space_dims <- h5_read(h5, "/header/dim")[2:5]
      space <- neuroim2::NeuroSpace(dim = space_dims)

      # Construct LatentNeuroVec
      LatentNeuroVec(
        basis = basis,
        loadings = loadings,
        space = space,
        mask = mask,
        offset = offset
      )
    },

    # Labeled volume sets
    "LabeledVolumeSet" = read_labeled_vec(file),

    # Unknown type
    stop("Unsupported type: ", type)
  )
}

#' Read Dataset from HDF5 File
#'
#' @description
#' Reads a dataset from an HDF5 file with automatic type detection.
#' The function examines the file structure to determine the appropriate
#' object type and returns the corresponding fmristore object.
#'
#' @param file Path to HDF5 file
#' @param type Optional character string specifying the expected type.
#'   If NULL (default), the type is auto-detected.
#' @param ... Additional arguments passed to the object constructor
#'
#' @return An fmristore object of the appropriate type:
#'   - \code{H5NeuroVol} for 3D volumes
#'   - \code{H5NeuroVec} for 4D time series
#'   - \code{H5NeuroVecSeq} for sequences of 4D data
#'   - \code{H5ParcellatedMultiScan} for parcellated experiments
#'   - \code{LatentNeuroVec} for latent representations
#'   - \code{LabeledVolumeSet} for labeled regions
#'
#' @export
#' @examples
#' \dontrun{
#' # Auto-detect type
#' obj <- read_dataset("data.h5")
#'
#' # Specify type explicitly
#' obj <- read_dataset("data.h5", type = "H5NeuroVec")
#'
#' # Pass additional arguments to constructor
#' obj <- read_dataset("data.h5", dataset_name = "processed")
#' }
setMethod(
  "read_dataset",
  signature(x = "character"),
  function(x, type = NULL, ...) {
    if (!file.exists(x)) {
      stop("File does not exist: ", x)
    }

    # Use explicit type if provided
    if (!is.null(type)) {
      message("Reading ", type, " from: ", x)
      return(read_dataset_typed(x, type, ...))
    }

    # Auto-detect type
    detected_type <- detect_h5_type(x)

    if (detected_type == "unknown") {
      stop(
        "Cannot determine dataset type for file: ", x,
        "\nTry specifying 'type' explicitly or check file structure",
        "\nYou can use detect_h5_type('", x, "') to investigate"
      )
    }

    message("Detected type '", detected_type, "' in file: ", x)
    read_dataset_typed(x, detected_type, ...)
  }
)

#' Read Dataset from H5File Object
#'
#' @rdname read_dataset
#' @export
setMethod(
  "read_dataset",
  signature(x = "H5File"),
  function(x, type = NULL, ...) {
    # For H5File objects, we need the file path for constructors
    file_path <- x$filename

    # Use explicit type if provided
    if (!is.null(type)) {
      message("Reading ", type, " from H5File")
      return(read_dataset_typed(file_path, type, ...))
    }

    # Auto-detect type using the H5File object
    detected_type <- detect_h5_type(x)

    if (detected_type == "unknown") {
      stop(
        "Cannot determine dataset type from H5File object",
        "\nTry specifying 'type' explicitly"
      )
    }

    message("Detected type '", detected_type, "' in H5File")
    read_dataset_typed(file_path, detected_type, ...)
  }
)
