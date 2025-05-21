Preface:

Okay, I've reviewed the codebase with an eye for simplifying HDF5 interactions using abstractions. Here are some key opportunities:

**1. Consistent Use of Existing `h5_utils.R` Helpers:**

*   **`h5_read()` for Metadata**:
    *   In `R/h5neurovec.R` (constructor `H5NeuroVec`) and `R/h5neurovol.R` (constructor `H5NeuroVol`), direct access like `h5obj[["space/dim"]][]` is used to read spatial metadata.
    *   **Suggestion**: Replace with `fmristore:::h5_read(h5obj, "space/dim")`. This centralizes error handling and makes the code more uniform.
        ```R
        # Example in H5NeuroVec constructor (R/h5neurovec.R)
        # Current:
        # sp <- NeuroSpace(
        #   dim    = h5obj[["space/dim"]][],
        #   origin = h5obj[["space/origin"]][],
        #   trans  = h5obj[["space/trans"]][,]
        # )
        # Suggested:
        sp_dim_val <- fmristore:::h5_read(h5obj, "space/dim", missing_ok = FALSE)
        sp_origin_val <- fmristore:::h5_read(h5obj, "space/origin", missing_ok = FALSE)
        sp_trans_val <- fmristore:::h5_read(h5obj, "space/trans", missing_ok = FALSE)
        sp <- NeuroSpace(
          dim    = sp_dim_val,
          origin = sp_origin_val,
          trans  = sp_trans_val
        )
        ```
    *   Similarly in `R/cluster_experiment.R` (constructor `H5ClusterExperiment`) for reading `/header/*` fields.

*   **`with_h5_dataset()` for Scoped Reads**:
    *   In `R/h5neurovol.R` (methods `linear_access` and `[`), there's manual `tryCatch/finally` for opening, reading, and closing datasets.
    *   **Suggestion**: Use `fmristore:::with_h5_dataset()` to simplify handle management.
        ```R
        # Example in H5NeuroVol's linear_access (R/h5neurovol.R)
        # Current complex tryCatch/finally for dset:
        # dset <- NULL; subvol <- NULL
        # tryCatch({
        #     dset <- x@h5obj[["data/elements"]]
        #     subvol <- dset[minx:maxx, miny:maxy, minz:maxz, drop=FALSE]
        # }, finally = { close_h5_safely(dset) })
        # if (is.null(subvol)) stop(...)

        # Suggested:
        subvol <- tryCatch({
            fmristore:::with_h5_dataset(x@h5obj, "data/elements", function(dset_handle) {
                dset_handle[minx:maxx, miny:maxy, minz:maxz, drop=FALSE]
            })
        }, error = function(e) {
            stop("Failed to read sub-volume data for linear_access. Original error: ", e$message)
        })
        ```
    *   This pattern also appears in `R/io_h5_helpers.R` (`read_h5_mask_to_LogicalNeuroVol`, `read_h5_clusters_to_ClusteredNeuroVol`) and could be refactored.
    *   In `R/cluster_array.R` (`as.matrix` for `H5ClusterRunSummary` and the deprecated `make_run_summary`), direct dataset reads with `tryCatch/finally` can also use `with_h5_dataset`.

*   **`h5_read_subset()` where applicable**:
    *   In `R/cluster_array.R` (`.get_cluster_timeseries_by_mask_index`), the line `cluster_data_subset <- ds[row_offsets_in_dataset, time_indices, drop = FALSE]` could potentially use `fmristore:::h5_read_subset(x@obj, dset_path, index = list(row_offsets_in_dataset, time_indices))` if the dataset handle `ds` wasn't already open and managed. However, `ds` is opened just before. If `h5_read_subset` were more flexible or if the logic was slightly restructured, it might fit. For now, the current use of `assert_h5_path` and direct subsetting with `close_h5_safely` is reasonable.

**2. New Abstraction: Reading HDF5 Attributes (`h5_read_attr`)**

*   **Problem**: Reading attributes often involves `try(hdf5r::h5attr(...), silent=TRUE)` or checking `h5attr_exists` first. This is seen in `H5NeuroVol()`, `H5NeuroVec()`, `LatentNeuroVecSource`'s `load_data`, and `H5ClusterExperiment()`.
*   **Suggestion**: Create a helper in `R/h5_utils.R`:
    ```R
    #' Read an HDF5 attribute safely
    #' @keywords internal
    h5_read_attr <- function(h5_obj_or_grp, attr_name, default = NULL, missing_ok = FALSE) {
      if (!inherits(h5_obj_or_grp, "H5RefClass")) { # H5File or H5Group inherit from H5RefClass
        stop("'h5_obj_or_grp' must be an H5File or H5Group object.")
      }
      if (!h5_obj_or_grp$attr_exists(attr_name)) {
        if (missing_ok) return(default)
        stop(sprintf("Attribute '%s' not found in HDF5 object/group '%s'", 
                     attr_name, h5_obj_or_grp$get_obj_name()))
      }
      tryCatch({
        hdf5r::h5attr(h5_obj_or_grp, attr_name)
      }, error = function(e) {
        if (missing_ok) return(default)
        stop(sprintf("Error reading attribute '%s': %s", attr_name, conditionMessage(e)))
      })
    }
    ```
    *   Usage: `rtype <- fmristore:::h5_read_attr(h5obj, "rtype", missing_ok = TRUE)`

**3. New Abstraction: Writing a List of Datasets/Attributes to a Group (`h5_write_group_from_list`)**

*   **Problem**: Writing structured metadata like NIfTI-style headers (`write_labeled_vec`, `.write_header` in `latent_vec.R`) or other grouped datasets (`H5ClusterExperiment` metadata) involves many repetitive calls to `h5_write()` or direct `hdf5r` calls.
*   **Suggestion**: A helper in `R/h5_utils.R`:
    ```R
    #' Write elements of a list as datasets within an HDF5 group
    #' @keywords internal
    h5_write_group_from_list <- function(h5_parent_obj, group_path, data_list, 
                                         create_missing_groups = TRUE, overwrite_group = FALSE) {
      if (overwrite_group && h5_parent_obj$exists(group_path)) {
        h5_parent_obj$link_delete(group_path)
      }
      if (create_missing_groups) {
        fmristore:::ensure_h5_groups(h5_parent_obj, group_path) # Ensures group_path itself exists
      } else if (!h5_parent_obj$exists(group_path)) {
        stop(sprintf("Target group '%s' does not exist and create_missing_groups is FALSE.", group_path))
      }
      
      for (name in names(data_list)) {
        dataset_path_in_group <- file.path(group_path, name) # Correctly forms full path
        # h5_write handles overwrite logic for individual datasets within the group
        fmristore:::h5_write(h5_parent_obj, dataset_path_in_group, data_list[[name]], 
                             overwrite = TRUE) # Assuming overwrite within the group is desired
      }
      invisible(NULL)
    }
    ```
    *   Usage: `.write_header` in `latent_vec.R` becomes:
        `fmristore:::h5_write_group_from_list(h5, "/header", hdr_fields)` (after ensuring `qfac` is in `hdr_fields` or handled separately if it's an attribute vs dataset).
    *   A similar helper, `h5_write_attributes_from_list(h5_obj_or_grp, attr_list)`, could be made for attributes.

**4. Potential Abstraction: Sparse Matrix I/O**

*   **Problem**: In `R/latent_vec.R` (`.write_basis` and `load_data`), writing and reading sparse matrices (CSC format) involves multiple specific HDF5 operations for `data`, `indices`, `indptr` and attributes like `storage` and `shape`.
*   **Suggestion**: If this pattern is (or becomes) common, create:
    *   `h5_write_sparse_matrix(h5_obj_or_grp, path_base_name, matrix_obj, format="csc", compression=0, ...)`
    *   `h5_read_sparse_matrix(h5_obj_or_grp, path_base_name)`
    *   These would create/read a group named `path_base_name` containing the datasets and attributes.
    *   This would significantly clean up `.write_basis` and the sparse reading part of `load_data` in `latent_vec.R`.

**5. Streamlining `n_time` Determination in `H5ClusterRun`**

*   **Problem**: The logic in `H5ClusterRun` constructor (and deprecated `make_run_full`) to determine `n_time` from attributes, then metadata datasets, then by inferring from a cluster dataset, is verbose and involves direct HDF5 calls.
*   **Suggestion**: Refactor this logic to use `h5_read_attr` (from suggestion #2) and `h5_read(missing_ok=TRUE)` to simplify the HDF5 access parts. A dedicated helper `h5_resolve_n_time(h5obj, scan_group_path, first_cluster_ref_path)` might be an option, but careful use of existing/new small helpers within the constructor might be sufficient.

**Example Refactor for `n_time` determination in `H5ClusterRun` (conceptual):**
```R
# Inside H5ClusterRun constructor (R/constructors.R)
if (is.null(determined_n_time)) {
  scan_group_path <- paste0("/scans/", scan_name)
  if (h5obj$exists(scan_group_path)) {
    # Try attribute first
    determined_n_time <- fmristore:::h5_read_attr(h5obj[[scan_group_path]], "n_time", missing_ok = TRUE) # Assumes h5_read_attr can take group handle
    
    if (is.null(determined_n_time)) {
      # Try metadata dataset
      determined_n_time <- fmristore:::h5_read(h5obj, file.path(scan_group_path, "metadata/n_time"), missing_ok = TRUE)
    }
    
    if (is.null(determined_n_time) && h5obj$exists(file.path(scan_group_path, "clusters"))) {
      # Try inferring (simplified logic here, full logic is more complex)
      first_cid <- clusters@clusters[1] # Assuming clusters object is available
      if (!is.null(first_cid) && !is.na(first_cid)) {
          dset_path_cid1 <- sprintf("%s/clusters/cluster_%d", scan_group_path, as.integer(first_cid))
          if (h5obj$exists(dset_path_cid1)) {
              dims_cid1 <- fmristore:::h5_dataset_dims(h5obj, dset_path_cid1) # New helper h5_dataset_dims
              if (length(dims_cid1) == 2) determined_n_time <- dims_cid1[2]
          }
      }
    }
  }
  # ... rest of validation ...
}
```
This uses the proposed `h5_read_attr` and an existing `h5_read`. A new `h5_dataset_dims` (which would use `with_h5_dataset` internally) is also shown.

By consistently applying these existing and new small, focused helpers, the HDF5 interaction code across the package can become more robust, readable, and maintainable.The codebase already has a good set of HDF5 utility functions in `R/h5_utils.R`. The primary opportunity is to ensure these are used consistently and to identify a few more targeted abstractions for common, verbose patterns.

Here are specific suggestions:

**1. Consistent Use of Existing `h5_utils.R` Helpers:**

*   **`h5_read()` for Metadata and Single Datasets:**
    *   **Location(s):**
        *   `R/h5neurovec.R`: `H5NeuroVec()` constructor (reading `space/dim`, `space/origin`, `space/trans`).
        *   `R/h5neurovol.R`: `H5NeuroVol()` constructor (reading `space/dim`, `space/origin`, `space/trans`).
        *   `R/cluster_experiment.R`: `H5ClusterExperiment()` constructor (reading `/header/*` fields).
        *   `R/latent_vec.R`: `load_data()` method for `LatentNeuroVecSource` (reading various components like `/mask`, `/voxel_coords`, `/basis/basis_matrix`, `/offset`, scan embedding).
    *   **Current Pattern:** Direct access like `h5obj[["path/to/dataset"]][]`.
    *   **Refactor:** Replace with `fmristore:::h5_read(h5obj, "path/to/dataset", missing_ok = FALSE)`. This standardizes error handling and dataset access.

*   **`with_h5_dataset()` for Scoped Dataset Reads/Operations:**
    *   **Location(s):**
        *   `R/h5neurovol.R`: `linear_access()` and `[` methods (manual `tryCatch/finally` for dataset handle).
        *   `R/io_h5_helpers.R`: `read_h5_mask_to_LogicalNeuroVol()` and `read_h5_clusters_to_ClusteredNeuroVol()` (manual `tryCatch/finally`).
        *   `R/cluster_array.R`: `as.matrix.H5ClusterRunSummary()` and deprecated `make_run_summary()` (direct dataset read with manual `tryCatch/finally`).
    *   **Current Pattern:**
        ```R
        dset <- NULL; # ...
        tryCatch({
            dset <- x@h5obj[["data/elements"]]
            # ... operations with dset ...
        }, finally = {
            if (!is.null(dset) && inherits(dset, "H5D") && dset$is_valid) {
                close_h5_safely(dset)
            }
        })
        ```
    *   **Refactor:**
        ```R
        result <- tryCatch({
            fmristore:::with_h5_dataset(h5_handle, "path/to/dataset", function(dset_handle) {
                # ... operations with dset_handle ...
                # return what's needed
            })
        }, error = function(e) { /* handle error */ })
        ```
        This simplifies resource management significantly.

**2. New Abstraction: Safe Attribute Reading (`h5_read_attr`)**

*   **Problem:** Reading HDF5 attributes often involves `try(hdf5r::h5attr(...), silent=TRUE)` or manual checks with `h5obj$attr_exists()` followed by `h5attr()`. This is verbose and error-prone.
*   **Location(s):**
    *   `R/h5neurovec.R`: `H5NeuroVec()` (reading `rtype`).
    *   `R/h5neurovol.R`: `H5NeuroVol()` (reading `rtype`).
    *   `R/latent_vec.R`: `load_data()` for `LatentNeuroVecSource` (reading `latent_spec_version`).
    *   `R/cluster_experiment.R`: `H5ClusterExperiment()` constructor (reading `/scans@summary_only`).
*   **Suggestion:** Add `h5_read_attr` to `R/h5_utils.R`:
    ```R
    #' Read an HDF5 attribute safely
    #'
    #' Reads an attribute from an HDF5 object (file or group), with options
    #' for handling missing attributes and providing default values.
    #'
    #' @param h5_obj An open H5File or H5Group object.
    #' @param attr_name The name of the attribute to read.
    #' @param default The value to return if the attribute is missing and `missing_ok = TRUE`.
    #' @param missing_ok If `TRUE`, returns `default` if the attribute is not found.
    #'                   If `FALSE` (default), stops with an error if missing.
    #' @return The attribute value, or `default` if applicable.
    #' @keywords internal
    h5_read_attr <- function(h5_obj, attr_name, default = NULL, missing_ok = FALSE) {
      if (!is(h5_obj, "H5RefClass")) { # Covers H5File, H5Group
        stop("'h5_obj' must be an H5File or H5Group object.")
      }
      if (!is.character(attr_name) || length(attr_name) != 1 || !nzchar(attr_name)) {
        stop("'attr_name' must be a single, non-empty character string.")
      }

      if (!h5_obj$attr_exists(attr_name)) {
        if (missing_ok) {
          return(default)
        } else {
          stop(sprintf("Attribute '%s' not found in HDF5 object '%s'",
                       attr_name, h5_obj$get_obj_name()))
        }
      }
      tryCatch({
        hdf5r::h5attr(h5_obj, attr_name)
      }, error = function(e) {
        # This case might be redundant if attr_exists works, but good for robustness
        if (missing_ok) return(default)
        stop(sprintf("Error reading attribute '%s' from HDF5 object '%s': %s",
                     attr_name, h5_obj$get_obj_name(), conditionMessage(e)))
      })
    }
    ```
    This centralizes existence checks and error handling for attribute reading.

**3. New Abstraction: Writing a List of Datasets to a Group (`h5_write_group_from_list`)**

*   **Problem:** Writing structured metadata, like NIfTI-style headers, involves many repetitive calls to `h5_write()` for each field.
*   **Location(s):**
    *   `R/labeled_vec.R`: `write_labeled_vec()` (writing header fields).
    *   `R/latent_vec.R`: Internal helper `.write_header()` (writing header fields for latent spec).
*   **Suggestion:** Add `h5_write_group_from_list` to `R/h5_utils.R`:
    ```R
    #' Write elements of a list as datasets/attributes within an HDF5 group
    #'
    #' Iterates over a named list and writes each element as a dataset (default)
    #' or attribute under the specified HDF5 group.
    #'
    #' @param h5_parent_obj An open H5File or H5Group object that will contain `group_path`.
    #' @param group_path The path (relative to `h5_parent_obj` if it's H5File, or direct name if H5Group)
    #'                   for the group to create/use.
    #' @param data_list A named list where names become dataset/attribute names and
    #'                  values are the data to write.
    #' @param as_attributes Logical, if `TRUE`, writes list elements as attributes of `group_path`
    #'                      instead of datasets within it. Default `FALSE`.
    #' @param create_missing_groups Logical, if `TRUE` (default), recursively creates `group_path`.
    #' @param overwrite_group Logical, if `TRUE`, deletes `group_path` if it exists before writing. Default `FALSE`.
    #'                        Note: Individual datasets/attributes within the group will be overwritten
    #'                        by `h5_write` or `h5attr<-` if they exist.
    #' @keywords internal
    h5_write_group_from_list <- function(h5_parent_obj, group_path, data_list,
                                         as_attributes = FALSE,
                                         create_missing_groups = TRUE,
                                         overwrite_group = FALSE) {
      if (!is(h5_parent_obj, "H5RefClass")) stop("'h5_parent_obj' must be H5File or H5Group.")
      if (!is.character(group_path) || !nzchar(group_path)) stop("'group_path' must be a non-empty string.")
      if (!is.list(data_list) || is.null(names(data_list))) stop("'data_list' must be a named list.")

      full_group_path <- if (group_path == "/" || startsWith(group_path, "/")) {
        group_path # Absolute path
      } else {
        file.path(h5_parent_obj$get_obj_name(), group_path) # Relative path
      }
      
      target_h5_obj_for_group_ops <- if (is(h5_parent_obj, "H5File")) h5_parent_obj else h5_parent_obj$get_file()

      if (overwrite_group && target_h5_obj_for_group_ops$exists(full_group_path)) {
        target_h5_obj_for_group_ops$link_delete(full_group_path)
      }

      if (create_missing_groups) {
        fmristore:::ensure_h5_groups(target_h5_obj_for_group_ops, full_group_path)
      } else if (!target_h5_obj_for_group_ops$exists(full_group_path)) {
        stop(sprintf("Target group '%s' does not exist and create_missing_groups is FALSE.", full_group_path))
      }

      group_handle_for_elements <- target_h5_obj_for_group_ops[[full_group_path]]
      on.exit(fmristore:::close_h5_safely(group_handle_for_elements), add = TRUE)

      for (name in names(data_list)) {
        value <- data_list[[name]]
        if (as_attributes) {
          hdf5r::h5attr(group_handle_for_elements, name) <- value
        } else {
          # Use h5_write, ensuring full path by operating on target_h5_obj_for_group_ops
          element_full_path <- file.path(full_group_path, name)
          fmristore:::h5_write(target_h5_obj_for_group_ops, element_full_path, value,
                               overwrite = TRUE) # Overwrite dataset if it exists
        }
      }
      invisible(NULL)
    }
    ```
    This would make header writing more concise.

**4. New Abstraction: Sparse Matrix I/O (If Widely Used)**

*   **Problem:** `R/latent_vec.R` contains bespoke logic for writing (`.write_basis`) and reading (`load_data`) sparse matrices in CSC format (attributes for `storage`, `shape`; datasets for `data`, `indices`, `indptr`).
*   **Location(s):** `R/latent_vec.R`.
*   **Suggestion:** If this sparse matrix storage format is intended to be a package standard or used in more than one place:
    *   `h5_write_sparse_matrix(h5_parent_obj, group_name_for_sparse, matrix_obj, format="csc", compression=0, ...)`
    *   `h5_read_sparse_matrix(h5_parent_obj, group_name_for_sparse)`
    These helpers would encapsulate the specific HDF5 structure for sparse matrices. Given it's currently localized to `latent_vec.R`, this might be a lower priority unless reuse is anticipated. The current implementation within `latent_vec.R` leverages `h5_write` for its components, which is good.

**5. Refactoring `n_time` Determination in H5ClusterRun and H5ClusterRunSummary Constructors**

*   **Problem:** The logic to determine `n_time` (from HDF5 attribute, then metadata dataset, then inferring from a cluster dataset) in `R/constructors.R` (for `H5ClusterRun`) and the deprecated `make_run_full` is somewhat verbose with direct `hdf5r` calls.
*   **Suggestion:** Refactor this logic internally using `h5_read_attr` (from suggestion #2) and `h5_read(missing_ok=TRUE)`. Add `fmristore:::h5_dataset_dims(h5_obj, path)` as a new helper in `R/h5_utils.R` (which internally uses `with_h5_dataset`) to cleanly get dataset dimensions for inference.
    ```R
    # In R/h5_utils.R
    #' Obtain the dimensions of an HDF5 dataset
    #' @keywords internal
    h5_dataset_dims <- function(h5, path) {
      fmristore:::with_h5_dataset(h5, path, function(ds) ds$dims)
    }
    ```
    This makes the `n_time` resolution logic in the constructors cleaner and more reliant on the established utility functions.

By implementing suggestions 1, 2, 3, and 5, the codebase would see significant improvements in clarity, robustness, and maintainability of its HDF5 interactions. Suggestion 4 (sparse matrix I/O) is valuable if that pattern becomes more common.
---

**Phase 1: Implement New Core Utility Functions**

**Ticket 1: Implement `h5_read_attr` Utility Function**

*   **Description:** Reading HDF5 attributes currently involves manual checks (`attr_exists`) and `tryCatch` blocks in multiple places. A centralized helper is needed.
*   **Task:**
    1.  Create a new function `h5_read_attr(h5_obj, attr_name, default = NULL, missing_ok = FALSE)` in `R/h5_utils.R`.
    2.  The function should:
        *   Validate `h5_obj` (must be `H5File` or `H5Group`).
        *   Validate `attr_name` (character string).
        *   If attribute exists, read and return it.
        *   If attribute does not exist:
            *   If `missing_ok = TRUE`, return `default`.
            *   If `missing_ok = FALSE`, stop with an informative error.
        *   Handle potential errors during attribute reading gracefully.
    3.  Add roxygen documentation and mark as `@keywords internal`.
    4.  Add unit tests for various scenarios (attribute exists, missing with `ok=TRUE/FALSE`, error during read).
*   **Files to Modify:** `R/h5_utils.R`, `tests/testthat/test-h5_utils.R` (new or existing test file for utils).
*   **Acceptance Criteria:** Function implemented, documented, and unit tested.

**Ticket 2: Implement `h5_write_group_from_list` Utility Function**

*   **Description:** Writing multiple datasets or attributes to an HDF5 group from a list is a common pattern (e.g., NIfTI headers). This is currently done with repetitive `h5_write` or direct `hdf5r` calls.
*   **Task:**
    1.  Create a new function `h5_write_group_from_list(h5_parent_obj, group_path, data_list, as_attributes = FALSE, create_missing_groups = TRUE, overwrite_group = FALSE)` in `R/h5_utils.R`.
    2.  The function should:
        *   Validate inputs.
        *   Handle `group_path` creation/deletion based on `create_missing_groups` and `overwrite_group`.
        *   Iterate through `data_list`:
            *   If `as_attributes = TRUE`, write each element as an attribute of `group_path` using `hdf5r::h5attr()<-`.
            *   If `as_attributes = FALSE`, write each element as a dataset under `group_path` using `fmristore:::h5_write()`, ensuring individual datasets are overwritten if they exist.
    3.  Add roxygen documentation and mark as `@keywords internal`.
    4.  Add unit tests for writing datasets and attributes, handling group creation/overwrite.
*   **Files to Modify:** `R/h5_utils.R`, `tests/testthat/test-h5_utils.R`.
*   **Acceptance Criteria:** Function implemented, documented, and unit tested.

**Ticket 3: Implement `h5_dataset_dims` Utility Function**

*   **Description:** A small helper to get dataset dimensions cleanly, often needed before reading or for validation.
*   **Task:**
    1.  Create a new function `h5_dataset_dims(h5_obj, path)` in `R/h5_utils.R`.
    2.  The function should use `fmristore:::with_h5_dataset(h5_obj, path, function(ds) ds$dims)` to retrieve dimensions.
    3.  Add roxygen documentation and mark as `@keywords internal`.
    4.  Add basic unit tests.
*   **Files to Modify:** `R/h5_utils.R`, `tests/testthat/test-h5_utils.R`.
*   **Acceptance Criteria:** Function implemented, documented, and unit tested.

---

**Phase 2: Refactor Existing Code to Use Utility Functions**

**Ticket 4: Refactor Metadata Reading in Constructors to Use `fmristore:::h5_read`**

*   **Description:** Constructors for `H5NeuroVec`, `H5NeuroVol`, `H5ClusterExperiment`, and the `load_data` method for `LatentNeuroVecSource` currently use direct HDF5 access (`h5obj[["path"]][]`) for reading metadata datasets.
*   **Task:**
    1.  Modify `R/h5neurovec.R` (`H5NeuroVec()` constructor): Replace direct reads of `space/dim`, `space/origin`, `space/trans` with `fmristore:::h5_read(..., missing_ok = FALSE)`.
    2.  Modify `R/h5neurovol.R` (`H5NeuroVol()` constructor): Replace direct reads of `space/dim`, `space/origin`, `space/trans` with `fmristore:::h5_read(..., missing_ok = FALSE)`.
    3.  Modify `R/cluster_experiment.R` (`H5ClusterExperiment()` constructor): Replace direct reads of `/header/*` fields with `fmristore:::h5_read(..., missing_ok = TRUE/FALSE)` as appropriate.
    4.  Modify `R/latent_vec.R` (`load_data` for `LatentNeuroVecSource`): Replace direct reads of `/mask`, `/voxel_coords`, `/basis/basis_matrix`, `/offset`, and scan embedding with `fmristore:::h5_read(..., missing_ok = TRUE/FALSE)`.
*   **Acceptance Criteria:**
    *   Code updated in specified files.
    *   Existing functionality and tests for these constructors/methods remain passing.
    *   Error handling for missing required datasets is now managed by `h5_read`.

**Ticket 5: Refactor Direct Dataset Operations to Use `fmristore:::with_h5_dataset`**

*   **Description:** Several places manually open, operate on, and close HDF5 dataset handles using `tryCatch/finally`.
*   **Task:**
    1.  Modify `R/h5neurovol.R` (`linear_access()` and `[` methods): Replace manual dataset handle management with `fmristore:::with_h5_dataset()`.
    2.  Modify `R/io_h5_helpers.R` (`read_h5_mask_to_LogicalNeuroVol()` and `read_h5_clusters_to_ClusteredNeuroVol()`): Refactor to use `fmristore:::with_h5_dataset()` for reading.
    3.  Modify `R/cluster_array.R` (`as.matrix.H5ClusterRunSummary()` and deprecated `make_run_summary`): Refactor dataset reading to use `fmristore:::with_h5_dataset()`.
*   **Acceptance Criteria:**
    *   Code updated in specified files.
    *   Resource management for dataset handles is simplified.
    *   Existing functionality and tests for these methods remain passing.

**Ticket 6: Refactor Attribute Reading to Use `fmristore:::h5_read_attr`**

*   **Description:** Several places read HDF5 attributes using direct `hdf5r` calls and manual error/existence checks.
*   **Task:** (Depends on Ticket 1)
    1.  Modify `R/h5neurovec.R` (`H5NeuroVec()` constructor): Use `fmristore:::h5_read_attr()` to read the `rtype` attribute.
    2.  Modify `R/h5neurovol.R` (`H5NeuroVol()` constructor): Use `fmristore:::h5_read_attr()` to read the `rtype` attribute.
    3.  Modify `R/latent_vec.R` (`load_data` for `LatentNeuroVecSource`): Use `fmristore:::h5_read_attr()` to read `latent_spec_version` and sparse matrix attributes (`storage`, `shape`).
    4.  Modify `R/cluster_experiment.R` (`H5ClusterExperiment()` constructor): Use `fmristore:::h5_read_attr()` to read the `/scans@summary_only` attribute.
    5.  Modify `R/constructors.R` (`H5ClusterRun()` constructor): If applicable, use `fmristore:::h5_read_attr()` for reading `n_time` from `/scans/<scan_name>@n_time`.
*   **Acceptance Criteria:**
    *   Code updated in specified files.
    *   Attribute reading is now centralized and robust.
    *   Existing functionality and tests for these methods remain passing.

**Ticket 7: Refactor Grouped Dataset Writing to Use `fmristore:::h5_write_group_from_list`**

*   **Description:** Writing NIfTI-style headers or other grouped metadata involves multiple `h5_write` calls.
*   **Task:** (Depends on Ticket 2)
    1.  Modify `R/labeled_vec.R` (`write_labeled_vec()`): Refactor the writing of `/header` fields to use `fmristore:::h5_write_group_from_list(..., as_attributes = FALSE)`. Handle `qfac` separately if it needs to be an attribute vs. dataset consistently.
    2.  Modify `R/latent_vec.R` (internal helper `.write_header` called by `to_h5_latentvec`): Refactor writing of `/header` fields to use `fmristore:::h5_write_group_from_list(..., as_attributes = FALSE)`.
*   **Acceptance Criteria:**
    *   Code updated in specified files.
    *   Header/metadata writing is more concise.
    *   Output HDF5 structure remains compatible and existing tests pass.

**Ticket 8: Refactor `n_time` Determination in `H5ClusterRun` and `H5ClusterRunSummary` Constructors**

*   **Description:** The logic for determining `n_time` (from attribute, metadata dataset, or inference) in `R/constructors.R` can be streamlined.
*   **Task:** (Depends on Ticket 1 and Ticket 3)
    1.  In `R/constructors.R` (`H5ClusterRun()`):
        *   Use `fmristore:::h5_read_attr(..., missing_ok=TRUE)` to attempt reading `n_time` from the scan group attribute.
        *   If not found, use `fmristore:::h5_read(..., missing_ok=TRUE)` to attempt reading from the `metadata/n_time` dataset.
        *   If still not found, use `fmristore:::h5_dataset_dims()` to get dimensions from the first cluster dataset for inference.
        *   Ensure the final `determined_n_time` is validated.
    2.  In `R/constructors.R` (`H5ClusterRunSummary()`):
        *   Use `fmristore:::h5_dataset_dims()` to get `n_time` from the summary dataset dimensions.
*   **Acceptance Criteria:**
    *   `n_time` determination logic in specified constructors is refactored and cleaner.
    *   The order of precedence for finding `n_time` is maintained.
    *   Existing functionality and tests for these constructors pass.

---

**Phase 3: Potential Future Enhancements (Lower Priority / If Pattern Repeats)**

**Ticket 9: Investigate Abstraction for Sparse Matrix I/O**

*   **Description:** Writing and reading sparse matrices in CSC/CSR format (as done in `R/latent_vec.R`) involves several specific HDF5 operations. If this pattern is reused or becomes complex, an abstraction could be beneficial.
*   **Task:**
    1.  Analyze the sparse matrix I/O logic in `R/latent_vec.R` (`.write_basis`, `load_data`).
    2.  Determine if the complexity and potential for reuse warrant dedicated helper functions like `h5_write_sparse_matrix()` and `h5_read_sparse_matrix()` in `R/h5_utils.R`.
    3.  If deemed beneficial, implement these helpers and refactor `R/latent_vec.R` to use them.
*   **Acceptance Criteria:**
    *   Analysis completed.
    *   If implemented, new helpers are documented and tested.
    *   `R/latent_vec.R` is refactored if helpers are created, and existing tests pass.
*   **Note:** This is more exploratory. The current implementation within `latent_vec.R` is localized and might be acceptable if not reused.

---

These tickets should provide a clear path to refactoring the HDF5 interactions. Each ticket aims to be a manageable chunk of work.