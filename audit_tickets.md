please see high-level audit.md for the code review.

## fmristore Code Review - Actionable Tickets

** When you complete a ticket, please update this audit_tickets.md file to reflect the changes. **

### 🟥 High Priority / Potential Bugs

**[x] Ticket 1: Resolve Multiple HDF5 Handles in `H5ClusteredExperiment` Constructor**
    *   **Issue**: `H5ClusteredExperiment` constructor, when given a file path, opens the HDF5 file. It then passes the *path* (not the handle) to `H5ClusteredRunFull`/`Summary` constructors, causing them to re-open the same file, leading to multiple handles.
    *   **File(s)**: `R/cluster_experiment.R` (main constructor), `R/constructors.R` (run constructors, to ensure they correctly use `open_h5` when passed an existing handle).
    *   **Recommendation**:
        1.  Modify `H5ClusteredExperiment` constructor: if it opens the file (because a path was provided), it should pass the *opened `H5File` object* (e.g., `fh$h5`) to the `H5ClusteredRunFull`/`Summary` constructors.
        2.  Ensure `H5ClusteredRunFull` and `H5ClusteredRunSummary` constructors correctly use `open_h5` internally so they can accept either a file path (and open it) or an existing `H5File` object (and use it, setting `owns = FALSE`).

**[x] Ticket 2: Correct Location of `matricized_access` for `LatentNeuroVec`**
    *   **Issue**: The `matricized_access` method for `LatentNeuroVec` (signature `i="integer"`) is incorrectly located in `R/reduced_clustered_vec.R`.
    *   **File(s)**: `R/reduced_clustered_vec.R`, `R/latent_vec.R`.
    *   **Recommendation**: Move the `setMethod("matricized_access", signature(x = "LatentNeuroVec", i = "integer"), ...)` definition from `R/reduced_clustered_vec.R` to `R/latent_vec.R`. Update `DESCRIPTION` Collate order if necessary.

**[-] Ticket 3: Review and Consolidate `.get_cluster_timeseries_by_mask_index` Helper (Obsolete due to H5ClusteredVec deprecation)**
    *   **Issue**: A deprecated version of `.get_cluster_timeseries_by_mask_index` exists in `R/clustered_vec.R`, while a more current version (using `.dataset_path` generic) is in `R/cluster_array.R`. The deprecated `H5ClusteredVec` methods might still use the old one.
    *   **File(s)**: `R/clustered_vec.R`, `R/cluster_array.R`.
    *   **Recommendation**:
        1.  Identify if any active (even if deprecated) code paths still call the old helper in `R/clustered_vec.R`.
        2.  If so, refactor them to use the modern helper from `R/cluster_array.R`.
        3.  Remove the outdated helper function from `R/clustered_vec.R` to prevent future misuse.

---

### 🟧 Medium Priority / Design & Modularity Improvements

**[x] Ticket 4: Clarify HDF5 Handle Lifecycle for `H5NeuroVol`/`Vec`/`LabeledVolumeSet`**
    *   **Issue**: Constructors like `H5NeuroVol(file_name)`, `H5NeuroVec(file_name)`, and `read_labeled_vec(file_path)` open HDF5 files but make the user responsible for calling `close(object)`. This needs very clear documentation.
    *   **File(s)**: `R/h5neurovol.R`, `R/h5neurovec.R`, `R/labeled_vec.R`, package documentation.
    *   **Recommendation**:
        1.  Add prominent documentation (e.g., in `@details` or a dedicated "Lifecycle" section in roxygen comments) for these classes and their constructors, explicitly stating the user's responsibility to `close()` objects created from file paths.
        2.  Briefly consider if `reg.finalizer` is a viable (though imperfect) safety net, but prioritize clear documentation.

**[x] Ticket 5: Centralize Duplicated Helper Functions**
    *   **Issue**:
        *   `%||%` operator is defined in `R/assertions.R` and `R/cluster_experiment.R`.
        *   `guess_h5_type` is defined in `R/h5_utils.R` and `R/read_write_hdf5.R`.
    *   **File(s)**: `R/assertions.R`, `R/cluster_experiment.R`, `R/h5_utils.R`, `R/read_write_hdf5.R`.
    *   **Recommendation**:
        1.  Move one definition of `%||%` (e.g., from `R/assertions.R`) to a central utility file (or keep in `R/assertions.R` if deemed most appropriate) and remove the other.
        2.  Move one definition of `guess_h5_type` (e.g., from `R/h5_utils.R`) to a central utility file and remove the other. Update all call sites.

**[x] Ticket 6: Simplify Mask Validation in `H5ClusteredExperiment` Constructor**
    *   **Issue**: The two-stage mask validation and the method for detecting if a mask was user-provided in the `H5ClusteredExperiment` constructor are complex and potentially fragile.
    *   **File(s)**: `R/cluster_experiment.R`.
    *   **Recommendation**:
        1.  Simplify the logic for determining if `mask` was originally non-NULL when the constructor was called.
        2.  Consider refactoring the stricter HDF5-based validation checks (currently inline) into `ensure_mask` or a new, clearly named helper function.

**[x] Ticket 7: Simplify `dtype` Argument in `write_labeled_vec`**
    *   **Issue**: `write_labeled_vec` accepts a list for `dtype` but errors if types in the list are not identical, as the NIfTI-like header supports only one datatype.
    *   **File(s)**: `R/labeled_vec.R`.
    *   **Recommendation**: Change the `dtype` argument to accept only a single `H5T` type object. Update documentation and internal logic accordingly.

**[x] Ticket 8: Validate Consistent Run Types for `/scans@summary_only` Attribute**
    *   **Issue**: The `/scans@summary_only` attribute in `write_clustered_experiment_h5` is set based on the `type` of the *first* run. This could be misleading if `runs_data` contains mixed run types.
    *   **File(s)**: `R/io_write_h5.R`.
    *   **Recommendation**: Add validation to `write_clustered_experiment_h5` to ensure all items in the `runs_data` list have the same `type` (either all "full" or all "summary").

---

### 🟦 Low Priority / Refinements & Observations

**[ ] Ticket 9: Review `obj = NULL` Prototype for HDF5 Slots**
    *   **Issue**: Classes like `H5NeuroVol` use `h5obj = NULL` in their prototype.
    *   **File(s)**: `R/all_class.R` (and relevant class definitions).
    *   **Recommendation**: Review if this pattern is ideal or if ensuring a valid (even if "closed" state) handle upon object creation is preferable. For now, ensure all methods robustly check for NULL and $is_valid.

**[x] Ticket 10: Verify DESCRIPTION Collate Order**
    *   **Issue**: The DESCRIPTION file has 'labeled_volume.R' which seems to be a typo for 'labeled_vec.R'.
    *   **File(s)**: `DESCRIPTION`.
    *   **Recommendation**: Correct the typo. Review the Collate order if Roxygen isn't managing it perfectly, especially given the S4 class dependencies.

**[x] Ticket 11: Robust nTime Inference in write_clustered_experiment_h5**
    *   **Issue**: For "full" data, nTime is inferred from ncol(sdata[[1]]).
    *   **File(s)**: `R/io_write_h5.R`.
    *   **Recommendation**:
        1.  If run$metadata$n_time is available, prioritize it.
        2.  If not, iterate through sdata (the list of cluster matrices) to find the first non-empty matrix to infer nTime, rather than relying solely on sdata[[1]].

**[x] Ticket 12: Ensure Full Cleanup of Deprecated Code Paths and Helpers**
    *   **Issue**: Deprecated classes/functions might have associated internal helpers that also become obsolete.
    *   **File(s)**: Primarily `R/cluster_vec_seq.R`, `R/clustered_vec.R`, `R/reduced_clustered_vec.R`.
    *   **Recommendation**: As part of phasing out deprecated functionality, ensure that their specific helper functions (e.g., the older .get_cluster_timeseries_by_mask_index) are removed or clearly marked if still essential for a very specific, documented transitional purpose.

**[x] Ticket 13: Enhance validate_latent_file Test Coverage**
    *   **Issue**: Test coverage for validate_latent_file could be improved.
    *   **File(s)**: `tests/testthat/test_latent_vec.R` (or a dedicated validation test file).
    *   **Recommendation**: Add specific tests that create intentionally malformed HDF5 files designed to trigger each of the different validation checks within validate_latent_file, ensuring it correctly identifies and reports these diverse issues.

**[x] Ticket 14: Investigate Redundant spacing in as_h5 for NeuroVol**
    *   **Issue**: In `R/h5neurovol.R`, `as_h5("NeuroVol", ...)` reconstructs `sp_read` using `NeuroSpace(dim=..., spacing=..., origin=..., trans=...)`. `neuroim2::NeuroSpace` might derive spacing from `trans` if both are provided.
    *   **File(s)**: `R/h5neurovol.R`.
    *   **Recommendation**: Briefly check `neuroim2::NeuroSpace` behavior: if `trans` is given, is spacing derived or does the explicit spacing argument override? Ensure the intended behavior. If spacing is redundant, it can be removed from the call.
    *   **Resolution Note**: The `as_h5` method explicitly writes `/space/spacing` and reads it back to pass to the `NeuroSpace` constructor, alongside `trans`. Debug messages are in place to compare `spacing(sp_read)` with `diag(trans(sp_read))[1:3]`. User should run and observe this output. If consistent, explicitly passing `spacing` might be redundant. If inconsistent, the `spacing` argument might be overriding `trans`.

**[x] Ticket 15: Ensure Graceful NULL Handling for clusters Slot in H5ClusteredRunSummary**
    *   **Issue**: H5ClusteredRunSummary objects might have a NULL @clusters slot.
    *   **File(s)**: `R/cluster_array.R` (methods for H5ClusteredRunSummary).
    *   **Recommendation**: Review all methods for H5ClusteredRunSummary (e.g., show, as.matrix, any others that might implicitly access @clusters) to ensure they gracefully handle cases where @clusters is NULL.
    *   **Resolution Note**: Methods for H5ClusteredRunSummary reviewed. `as.matrix` and `as.data.frame` do not use the `@clusters` (map) slot. Voxel-level accessors (`[`, `series`, `linear_access`) are disabled. The `show` method correctly checks for `NULL` and validity of `@clusters` before use. No changes needed.
