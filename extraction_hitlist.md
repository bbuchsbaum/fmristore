# Refactoring Hit List for fmristore Boilerplate Extraction

This list outlines the proposed refactoring tasks to extract common boilerplate code into helper functions, improving readability and maintainability without changing behaviour.

**Plan for Identifying Patterns:**

*   Use `grep` for specific code patterns (`H5File$new`, `on.exit`, `is.null(mask)`, `read_h5_mask...`, `h5obj$exists`, `$read()`, `$create_dataset`, `identical(dim(`).
*   Use semantic search if necessary for more complex patterns like the assertion guards.
*   Target files mentioned in the suggestions and relevant test files.

**Refactoring Tickets:**

- [X] **1. HDF5 Open/Close Boilerplate:**
    - [X] Create `R/h5_utils.R` file.
    - [X] Implement `open_h5(src, mode = "r")` helper function in `R/h5_utils.R`.
    - [X] Replace manual `H5File$new`, `file.exists` checks, and `on.exit` logic with `open_h5()` and `withr::defer(if (fh$owns) fh$h5$close_all())` in:
        - [X] `R/cluster_array.R` (make_run_full, make_run_summary)
        - [X] `R/cluster_experiment.R` (H5ClusteredExperiment constructor)
        - [ ] `R/h5neurovol.R` (*Not found by search*)
        - [ ] `R/h5neurovec.R` (*Not found by search*)
        - [X] `R/latent_vec.R` (load_data, validate_latent_file)
        - [X] `R/labeled_volume.R` (read_labeled_vec)
        - [X] `R/read_write_hdf5.R` (read_clustered_dataset)
        - [ ] Relevant test files (identify via search)

- [ ] **2. Mask Load/Validate Idiom:**
    - [X] Implement `ensure_mask(mask, h5, space, path = "/mask")` helper function in `R/h5_utils.R` (placeholder for dependency `read_h5_mask_to_LogicalNeuroVol` added).
    - [ ] Replace the 20-line `if (is.null(mask)) ... else ...` block with a call to `ensure_mask()` in:
        - [ ] `R/cluster_array.R` (*Pattern not found*)
        - [X] `R/cluster_experiment.R` (Stricter validation for user-provided masks kept separate for now)
        - [ ] `R/clustered_vec.R` (*Pattern not found*)
        - [ ] `R/reduced_clustered_vec.R` (*Pattern not found*)

- [X] **3. Assertion Wrappers:**
    - [X] Create `R/assertions.R` file.
    - [X] Implement `assert_non_empty_numeric(x, arg, fn)` in `R/assertions.R`.
    - [X] Identify and replace `is.null / !is.numeric / length == 0L` guards (approx. 16 places) with calls to `assert_non_empty_numeric()`.
        - [X] `R/cluster_experiment.R` (`series_concat`)
        - [ ] *(Other instances not found with specific grep/semantic search; further review potentially needed but deferred)*

- [X] **4. Read/Write Dataset Helpers:**
    - [X] Implement `h5_read(h5, path, missing_ok = FALSE)` in `R/h5_utils.R`.
    - [X] Implement `h5_write(h5, path, data, ..., missing_ok = FALSE)` in `R/h5_utils.R`.
    - [X] Replace manual `$exists()`, `$read()`, `$create_dataset()`, etc., calls with `h5_read()` and `h5_write()`.
        - [X] `R/labeled_volume.R` (read_labeled_vec, write_labeled_vec)
        - [X] `R/latent_vec.R` (load_data, validate_latent_file, to_h5_latentvec)
        - [X] `R/read_write_hdf5.R` (write_clustered_dataset)
        - [X] `R/io_write_h5.R` (write_clustered_experiment_h5)
        - [ ] `R/h5neurovol.R` (needs review)
        - [ ] `R/h5neurovec.R` (needs review)
        - [ ] *(Further review needed, deferred)*

- [X] **5. Shared Dimension Check Utility:**
    - [X] Implement `check_same_dims(a, b, dims_to_compare = NULL, msg = NULL)` in `R/assertions.R`.
    - [ ] Replace manual `identical(dim(...))` checks (approx. 6 places) with calls to `check_same_dims()`.
        - [X] `R/h5_utils.R` (ensure_mask)
        - [X] `R/cluster_array.R` (make_run_full, make_run_summary)
        - [X] `R/io_write_h5.R` (write_clustered_experiment_h5)
        - [ ] `R/latent_vec.R` (LatentNeuroVec constructor, .validate_LatentNeuroVec)
        - [ ] *(Identify remaining locations)*

- [ ] **6. Add `withr` Dependency:** Add `withr` to `Imports` in DESCRIPTION if `defer` is used.
- [ ] **7. Testing:** Run tests (`devtools::test()`) to ensure no behavior changed.

**Post-Refactoring Steps:**

- [ ] **8. Styling:**
    - [ ] Run `styler::style_pkg()` after all mechanical edits are complete. 

---

# API Consistency Sprint Tickets

- [x] **T-01: Introduce canonical constructors for clustered classes**
    - [x] Create `R/constructors.R` (or similar file).
    - [x] Implement `H5ClusteredRunFull()` factory function wrapping `new()`, `open_h5()`, `withr::defer()`.
    - [x] Implement `H5ClusteredRunSummary()` factory function wrapping `new()`, `open_h5()`, `withr::defer()`.
    - [x] Mark `make_run_full()` with `lifecycle::deprecate_warn()` pointing to `H5ClusteredRunFull()`.
    - [x] Mark `make_run_summary()` with `lifecycle::deprecate_warn()` pointing to `H5ClusteredRunSummary()`.
    - [x] Update unit tests (`test-cluster_run_full.R`, `test-cluster_run_summary.R`) to use `H5ClusteredRunFull()` / `H5ClusteredRunSummary()` factories.

- [ ] **T-02: Standardise parameter names**
    - [ ] Replace `file_source`, `file_name`, `fname` arguments with `file` across all relevant functions.
    - [ ] Harmonise `mask_idx` / `voxel_idx` parameter names (choose one standard name).
    - [ ] Harmonise `compression` / `compressed` parameter names (choose one standard name).
    - [ ] Verify removal: Run `grep -R '(file_source|file_name|fname)' R/` and ensure zero hits.

- [ ] **T-03: Implement as_h5() generic and methods**
    - [ ] Create `R/io_h5_generic.R` file.
    - [ ] Define `setGeneric("as_h5", ...)` in `R/io_h5_generic.R`.
    - [ ] Implement `setMethod("as_h5", "NeuroVol", ...)` to replace `to_nih5_vol()`.
    - [ ] Implement `setMethod("as_h5", "NeuroVec", ...)` to replace `to_nih5_vec()`.
    - [ ] Implement `setMethod("as_h5", "H5ClusteredRunFull", ...)` (if needed, or handle via internal dispatch).
    - [ ] Implement `setMethod("as_h5", "H5ClusteredRunSummary", ...)` (if needed, or handle via internal dispatch).
    - [ ] Implement `setMethod("as_h5", "H5ClusteredExperiment", ...)` (if needed, or handle via internal dispatch).
    - [ ] Implement `setMethod("as_h5", "LatentNeuroVec", ...)` to replace `to_h5_latentvec()`.
    - [ ] Implement `setMethod("as_h5", "LabeledVolume", ...)` to replace `write_labeled_vec()`.
    - [ ] Implement `setMethod("as_h5", "ClusteredDataset", ...)` to replace `write_clustered_dataset()`.
    - [ ] Mark `to_nih5_vol()` deprecated, forward to `as_h5()`.
    - [ ] Mark `to_nih5_vec()` deprecated, forward to `as_h5()`.
    - [ ] Mark `write_clustered_dataset()` deprecated, forward to `as_h5()`.
    - [ ] Mark `write_labeled_vec()` deprecated, forward to `as_h5()`.
    - [ ] Mark `to_h5_latentvec()` deprecated, forward to `as_h5()`.
    - [ ] Add unit tests for round-trip `object -> as_h5(tmpfile) -> reload -> identical(object, reloaded)` for relevant classes.
    - [ ] Verify deprecation warnings appear correctly during tests/checks.

- [ ] **T-04: Remove duplicate helper definitions**
    - [ ] Define canonical `%||%` helper in `R/utils.R` (or `R/h5_utils.R`) and export/import appropriately.
    - [ ] Remove duplicate `%||%` definition from `R/assertions.R`.
    - [ ] Remove duplicate `%||%` definition from `R/cluster_experiment.R`.
    - [ ] Define canonical `is_h5file()` helper in `R/utils.R` (or `R/h5_utils.R`) and export/import appropriately.
    - [ ] Remove duplicate `is_h5file()` definition from `R/cluster_array.R`.
    - [ ] Verify no duplicate definitions remain (`grep` for function definition).

- [ ] **T-05: Harmonise slot & attribute names (deep audit)**
    - [ ] Create `docs/naming_schema.md` defining canonical slot/attribute names (e.g., `n_time`, `compressed`, `space`, `mask`, `clusters`).
    - [ ] Audit `setClass` definitions in `R/all_class.R` (and potentially others) and refactor slots to use canonical names.
    - [ ] Update `as_h5` methods (and any remaining writers) to write HDF5 attributes matching the canonical slot names.
    - [ ] Implement migration shims in reading functions/constructors (e.g., `H5ClusteredExperiment`, `read_labeled_vec`, `load_data` for LatentVec) to recognize and map old HDF5 attribute names to new canonical slot names.
    - [ ] Add regression tests specifically for loading HDF5 files saved with older attribute naming conventions.
    - [ ] Verify name changes (e.g., `grep -R '(compression\b|compressed_flag)' R/` should only show `compressed` or similar target names).

- [ ] **T-06: Promote generic as_matrix() (Backlog)**
    - [ ] Define `as_matrix()` S3/S4 generic.
    - [ ] Implement methods for relevant vector/clustered classes.
    - [ ] Deprecate/remove older `series_concat()`, `matricized_access()`, etc.

- [ ] **T-07: Update documentation & pkgdown (Backlog)**
    - [ ] Revise Roxygen documentation for changed constructors, functions, and parameters.
    - [ ] Regenerate package documentation (Rd files).
    - [ ] Run `pkgdown::build_site()`. 