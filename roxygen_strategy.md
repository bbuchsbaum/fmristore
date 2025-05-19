You are absolutely right! CRAN is very strict about `\dontrun{}` and prefers examples that can be fully run and tested. The only general exception is for examples that would take too long (e.g., > 5 seconds), for which `\donttest{}` is the preferred wrapper. For examples that would write files or modify options, those need to be handled carefully, usually by writing to `tempdir()` and cleaning up, or resetting options.

Here's the revised strategy and ticketing system, incorporating the need for runnable examples and a strategy for test data.

## Overall Strategy for CRAN-Compliant Roxygen Documentation (Revised for Runnable Examples)

1.  **Identify All Exports:** Same as before.
2.  **Develop a Test Data Strategy:** This is crucial. We need small, self-contained, and quickly generated data/file fixtures that can be used in examples.
3.  **Prioritize S3 Generics/Methods with `@rdname`:** Same as before.
4.  **Address S4 Generics/Methods:** Same as before.
5.  **Document Standalone Exported Functions:** Same as before.
6.  **Document S4 Classes:** Same as before.
7.  **Iterative `R CMD check --as-cran`:** Same as before.
8.  **Content Focus:**
    *   `@return` specifies *what* is returned.
    *   `@examples` are *fully runnable* (using `tempdir()` for file outputs, creating minimal reproducible data on the fly, or using small, included datasets if appropriate). Use `\donttest{}` for examples that are correct but might exceed CRAN's typical time limits for examples (usually a few seconds). `\dontrun{}` should be avoided entirely.
    * We will proceed starting with all_class.R and all_generic.R. When a generic method is identified, you must locate it's s3 implementations, and add the to a special ticekt (that you will add to the end of this document, one per generic).
    * Once the generic ticket is defined, you will document the group as a unit as described in the strategy in this document.

---

### Phase 0: Test Data and Example Infrastructure

**[x] Ticket DOC-DATA-001: Design and Implement Test Data Generation Helpers**
    *   **Goal**: Create internal (non-exported) helper functions within the package (e.g., in `R/zzz_test_helpers.R` or a similar file not directly part of the core logic, or within `tests/testthat/helper-*.R` files if examples will call them via `fmristore:::helper_function_name`) to generate minimal, valid objects needed for examples.
    *   **Considerations for Helper Functions**:
        *   **Minimalism**: Create the smallest possible valid objects (e.g., a 2x2x2x1 `NeuroVol`, a mask with a few TRUE voxels, tiny basis/loading matrices).
        *   **Speed**: They must run very quickly.
        *   **Self-Contained**: They should not rely on external files unless those files are small and included in the package (e.g., in `inst/extdata`).
        *   **HDF5 Files**: For examples involving HDF5, helpers should create *temporary* HDF5 files in `tempdir()`. Examples will need to show this and ideally clean up, though `tempdir()` files are session-specific.
        *   **Object Types**: Create helpers for:
            *   `LogicalNeuroVol`
            *   `ClusteredNeuroVol`
            *   `DenseNeuroVec` (small 4D)
            *   `LatentNeuroVec` (minimal components)
            *   `H5ClusteredRunFull` (backed by a temp HDF5 with one cluster)
            *   `H5ClusteredRunSummary` (backed by a temp HDF5 with one cluster summary)
            *   `H5ClusteredExperiment` (backed by a temp HDF5 with one or two minimal runs)
            *   `LabeledVolumeSet` (backed by a temp HDF5 with one or two labels)
            *   `H5NeuroVol` / `H5NeuroVec` (backed by temp HDF5)
    *   **Action**:
        1.  Define requirements for each type of test object needed by example sections.
        2.  Implement helper functions. For HDF5-backed objects, these helpers will likely need to:
            *   Create a `tempfile(fileext = ".h5")` or `file.path(tempdir(), "example.h5")`.
            *   Write minimal valid data to this HDF5 file (e.g., using simplified versions of your `write_*` functions or `h5_write`).
            *   Return the path to the temp HDF5 file or the constructed R object that points to it.
        3.  Test these helper functions independently to ensure they are fast and correct.
    *   **Output**: A set of helper functions.
    *   **LLM Instruction**: "Develop a set of internal R helper functions (e.g., in `R/zzz_example_helpers.R`) for the `fmristore` package. These functions will generate minimal, valid instances of the core S4 classes (`LogicalNeuroVol`, `ClusteredNeuroVol`, `DenseNeuroVec`, `LatentNeuroVec`, `H5ClusteredRunFull`, `H5ClusteredRunSummary`, `H5ClusteredExperiment`, `LabeledVolumeSet`, `H5NeuroVol`, `H5NeuroVec`). For HDF5-backed objects, these helpers must create temporary HDF5 files in `tempdir()` and write the necessary minimal data, then return either the path or the R object referencing the temp file. The helpers must be fast and self-contained."

**[x] Ticket DOC-DATA-002: Define Standard Example Structure for HDF5 Operations**
    *   **Goal**: Establish a template for examples that involve creating and using HDF5 files.
    *   **Action**: The following template has been defined for examples involving HDF5 files:
        ```R
        # @examples

            
          # Setup: Create a temporary HDF5 file using a helper
          temp_h5_path <- NULL
          example_object <- NULL
          tryCatch({
            # Create temporary H5 file with minimal data
            temp_h5_path <- fmristore:::create_minimal_h5_for_example()
            
            # Load or create object using the temp file
            example_object <- fmristore::ConstructorFunction(file_name = temp_h5_path)
            
            # Demonstrate the function/method being documented
            result <- fmristore::FunctionBeingDocumented(example_object, other_params)
            print(result)
            
          }, error = function(e) {
            message("Example failed: ", e$message)
          }, finally = {
            # Clean up: Close HDF5 handles
            if (!is.null(example_object)) {
              try(close(example_object), silent = TRUE)  # For objects with a close method
              # OR for other types that need explicit handle closing:
              # if (inherits(example_object, "H5ClusteredExperiment") && !is.null(h5file(example_object))) {
              #   try(h5file(example_object)$close_all(), silent = TRUE)
              # }
            }
            
            # Delete temporary file
            if (!is.null(temp_h5_path) && file.exists(temp_h5_path)) {
              unlink(temp_h5_path)
            }
          })
        # } else { # This part was removed by the user
        #   message("Skipping example: required packages/functions not available.")
        # }
        ```
    *   **Key Features**:
        1. Uses `tryCatch` to handle errors properly and ensure cleanup.
        2. Uses helper functions to create minimal temporary HDF5 files.
        3. Shows proper resource management (closing handles, deleting temp files).
        4. Assumes all necessary dependencies are available for examples to run.
    *   **Output**: A documented template that ensures consistent, robust examples for HDF5 operations.

---

### Phase 1: Setup and Initial Identification (Same as before)

**[ ] Ticket DOC-001: Generate List of All Exported Objects** (No change)
**[ ] Ticket DOC-002: Establish `@rdname` Naming Conventions** (No change)

---

### Phase 2: Documenting S3 Generics and Methods (Revising `@examples`)

*(Repeat tickets DOC-S3-GEN and DOC-S3-METH for each S3 generic and its associated methods)*

**[ ] Ticket DOC-S3-GEN-Template: Document S3 Generic `[generic_name]` (Runnable Examples)**
    *   **Goal**: Create CRAN-compliant Roxygen documentation for S3 generic `[generic_name]`, with fully runnable examples.
    *   **File(s) Involved**: As before.
    *   **Action**:
        1.  (Steps 1-7, 9-10 from previous DOC-S3-GEN-Template are the same).
        2.  **Revise `@examples` (Step 8)**:
            *   Ensure all examples are runnable without `\dontrun{}`.
            *   Use helper functions (from DOC-DATA-001) to create necessary minimal test objects (e.g., a small data frame, a minimal `NeuroVol`).
            *   If the generic operates on HDF5-backed objects, follow the standard HDF5 example structure (from DOC-DATA-002).
            *   Use `\donttest{}` if an example, despite using minimal data, might still be too slow for CRAN's automated checks (target < 5 seconds).
    *   **LLM Instruction**: "For the S3 generic `[generic_name]`, update its Roxygen documentation. The `@examples` section must now be fully runnable. Use internal helper functions (e.g., `fmristore:::create_minimal_data_for_generic()`) to generate any necessary input objects. For examples involving HDF5, follow the standard pattern of creating temp HDF5 files and showing handle closure. Avoid `\dontrun{}`; use `\donttest{}` only if an example is computationally intensive even with minimal data."

**[ ] Ticket DOC-S3-METH-Template: Document S3 Method `[method_name].[class_name]` (Runnable Examples)**
    *   **Goal**: Document S3 method `[method_name].[class_name]`, ensuring examples (if any specific to this method) are runnable.
    *   **File(s) Involved**: As before.
    *   **Action**:
        1.  (Steps 1-6 from previous DOC-S3-METH-Template are the same).
        2.  **Revise `@examples` (Step 7)**: If this method has specific examples not covered by the generic:
            *   Ensure they are runnable, using helpers from DOC-DATA-001.
            *   Follow HDF5 example structure from DOC-DATA-002 if applicable.
            *   Avoid `\dontrun{}`; use `\donttest{}` only if essential.
    *   **LLM Instruction**: "For the S3 method `[method_name].[class_name]`, if you add specific `@examples` not covered by the generic, ensure they are fully runnable. Use internal helper functions and the standard HDF5 example pattern as needed. Avoid `\dontrun{}`."

---

### Phase 3: Documenting S4 Generics and Methods (Revising `@examples`)

*(Repeat tickets DOC-S4-GEN and DOC-S4-METH for each S4 generic and its methods)*

**[ ] Ticket DOC-S4-GEN-Template: Document S4 Generic `[generic_name]` (Runnable Examples)**
    *   **Goal**: Document S4 generic `[generic_name]`, with fully runnable examples for key signatures.
    *   **File(s) Involved**: As before.
    *   **Action**:
        1.  (Steps 1-6, 8-9 from previous DOC-S4-GEN-Template are the same).
        2.  **Revise `@examples` (Step 7)**:
            *   Ensure examples for common/important method signatures are runnable.
            *   Use helper functions (DOC-DATA-001) to create S4 objects of the required classes.
            *   Follow HDF5 example structure (DOC-DATA-002) when demonstrating methods on HDF5-backed classes.
            *   Avoid `\dontrun{}`; use `\donttest{}` only if necessary.
    *   **LLM Instruction**: "For the S4 generic `[generic_name]`, update its Roxygen documentation. The `@examples` section covering key method signatures must be fully runnable. Use internal helper functions to create necessary S4 objects. For examples involving HDF5, follow the standard pattern. Avoid `\dontrun{}`."

**[ ] Ticket DOC-S4-METH-Template: Review/Add Documentation for S4 Method `[generic_name],[signature]` (Runnable Examples)**
    *   **Goal**: Ensure S4 methods are correctly documented, and any method-specific examples are runnable.
    *   **File(s) Involved**: As before.
    *   **Action**:
        1.  (Steps 1-3 from previous DOC-S4-METH-Template are the same).
        2.  **Revise `@examples` (Step 4)**: If this method has specific examples:
            *   Make them runnable using helpers (DOC-DATA-001) and standard HDF5 structure (DOC-DATA-002).
            *   Avoid `\dontrun{}`.
    *   **LLM Instruction**: "For the `setMethod` call for `[generic_name]` with signature `[signature]`, if you add specific `@examples` not covered by the generic, ensure they are fully runnable. Use internal helper functions and the standard HDF5 example pattern as needed. Avoid `\dontrun{}`."

**[x] Ticket DOC-ALLGENERIC-S4-COMPLETE: Documented all exported S4 generics in R/all_generic.R with runnable examples (as of 2024-07-26)**
    *   **Goal**: Ensure all exported S4 generics in `R/all_generic.R` have CRAN-compliant, runnable examples.
    *   **Action**: Reviewed and updated examples for n_scans, scan_names, scan_metadata, cluster_metadata, h5file, basis, loadings, offset, mask, map, clusters, series_concat, matrix_concat, as_h5.
    *   **Status**: Completed.

**[x] Ticket DOC-S4-GEN-AS_H5: Document as_h5 Generic and Methods (as of 2024-07-27)**
    *   **Goal**: Document the `as_h5` generic function and its method implementations with CRAN-compliant, runnable examples.
    *   **File(s) Involved**: 
        - `R/all_generic.R` (generic definition) 
        - `R/io_h5_generic.R` (implementations)
    *   **Action(s)**: 
        - Used `@rdname as_h5-methods` approach to group method docs
        - Enhanced documentation of generic with comprehensive parameters, return values, and examples
        - Removed redundant docs from implementation methods
        - Added CRAN-compliant runnable examples with proper error handling and resource cleanup
    *   **Status**: Completed, tested examples are working.

**[x] Ticket DOC-S4-IMPORTED-WRITE_VEC: Fix write_vec Method Documentation (as of 2024-07-27)**
    *   **Goal**: Correctly document the `write_vec` method for LatentNeuroVec that implements the imported `write_vec` generic from neuroim2.
    *   **File(s) Involved**: 
        - `R/latent_vec.R` (write_vec method for LatentNeuroVec)
    *   **Action(s)**: 
        - Added proper `@importFrom neuroim2 write_vec` to specify the source of the generic
        - Replaced non-CRAN-compliant `\dontrun{}` example with a proper runnable example
        - Added proper error handling and resource cleanup in the example
        - Retained existing `@rdname write_vec-methods` for method grouping
    *   **Status**: Completed, example uses helper function to create minimal test data.

**[ ] Ticket DOC-FUNC-WRITE_CLUST_EXP: Document write_clustered_experiment_h5 (as of 2024-07-27)**
    *   **Goal**: Add a CRAN-compliant runnable example for the standalone function `write_clustered_experiment_h5`.
    *   **File(s) Involved**: 
        - `R/io_write_h5.R`
    *   **Action(s)**: 
        - Created a runnable example that uses `fmristore:::create_minimal_LogicalNeuroVol` and `fmristore:::create_minimal_ClusteredNeuroVol` helpers.
        - Manually constructed a minimal `runs_data` list with both "full" and "summary" type runs.
        - Included optional `cluster_metadata` generation.
        - Ensured the example writes to a temporary file and cleans up.
        - Incorporated standard checks for dependencies and error handling.
    *   **Status**: Pending review of example execution and CRAN compliance.

**[ ] Ticket DOC-S4-IMPORTED-LINEAR_ACCESS-HCRF: Document linear_access for H5ClusteredRunFull (as of 2024-07-27)**
    *   **Goal**: Add CRAN-compliant documentation and runnable example for the `linear_access` method for `H5ClusteredRunFull` objects, which implements an imported generic from `neuroim2`.
    *   **File(s) Involved**: 
        - `R/cluster_array.R` (linear_access method for H5ClusteredRunFull)
    *   **Action(s)**: 
        - Added `@rdname linear_access-methods` for documentation grouping.
        - Provided method-specific `@description`, `@param`, and `@return` tags.
        - Created a runnable `@examples` block using `fmristore:::create_minimal_h5_for_H5ClusteredExperiment` to set up an `H5ClusteredRunFull` object.
        - Ensured example includes dependency checks, error handling, and resource cleanup.
        - Verified that `@importFrom neuroim2 linear_access` is present at the file level.
    *   **Status**: Pending review of example execution and CRAN compliance.

---

### Phase 4: Documenting Standalone Exported Functions (Revising `@examples`)

*(Repeat ticket DOC-FUNC for each standalone exported function)*

**[ ] Ticket DOC-FUNC-Template: Document Standalone Function `[function_name]` (Runnable Examples)**
    *   **Goal**: Document standalone function `[function_name]` with fully runnable examples.
    *   **File(s) Involved**: As before.
    *   **Action**:
        1.  (Steps 1-6, 8 from previous DOC-FUNC-Template are the same).
        2.  **Revise `@examples` (Step 7)**:
            *   Ensure all examples are runnable.
            *   Use helpers (DOC-DATA-001) to create any complex input objects.
            *   Follow HDF5 example structure (DOC-DATA-002) if the function operates on or produces HDF5 files/objects.
            *   Avoid `\dontrun{}`.
    *   **LLM Instruction**: "For the exported function `[function_name]`, update its Roxygen documentation. The `@examples` section must be fully runnable. Use internal helper functions for complex inputs. If the function interacts with HDF5, follow the standard HDF5 example pattern. Avoid `\dontrun{}`."

---

### Phase 5: Documenting S4 Classes (Revising `@examples`)

*(Repeat ticket DOC-CLASS for each S4 class)*

**[ ] Ticket DOC-CLASS-Template: Document S4 Class `[ClassName]` (Runnable Examples)**
    *   **Goal**: Document S4 class `[ClassName]` with runnable examples for its constructor/`new()`.
    *   **File(s) Involved**: As before.
    *   **Action**:
        1.  (Steps 1-7, 9-11 from previous DOC-CLASS-Template are the same).
        2.  **Revise `@examples` (Step 8)**:
            *   Examples showing object creation (via `new()` or dedicated constructor like `H5ClusteredExperiment()`) must be runnable.
            *   Use helpers (DOC-DATA-001) if construction requires complex pre-existing objects or HDF5 files. Follow HDF5 example structure (DOC-DATA-002).
            *   Demonstrate basic usage of the created object.
            *   Avoid `\dontrun{}`.
    *   **LLM Instruction**: "For the S4 class `[ClassName]`, update its Roxygen documentation. The `@examples` section showing object creation (via `new()` or a constructor function) and basic usage must be fully runnable. Use internal helper functions if complex inputs or HDF5 files are needed for construction, following the standard HDF5 example pattern. Avoid `\dontrun{}`."

**[x] Ticket DOC-ALLCLASS-COMPLETE: Documented all S4 classes in R/all_class.R with runnable examples (as of 2024-07-26)**
    *   **Goal**: Ensure all S4 class definitions in `R/all_class.R` have CRAN-compliant, runnable examples.
    *   **Action**: Reviewed and updated examples for H5NeuroVol, H5NeuroVec, LatentNeuroVec, H5ClusteredVec (deprecated), H5ClusteredVecSeq (deprecated), H5ReducedClusteredVec (deprecated), H5ReducedClusteredVecSeq (deprecated), LabeledVolumeSet, H5ClusteredArray (virtual), H5ClusteredRunFull, H5ClusteredRunSummary, H5ClusteredExperiment.
    *   **Status**: Completed.

---

### Phase 6: Build, Check, and Refine (Same as before)

**[ ] Ticket DOC-BUILD-CHECK: Perform Full Package Build and CRAN Check** (No change in action, but the expectation is fewer `\dontrun` related errors and more focus on example correctness.)
**[ ] Ticket DOC-REVIEW-RD: Review Generated Rd Files** (No change)