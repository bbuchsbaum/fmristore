# TODO List for LatentNeuroVec Implementation and HDF5 Spec Alignment (Revised 2)

This list tracks fixes and improvements needed for `LatentNeuroVec` based on review against `BasisEmbeddingSpec.yaml`. Items marked [x] are considered complete *unless* subsequent feedback requires revision.

## I. HDF5 Spec Alignment (I/O Functions)

These items focus on making `write_vec`/`to_h5_latentvec` and `load_data`/`LatentNeuroVecSource` compliant with `BasisEmbeddingSpec.yaml`.

- [ ] **Write NIFTI Header:** Modify `to_h5_latentvec` to create `/header`. Populate ALL required NIfTI-1 fields, zero-filling if necessary (`dim_info`, `xyzt_units`, `slice_*`, etc.). **Notes:**
    - Include `magic`, `datatype`, `bitpix`.
    - Ensure `pixdim[4]` (TR) is populated correctly (from `NeuroSpace` attr or argument).
- [ ] **Write 3D Mask:** Modify `to_h5_latentvec` to write `/mask` (uint8, [X,Y,Z]). **Note:** Use default chunking `c(32,32,32)` and gzip.
- [ ] **Write Voxel Coords:** Modify `to_h5_latentvec` to calculate and write `/voxel_coords` ([nVox, 3], int32) from `which(mask == 1)`, ensuring order matches basis columns. Read it back in `load_data`.
- [ ] **Write Basis Matrix (Spec Compliant):** Modify `to_h5_latentvec`: create `/basis`, write `t(object@loadings)` to `/basis/basis_matrix` ([k, nVox]). **Note:** Ensure default chunking `c(k, chunk_size)` & gzip.
- [ ] **Write Offset (Spec Compliant):** Modify `to_h5_latentvec` to write `object@offset` to `/offset` ([nVox]).
- [ ] **Write Embeddings (Multi-Scan):** Modify `to_h5_latentvec`: create `/scans/<scan_name>/`, write `object@basis` to `/scans/<scan_name>/embedding` ([runLength, k]). **Note:** Use default chunking `c(chunk_size, k)` & gzip.
- [ ] **Write Scan Metadata:** Modify `to_h5_latentvec`: write metadata (`run_length`, `TR`) into `/scans/<scan_name>/metadata`. Write placeholders (empty strings/NA) for `subject_id`, `task`, `session`.
- [ ] **Write Root Attributes:** Ensure `to_h5_latentvec` writes `latent_spec_version="1.0"` and `rtype="LatentNeuroVec"` attributes.
- [ ] **Read NIFTI Header:** Modify `load_data` to read fields from `/header` to reconstruct `NeuroSpace`.
- [ ] **Read 3D Mask:** Modify `load_data` to read `/mask` dataset.
- [ ] **Read Voxel Coords:** Modify `load_data` to read `/voxel_coords`.
- [ ] **Read Basis Matrix (Spec Compliant):** Modify `load_data` to read `/basis/basis_matrix` ([k, nVox]) and assign `t()` to `@loadings`.
- [ ] **Read Offset (Spec Compliant):** Modify `load_data` to read `/offset`, default to zeros if absent.
- [ ] **Read Embeddings (Multi-Scan):** Modify `load_data` to load single specified scan via `scan_name`, reading `/scans/<scan_name>/embedding` into `@basis`.
- [ ] **Support `basis_reference` (Low Priority):** Modify I/O to optionally handle `basis_reference` structure.

## II. Code Bugs & Omissions

- [x] **B1: Add Missing Slots:** Added `map` and `label` to `setClass` with prototypes.
- [x] **B2: Remove `browser()`:** Removed from `linear_access`.
- [x] **B3: Guard Matricized Access (Matrix):** Added `ncol` check in `matricized_access,matrix-method`.
- [x] **Add Matricized Access Guard (Integer):** Added `stopifnot(ncol(x@basis)==ncol(x@loadings))` check to `matricized_access,integer-method`.
- [x] **B5: Close HDF5 Datasets:** Added `on.exit` for dataset handles in `to_h5_latentvec`.
- [x] **Rewrite `linear_access,LatentNeuroVec`:** Refactored `linear_access` for efficiency.
- [ ] **Fix Dataset Property Creation:** Correctly create and use the dataset property list object (`dtype_prop`) in `to_h5_latentvec` for setting chunking/compression.
- [ ] **Fix `linear_access` Recursion:** Modify `linear_access, i="numeric"` method to use `callNextMethod(x = x, i = as.integer(i))`.
- [ ] **Fix `space_for_map` Definition:** Initialize `space_for_map` before the `if` block in `LatentNeuroVec` constructor to avoid undefined variable error when mask is an array.
- [x] **Fix `%||%` Operator Usage:** Replaced with base R `if(is.null(...))`.
- [x] **Fix `IndexLookupVol` Space:** Ensured `IndexLookupVol` uses correct space, warns on mismatch.
- [ ] **Remove `write_vec` Double Close:** Remove explicit `obj$close()` call in `write_vec` method.

## III. Performance

- [ ] **Vectorize `[]` accessor:** Refactor `[` method using matrix math and `array()` reshape to avoid looping through time slices.
- [x] **Vectorize `series()` method:** Refactored using matrix operations.
- [ ] **Optimize `linear_access` Block:** Re-evaluate `linear_access` implementation to avoid creating potentially huge intermediate `data_block` if memory becomes an issue.
- [ ] **Add Sparse Coercion Checks:** Add `if (!is(matrix_obj, "Matrix"))` guards before `Matrix::Matrix()` coercion.

## IV. Naming & Clarity (Optional)

- [ ] **Consider Internal Slot Renaming:** Evaluate renaming internal slots `@basis` -> `@embedding` and `@loadings` -> `@spatial_basis`.

## V. Cross-Cutting Tasks

- [x] **File Version Attribute:** Implemented writing/reading `latent_spec_version` attribute.
- [x] **Validator Function:** Created `validate_latent_file()`.
- [ ] **Refine Validator:** Update `validate_latent_file()` to check `/voxel_coords` length vs basis columns.
- [x] **Add `setValidity` Method:** Implemented S4 `setValidity` for `LatentNeuroVec`.
- [ ] **Add Unit Tests:** Create comprehensive tests, including round-trip tests and validator tests. 