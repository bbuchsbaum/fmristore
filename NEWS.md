# fmristore 0.1.0.9000 (Development)

## Breaking changes

* Removed the legacy surface that no package uses any more:
  `H5NeuroVol`, `H5NeuroVec`, `H5NeuroVecSeq`, `to_nih5_vec()`,
  `neurovecseq_to_h5()`, `LatentNeuroVec` I/O (`write_vec()`, `load_data()`
  for latent files, `validate_latent_file()`) and the `LatentNeuroVec`,
  `basis()`, `loadings()`, `offset()` and `map()` re-exports, `as_h5()`,
  `write_dataset()`, `write_cluster_result()`, `write_parcellated_scan_h5()`,
  the deprecated `write_clustered_experiment_h5()` and the deprecated
  `make_run_full()` and `make_run_summary()`. `read_dataset()` on those file
  types now fails with a message pointing to `write_frame_h5()`. fmristore no
  longer depends on fmrilatent.
* `H5ParcellatedMultiScan`, `H5ParcellatedScan`, `H5ParcellatedScanSummary`,
  `write_parcellated_experiment_h5()`, `LabeledVolumeSet` with
  `read_labeled_vec()` and `write_labeled_vec()`, `read_dataset()` and
  `detect_h5_type()` remain for neurotabs and fmrigds. They are frozen and
  will be removed once those packages store frames.
* Requires fmridataset >= 0.11.0.
* `h5vol[i]` on an `H5NeuroVol` became linear indexing before the class was
  removed; `write_parcellated_experiment_h5()` and the parcellated writers
  match full data to clusters by `cluster_<id>` name and validate coverage,
  row counts, column counts and `n_time`; `H5ParcellatedMultiScan(keep_handle_open = FALSE)`
  closes its handle when the constructor returns.
* The frame codec rejects integer values a declared dtype cannot store
  (out of range, fractional, or `NA` where the dtype has none) instead of
  silently converting them, and stores logicals with `NA` preserved.

## Fixes and improvements

* Closing a handle no longer closes every other handle open on the same file.
* Parcellated experiments return scans in the order they were written
  (`/scans@scan_order`).
* `H5ParcellatedMultiScan` gains `[[` (by scan name or position) and
  `names()`; `series_concat()` and `matrix_concat()` accept scan names.
* Summary scans label columns with the cluster IDs they hold.
* Parcellated headers store the real quaternion, offsets and `n_time`, so
  origin and orientation survive a round trip.
* `overwrite = TRUE` no longer deletes the existing file before validating
  the new input.
* Frame recovery works on healthy stores, only cleans directories with a
  valid root manifest, and finds the partials of interrupted shard appends.
  Zero-row and zero-feature frames can be written; read plans no longer
  overflow on large sources.
* Legacy `0.1-provisional` frames with free-form provenance migrate instead of
  failing.
* The test suite now runs under `R CMD check` (`tests/testthat.R` was
  missing), against current fmridataset.
* Replaced the vignettes with one on storing frames, `vignette("fmristore")`.

## Earlier development changes

* Added metadata-only frame-store inspection, conservative rollback of proven
  uncommitted write artifacts, and an atomic streaming migration from the real
  `0.1-provisional` HDF5 frame layout to FDS v1.

* Added atomic sharded-frame directories. New observation shards are committed
  as immutable FDS v1 HDF5 files before a small root-manifest swap, so appends
  never rewrite prior assay data and reopen as canonical row-sharded sources.

* HDF5 frame writers now support balanced, imagewise, and featurewise chunk
  layouts under explicit byte targets. `h5_read_plan()` and
  `validate_h5_read_amplification()` expose enforceable, metadata-only I/O
  amplification gates, including exact partial-edge chunk accounting.

* Added reconstructible `h5_array_source()` descriptors implementing the
  `fmridataset` lazy source protocol without serialized file handles.
* Replaced the provisional frame format with an atomic FDS v1 codec.
  `write_frame_h5()` now streams separately bound assays and aligned blocks
  under a hard memory budget; `open_frame_h5()` validates semantic and storage
  manifests while leaving every numerical array lazy.

# fmristore 0.1.0

* Initial CRAN release.

## Features

* HDF5-backed neuroimaging data structures:

  - `H5NeuroVol`: 3D brain volumes
  - `H5NeuroVec`: 4D time series data
  - `H5NeuroVecSeq`: Sequences of 4D scans

* Parcellated data storage:

  - `H5ParcellatedScan`: Single scan with cluster-based organization
  - `H5ParcellatedScanSummary`: Summary statistics per cluster
  - `H5ParcellatedMultiScan`: Multi-run experiments

* Latent representation I/O:

  - Read/write `LatentNeuroVec` objects to HDF5
  - Spec-compliant BasisEmbedding format

* Labeled volume sets:


  - `LabeledVolumeSet`: Named brain region collections
  - Efficient HDF5 storage with compression

* Utility functions:

  - `read_dataset()`: Auto-detect and load HDF5 neuroimaging data
  - `write_dataset()`: Generic writing interface
  - `as_h5()`: Convert in-memory objects to HDF5 format

## Dependencies

* Integrates with `neuroim2` for standard neuroimaging data structures
* Uses `hdf5r` for HDF5 file operations
* Requires `fmrilatent` for latent representation support
