# fmristore 0.2.0

fmristore now does one thing: store `fmridataset::fmri_frame` objects in
HDF5. Its only imports are fmridataset and hdf5r.

## Breaking changes

* Removed every pre-frame class and reader: `H5ParcellatedMultiScan`,
  `H5ParcellatedScan`, `H5ParcellatedScanSummary`,
  `write_parcellated_experiment_h5()`, `LabeledVolumeSet` with
  `read_labeled_vec()` and `write_labeled_vec()`, `read_dataset()`,
  `detect_h5_type()`, `H5NeuroVol`, `H5NeuroVec`, `H5NeuroVecSeq`,
  `to_nih5_vec()`, `neurovecseq_to_h5()`, `LatentNeuroVec` I/O, `as_h5()`,
  `write_dataset()`, `write_cluster_result()`, `write_parcellated_scan_h5()`,
  `make_run_full()`, `make_run_summary()`, the `H5_PATHS`, `H5_ATTRS` and
  `H5_DSETS` constants, and the memory-warning settings. Store data as frames
  with `write_frame_h5()`; parcellations are a `fmridataset::parcel_space()`
  and labeled volumes an observation x voxel frame. Packages that still need
  the old classes can pin fmristore to commit `4d7adb7`.
* No longer imports neuroim2, fmrilatent, Matrix, assertthat, cli, lifecycle
  or withr.
* Requires fmridataset >= 0.11.0.
* The frame codec rejects integer values a declared dtype cannot store
  (out of range, fractional, or `NA` where the dtype has none) instead of
  silently converting them, and stores logicals with `NA` preserved.

## Fixes and improvements

* Closing a handle no longer closes every other handle open on the same file.
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
