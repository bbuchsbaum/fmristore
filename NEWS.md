# fmristore 0.1.0.9000 (Development)

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
