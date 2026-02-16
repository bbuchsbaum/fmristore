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
