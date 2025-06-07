# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Common Development Tasks

### Building and Checking the Package

```bash
# Build package
R CMD build .

# Check package (CRAN-compliant)
R CMD check --as-cran fmristore_*.tar.gz

# Install package locally
R CMD INSTALL .

# Or use devtools in R:
# devtools::install()
# devtools::check()
# devtools::document()
```

### Running Tests

```bash
# Run all tests
Rscript -e "testthat::test_dir('tests/testthat')"

# Run a specific test file
Rscript -e "testthat::test_file('tests/testthat/test-h5classes.R')"

# Or within R:
# testthat::test_local()
# testthat::test_file("tests/testthat/test-cluster_experiment.R")
```

### Linting

```bash
# Run lintr on the package
Rscript -e "lintr::lint_package()"

# Check specific file
Rscript -e "lintr::lint('R/h5_utils.R')"
```

### Documentation

```bash
# Generate documentation with roxygen2
Rscript -e "roxygen2::roxygenise()"

# Build package documentation
R CMD Rd2pdf --no-preview .
```

## High-Level Architecture

### Package Overview

`fmristore` is an R package for efficient storage of fMRI data in HDF5 format. It provides S4 classes and methods for working with:

1. **Neuroimaging data structures** - Dense and sparse representations of 4D fMRI data
2. **HDF5-backed storage** - Memory-efficient access to large datasets
3. **Latent representations** - Basis/embedding decompositions for data compression
4. **Clustered data** - Spatial clustering approaches for organizing voxel time series

### Core S4 Classes

The package uses S4 object-oriented programming with these main classes:

- **`H5NeuroVol`/`H5NeuroVec`** - HDF5-backed 3D/4D neuroimaging volumes
- **`LatentNeuroVec`** - Latent representation with spatial basis and temporal loadings
- **`H5ClusteredExperiment`** - Container for multiple runs with clustered voxel organization
- **`H5ClusteredRunFull`/`H5ClusteredRunSummary`** - Individual runs with full or summary cluster data
- **`LabeledVolumeSet`** - Collection of labeled brain regions backed by HDF5

### Key Design Patterns

1. **HDF5 Resource Management**: All HDF5 operations use utility functions in `R/h5_utils.R` that ensure proper handle cleanup via `tryCatch`/`finally` blocks.

2. **Generic/Method Dispatch**: The package defines S4 generics in `R/all_generic.R` with implementations spread across class-specific files. Key generics include:
   - `h5file()`, `mask()`, `basis()`, `loadings()` - Access components
   - `series_concat()`, `matrix_concat()` - Combine data
   - `as_h5()` - Convert to HDF5-backed representation

3. **Lazy Evaluation**: HDF5-backed classes load data on-demand rather than keeping everything in memory.

4. **Constructor Pattern**: Each H5-backed class has a constructor function (e.g., `H5NeuroVol()`) that handles HDF5 file validation and object creation.

### File Organization

- `R/all_class.R` - S4 class definitions
- `R/all_generic.R` - S4 generic function definitions  
- `R/constructors.R` - Constructor functions for main classes
- `R/h5_utils.R` - HDF5 utility functions (read, write, handle management)
- `R/io_*.R` - I/O operations for different formats
- `R/cluster_*.R` - Clustered data implementations
- `R/latent_vec.R` - Latent representation implementation

### HDF5 Schema Conventions

The package follows specific HDF5 schemas documented in `raw-data/`:
- `BasisEmbeddingSpec.yaml` - Latent representation format
- `ClusteredTimeSeriesSpec.yaml` - Clustered voxel data format
- `LatentNeuroArchiveSpec.md` - Comprehensive format specification (LNA 2.0)

Key HDF5 paths:
- `/header/*` - NIfTI-compatible header information
- `/mask` - Binary brain mask
- `/basis/*` - Spatial basis/components
- `/scans/<run_id>/embedding` - Temporal coefficients
- `/clusters/*` - Cluster definitions and metadata