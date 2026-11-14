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
Rscript -e "testthat::test_file('tests/testthat/test-frame-h5.R')"

# Or within R:
# testthat::test_local()
# testthat::test_file("tests/testthat/test-frame-h5.R")
```

### Linting

```bash
# Run lintr on the package
Rscript -e "lintr::lint_package()"

# Check specific file
Rscript -e "lintr::lint('R/frame_h5.R')"
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

`fmristore` stores `fmridataset::fmri_frame` objects in HDF5 (the "FDS v1"
frame store): atomic single-file snapshots, append-only sharded stores,
recovery and migration, and lazy HDF5 array sources. `fmridataset` owns frame
semantics; `fmristore` owns the physical storage. See
`design/fmristore-1.0-plan.md` and `design/deficiency-remediation-plan.md`.

### Frame store

- `R/frame_h5.R` - `write_frame_h5()` / `open_frame_h5()` (atomic write, lazy open)
- `R/frame_shards_h5.R` - `write_sharded_frame_h5()`, `append_frame_shard_h5()`, `open_sharded_frame_h5()`
- `R/frame_recovery.R` - `frame_store_info()`, `recover_frame_store()`, `migrate_frame_h5()`
- `R/array_source.R` - `h5_array_source()`, `h5_read_plan()`, `validate_h5_read_amplification()`
- `R/utils.R` - `close_h5_safely()` and `%||%`

The only imports are fmridataset and hdf5r. The legacy parcellated, labeled
and `H5NeuroVec` classes were removed in 0.2.0; consumers such as neurotabs and fmrigds
store frames instead. Do not reintroduce neuroimaging classes here: frame semantics
belong to fmridataset.

### Key Design Patterns

1. **HDF5 resource management**: close only the handle you own with
   `close_h5_safely()`. Never call `$close_all()` on a handle another object
   might share: HDF5 shares one file structure across handles on a path.
2. **Testing**: run tests against the installed package with `NOT_CRAN=true`
   and `_R_CHECK_MATRIX_DATA_=TRUE`, as CI does.
3. **Docs**: regenerate with roxygen2 at the version in `RoxygenNote`.

### HDF5 layout

The FDS v1 layout (manifests, assay arrays, schema version and migrations) is
defined in `R/frame_h5.R` and `design/fmristore-1.0-plan.md`.

## Issue Tracking

This project uses **beads** for issue tracking. See AGENTS.md for full workflow.

```bash
bd list                   # List open issues
bd ready                  # Find tasks without blockers
bd create "title" -p 1    # Create issue (P0-P3 priority)
bd show <id>              # View issue details
bd close <id>             # Close issue
bd sync                   # Sync to git
```

## Session Management

When ending a work session, always "land the plane":

1. File issues for remaining work
2. Run quality gates (tests, linting)
3. Update issue status
4. `bd sync && git push`
5. Provide handoff context for next session