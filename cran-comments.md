## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

* local macOS (aarch64-apple-darwin), R 4.4.x
* win-builder (devel and release)
* R-hub (various platforms)

## Downstream dependencies

This is a new package with no reverse dependencies.

## Notes

This package provides HDF5-backed storage for fMRI neuroimaging data,
integrating with the `neuroim2` package for standard data structures
and `fmrilatent` for latent space representations.

The package requires the HDF5 library, which is available on all major
platforms and is handled by the `hdf5r` package dependency.
