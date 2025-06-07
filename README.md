# fmristore

<!-- badges: start -->
[![R-CMD-check](https://github.com/USER/fmristore/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/USER/fmristore/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/USER/fmristore/branch/main/graph/badge.svg)](https://codecov.io/gh/USER/fmristore?branch=main)
[![Lint](https://github.com/USER/fmristore/actions/workflows/lint.yaml/badge.svg)](https://github.com/USER/fmristore/actions/workflows/lint.yaml)
[![Dependencies](https://github.com/USER/fmristore/actions/workflows/dependencies.yaml/badge.svg)](https://github.com/USER/fmristore/actions/workflows/dependencies.yaml)
[![pkgcheck](https://github.com/USER/fmristore/actions/workflows/pkgcheck.yaml/badge.svg)](https://github.com/USER/fmristore/actions/workflows/pkgcheck.yaml)
<!-- badges: end -->

> **Note**: Replace `USER` with your actual GitHub username in the badge URLs above.

An R package for efficient storage and retrieval of neuroimaging data using HDF5 format with specialized data structures for fMRI analysis.

## Features

- **Efficient Storage**: HDF5-based storage for large neuroimaging datasets
- **Memory-Optimized**: Latent vector representations for reduced memory footprint
- **High Performance**: Optimized data access patterns for fMRI analysis workflows
- **Rich Data Structures**: Support for clustered, dense, and sparse neuroimaging vectors
- **Comprehensive Testing**: Extensive test suite ensuring reliability and correctness

## Installation

### From GitHub

```r
# Install development version from GitHub
# install.packages("devtools")
devtools::install_github("USER/fmristore")
```

### Dependencies

The package requires HDF5 system libraries. Installation instructions by platform:

#### Ubuntu/Debian
```bash
sudo apt-get install libhdf5-dev
```

#### macOS
```bash
# Using Homebrew
brew install hdf5

# Using MacPorts
sudo port install hdf5
```

#### Windows
HDF5 libraries are automatically handled during R package installation.

## Quick Start

```r
library(fmristore)

# Load neuroimaging data
# Example workflow here...
```

## Documentation

- [Package documentation](https://USER.github.io/fmristore/)
- [Getting started vignette](vignettes/README.md)
- [API reference](https://USER.github.io/fmristore/reference/)

## Development

### Running Tests

```r
# Install test dependencies
devtools::install_dev_deps()

# Run tests
devtools::test()

# Check test coverage
covr::package_coverage()
```

### Code Style

This package follows the [tidyverse style guide](https://style.tidyverse.org/). Use the following to check and fix style issues:

```r
# Check code style
lintr::lint_package()

# Auto-format code (if using styler)
styler::style_pkg()
```

## Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Make your changes
4. Ensure tests pass (`devtools::test()`)
5. Check code style (`lintr::lint_package()`)
6. Commit your changes (`git commit -m 'Add amazing feature'`)
7. Push to the branch (`git push origin feature/amazing-feature`)
8. Open a Pull Request

## License

This project is licensed under the [license specified in DESCRIPTION file].

## Citation

To cite this package in publications, use:

```r
citation("fmristore")
```

## Getting Help

- [Issue tracker](https://github.com/USER/fmristore/issues) for bug reports and feature requests
- [Discussions](https://github.com/USER/fmristore/discussions) for questions and community support 