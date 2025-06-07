# fmristore

<!-- badges: start -->
[![R-CMD-check](https://github.com/bbuchsbaum/fmristore/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bbuchsbaum/fmristore/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/bbuchsbaum/fmristore/branch/main/graph/badge.svg)](https://codecov.io/gh/bbuchsbaum/fmristore?branch=main)
[![Code Style](https://github.com/bbuchsbaum/fmristore/actions/workflows/style.yaml/badge.svg)](https://github.com/bbuchsbaum/fmristore/actions/workflows/style.yaml)
[![Lint](https://github.com/bbuchsbaum/fmristore/actions/workflows/lint.yaml/badge.svg)](https://github.com/bbuchsbaum/fmristore/actions/workflows/lint.yaml)
[![Dependencies](https://github.com/bbuchsbaum/fmristore/actions/workflows/dependencies.yaml/badge.svg)](https://github.com/bbuchsbaum/fmristore/actions/workflows/dependencies.yaml)
[![pkgcheck](https://github.com/bbuchsbaum/fmristore/actions/workflows/pkgcheck.yaml/badge.svg)](https://github.com/bbuchsbaum/fmristore/actions/workflows/pkgcheck.yaml)
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
devtools::install_github("bbuchsbaum/fmristore")
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

- [Package documentation](https://bbuchsbaum.github.io/fmristore/)
- [Getting started vignette](vignettes/README.md)
- [API reference](https://bbuchsbaum.github.io/fmristore/reference/)

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

This package strictly follows the [tidyverse style guide](https://style.tidyverse.org/). Style compliance is enforced through automated checks:

**Automated Style Enforcement:**
- GitHub Actions automatically check style on every push/PR
- Pre-commit hooks available for local development
- Configured with strict tidyverse linting rules

**Checking and Fixing Style Issues:**

```r
# Check for style violations
lintr::lint_package()

# Automatically fix style issues
styler::style_pkg()

# Check what styler would change without modifying files
styler::style_pkg(dry = "on")
```

**Setting Up Pre-commit Hooks (Recommended for Contributors):**

```bash
# Install pre-commit (requires Python)
pip install pre-commit

# Install the git hook scripts
pre-commit install

# (Optional) Run against all files
pre-commit run --all-files
```

The pre-commit hooks will automatically:
- Format code with `styler` according to tidyverse style
- Run `lintr` to catch style violations
- Check roxygen documentation
- Ensure DESCRIPTION file follows tidy format

## Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. (Recommended) Set up pre-commit hooks: `pre-commit install`
4. Make your changes
5. Ensure tests pass (`devtools::test()`)
6. Ensure code follows tidyverse style:
   ```r
   styler::style_pkg()      # Auto-fix style issues
   lintr::lint_package()    # Check for remaining violations
   ```
7. Commit your changes (`git commit -m 'Add amazing feature'`)
8. Push to the branch (`git push origin feature/amazing-feature`)
9. Open a Pull Request

**Style Requirements:**
- All code must pass `lintr::lint_package()` with no violations
- Code should be formatted with `styler::style_pkg()` using tidyverse style
- Pre-commit hooks will help ensure compliance automatically

## License

This project is licensed under the [license specified in DESCRIPTION file].

## Citation

To cite this package in publications, use:

```r
citation("fmristore")
```

## Getting Help

- [Issue tracker](https://github.com/bbuchsbaum/fmristore/issues) for bug reports and feature requests
- [Discussions](https://github.com/bbuchsbaum/fmristore/discussions) for questions and community support 