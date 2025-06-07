# Development Guide

This document outlines the development workflow, quality assurance tools, and continuous integration setup for the fmristore package.

## Continuous Integration & Quality Assurance

### GitHub Actions Workflows

The following automated workflows run on every push and pull request:

1. **R-CMD-check** (`.github/workflows/R-CMD-check.yaml`)
   - Tests package across multiple R versions (devel, release, oldrel-1)
   - Tests on multiple OS (Ubuntu, macOS, Windows)
   - Includes system dependency installation (HDF5, etc.)

2. **Test Coverage** (`.github/workflows/test-coverage.yaml`)
   - Runs comprehensive test suite with coverage reporting
   - Uploads results to Codecov
   - Configured with 80% target coverage

3. **Code Style** (`.github/workflows/style.yaml`)
   - Enforces tidyverse style guide using styler
   - Fails if code formatting is inconsistent
   - Provides clear instructions for fixing issues

4. **Linting** (`.github/workflows/lint.yaml`)
   - Runs lintr with tidyverse-focused rules
   - Catches style violations and potential issues
   - Configured for scientific R packages

5. **Dependencies** (`.github/workflows/dependencies.yaml`)
   - Weekly check for security vulnerabilities
   - License compliance checking
   - Outdated package detection

6. **Package Check** (`.github/workflows/pkgcheck.yaml`)
   - rOpenSci package standards compliance
   - Additional quality checks

7. **Documentation** (`.github/workflows/pkgdown.yaml`)
   - Builds and deploys package documentation website
   - Automatic updates on main branch changes

### Configuration Files

- **`.lintr`**: Linting configuration with tidyverse rules
- **`codecov.yml`**: Coverage reporting configuration
- **`.styler.R`**: Code formatting configuration
- **`.pre-commit-config.yaml`**: Pre-commit hooks for local development

## Development Workflow

### Setting Up Local Environment

1. **Install dependencies**:
   ```r
   devtools::install_dev_deps()
   ```

2. **Set up pre-commit hooks** (recommended):
   ```bash
   pip install pre-commit
   pre-commit install
   ```

### Code Quality Checks

Before committing changes, ensure:

1. **Tests pass**:
   ```r
   devtools::test()
   ```

2. **Code follows tidyverse style**:
   ```r
   styler::style_pkg()      # Auto-fix formatting
   lintr::lint_package()    # Check for violations
   ```

3. **Coverage is adequate**:
   ```r
   covr::package_coverage()
   ```

### Style Guide Enforcement

This package strictly enforces the [tidyverse style guide](https://style.tidyverse.org/):

- **Automated formatting**: `styler::style_pkg()` fixes most issues
- **Linting rules**: Comprehensive tidyverse-focused linting
- **CI enforcement**: Style checks block merging of non-compliant code
- **Pre-commit hooks**: Catch issues before commit

### Contributing Process

1. Fork and create feature branch
2. Make changes following tidyverse style
3. Ensure all tests pass and coverage is maintained
4. Run style checks and fix any violations
5. Commit and push (pre-commit hooks will verify quality)
6. Open pull request (CI will run all checks)

## Badge Status

The README includes badges showing the status of:
- R-CMD-check across multiple platforms
- Test coverage percentage
- Code style compliance
- Linting status
- Dependency security

## Package Build Configuration

- **`.Rbuildignore`**: Excludes CI/development files from package build
- **System dependencies**: Properly configured for HDF5 across platforms
- **Documentation**: Automated pkgdown site generation

## Troubleshooting

### Style Issues
```r
# See what styler would change
styler::style_pkg(dry = "on")

# Apply fixes
styler::style_pkg()

# Check remaining violations
lintr::lint_package()
```

### Test Failures
```r
# Run specific test file
devtools::test_file("tests/testthat/test-file.R")

# Run with debugging
devtools::test(filter = "test_name")

# Check coverage for specific file
covr::file_coverage("R/file.R", "tests/testthat/test-file.R")
```

### CI Failures
- Check GitHub Actions logs for specific failure details
- Most common issues: style violations, test failures, missing dependencies
- Local reproduction: ensure you can reproduce CI environment locally

## Quality Metrics

Current quality standards:
- **Test Coverage**: Target 80%, threshold ±1%
- **Patch Coverage**: Target 75%, threshold ±5%
- **Style Compliance**: 100% (enforced)
- **Linting**: Zero violations (enforced)
- **R-CMD-check**: Pass on all platforms/R versions 