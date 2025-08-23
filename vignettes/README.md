# fmristore Package Vignettes

This directory contains vignettes demonstrating the key functionality of the fmristore package.

## Quick Start

New to fmristore? Start with the **overview vignette** for a comprehensive introduction and guidance on choosing the right storage format for your data.

## Available Vignettes

### 1. fmristore-overview.Rmd ⭐ NEW
**fmristore: Efficient Storage for Neuroimaging Data**
- Comprehensive package overview and quick start guide
- Data conversion workflows from neuroim2 structures
- Converting ClusteredNeuroVec objects to HDF5 format
- Decision guide for choosing storage formats
- Performance comparisons and best practices
- Common use cases and troubleshooting

### 2. H5Neuro.Rmd
**Getting Started with HDF5 Neuroimaging Data**
- Introduction to HDF5-backed neuroimaging data structures
- Working with H5NeuroVol (3D volumes) and H5NeuroVec (4D time series)
- Memory-efficient data access patterns
- Basic operations and subsetting

### 3. LabeledVolumeSet.Rmd
**Working with Multiple Brain Maps using LabeledVolumeSet**
- Storing collections of related brain volumes
- Organizing statistical maps, contrasts, and group results
- Efficient access by label or index
- Best practices for data organization

### 4. H5ParcellatedMultiScan.Rmd
**Working with Clustered fMRI Data: H5ParcellatedMultiScan**
- Overview of clustered/parcellated fMRI data storage
- Understanding single scan vs multi-scan containers
- Converting ClusteredNeuroVec objects with as_multiscan parameter
- Creating and managing multi-run experiments
- Memory-efficient access to clustered time series
- Working with both full and summary data

### 5. H5ParcellatedScan.Rmd
**Understanding H5ParcellatedScan and H5ParcellatedScanSummary**
- Converting neuroim2 ClusteredNeuroVec objects
- Detailed comparison of full vs summary cluster storage
- Single scan creation and management
- When to use each type of representation
- Storage efficiency and memory considerations
- Practical examples with visualizations

## Building Vignettes

To build all vignettes:
```r
devtools::build_vignettes()
```

To build a specific vignette:
```r
rmarkdown::render("vignettes/H5Neuro.Rmd")
```

## Viewing Vignettes

After installation, view vignettes with:
```r
browseVignettes("fmristore")
```