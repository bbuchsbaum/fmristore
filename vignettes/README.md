# fmristore Package Vignettes

This directory contains vignettes demonstrating the key functionality of the fmristore package.

## Available Vignettes

### 1. H5Neuro.Rmd
**Getting Started with HDF5 Neuroimaging Data**
- Introduction to HDF5-backed neuroimaging data structures
- Working with H5NeuroVol (3D volumes) and H5NeuroVec (4D time series)
- Memory-efficient data access patterns
- Basic operations and subsetting

### 2. LabeledVolumeSet.Rmd
**Working with Multiple Brain Maps using LabeledVolumeSet**
- Storing collections of related brain volumes
- Organizing statistical maps, contrasts, and group results
- Efficient access by label or index
- Best practices for data organization

### 3. H5ClusterExperiment.Rmd
**Working with Clustered fMRI Data: H5ClusterExperiment**
- Overview of clustered/parcellated fMRI data storage
- Creating and managing multi-run experiments
- Memory-efficient access to clustered time series
- Working with both full and summary data

### 4. H5ClusterRun.Rmd
**Understanding H5ClusterRun and H5ClusterRunSummary**
- Detailed comparison of full vs summary cluster storage
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