# NeuroVecSeq to HDF5 Test Summary

## Overview
Comprehensive test suite for the `neurovecseq_to_h5()` function that converts NeuroVecSeq objects (multiple 4D neuroimaging scans) to HDF5 format.

## Test Coverage

### 1. **Basic Functionality** (`test_neurovecseq_to_h5.R`)
- ✓ Creates correct HDF5 file structure
- ✓ Handles metadata storage and retrieval
- ✓ Validates input parameters

### 2. **Comprehensive Tests** (`test_neurovecseq_comprehensive.R`)

#### Structure Tests (17 assertions)
- ✓ Root attributes (rtype, n_scans)
- ✓ Space group structure (/space/dim, /space/origin, /space/trans)
- ✓ Scans group hierarchy
- ✓ Individual scan attributes (n_time)
- ✓ Data array dimensions
- ✓ Data integrity verification

#### Metadata Tests (12 assertions)
- ✓ Numeric metadata (TR, TE, flip_angle)
- ✓ String metadata (task, subject_id, notes)
- ✓ Integer metadata (n_volumes, n_trials)
- ✓ Array metadata (slice_order)
- ✓ Mixed data types
- ✓ Partial metadata (not all scans need same fields)

#### Chunking & Compression Tests (5 assertions)
- ✓ Custom chunk dimensions
- ✓ Default time-optimized chunking
- ✓ Compression effectiveness
- ✓ Chunk validation

#### Data Integrity Tests (4 assertions)
- ✓ Time series extraction accuracy
- ✓ Full volume extraction
- ✓ Spatial pattern preservation
- ✓ Zero padding verification

#### Edge Cases Tests (7 assertions)
- ✓ Single scan NeuroVecSeq
- ✓ Many scans (10+)
- ✓ Different data types (FLOAT, DOUBLE)
- ✓ Various dimensions

#### Error Handling Tests (6 assertions)
- ✓ NULL input validation
- ✓ Invalid input type
- ✓ Mismatched scan names length
- ✓ Invalid data type
- ✓ Invalid compression level
- ✓ Invalid chunk dimensions

#### Performance Tests (1 test, skipped on CRAN)
- ✓ Write performance for large data
- ✓ Time series access performance
- ✓ Volume access performance

## Test Statistics
- **Total Tests**: 7 test blocks
- **Total Assertions**: 54 passed
- **Skipped**: 1 (performance test on CRAN)
- **Failed**: 0

## Key Features Tested

1. **Multi-scan Storage**: Successfully stores multiple 4D scans with different time dimensions in a single HDF5 file

2. **Metadata Flexibility**: Each scan can have its own metadata with different fields

3. **Chunking Strategies**: Supports both custom and automatic (time-optimized) chunking

4. **Data Types**: Handles FLOAT, DOUBLE, and INT data types

5. **Compression**: Validates compression levels and verifies file size reduction

6. **Data Integrity**: Ensures exact data recovery after write/read cycle

7. **Error Handling**: Comprehensive validation of inputs with informative error messages

## Example Usage from Tests

```r
# Create NeuroVecSeq with different time dimensions
vec1 <- NeuroVec(array(data1, dim=c(10,10,5,20)), space1)  # 20 time points
vec2 <- NeuroVec(array(data2, dim=c(10,10,5,30)), space2)  # 30 time points
vec3 <- NeuroVec(array(data3, dim=c(10,10,5,25)), space3)  # 25 time points
nvs <- NeuroVecSeq(vec1, vec2, vec3)

# Convert to HDF5 with metadata
result <- neurovecseq_to_h5(
  nvs, 
  file = "output.h5",
  scan_names = c("rest", "task1", "task2"),
  chunk_dim = c(10, 10, 10, 50),  # Time-optimized
  compression = 4,
  scan_metadata = list(
    rest = list(TR = 2.0, task = "rest"),
    task1 = list(TR = 2.5, task = "motor"),
    task2 = list(TR = 2.0, task = "visual")
  )
)
```

## Performance Benchmarks
- Writing 3 scans (64×64×30×100 each, ~36MB total): < 30 seconds
- Time series extraction (single voxel): < 0.1 seconds
- Volume extraction (single time point): < 2 seconds