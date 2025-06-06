# Example: Using fmristore to store multiple 4D neuroimaging scans without clustering

library(fmristore)
library(neuroim2)

# Scenario: We have 3 fMRI runs from an experiment
# Each has the same spatial dimensions but different numbers of time points
# We want to store them in a single HDF5 file for efficient access

# Create example data - in practice these would be loaded from NIfTI files
# Run 1: Resting state scan (200 time points)
run1_data <- array(rnorm(64*64*30*200), dim = c(64, 64, 30, 200))
run1_space <- NeuroSpace(c(64, 64, 30, 200), spacing = c(3, 3, 3, 2))  # 3mm voxels, TR=2s
run1 <- NeuroVec(run1_data, run1_space)

# Run 2: Task scan (150 time points)
run2_data <- array(rnorm(64*64*30*150), dim = c(64, 64, 30, 150))
run2_space <- NeuroSpace(c(64, 64, 30, 150), spacing = c(3, 3, 3, 2))
run2 <- NeuroVec(run2_data, run2_space)

# Run 3: Another task scan (180 time points)
run3_data <- array(rnorm(64*64*30*180), dim = c(64, 64, 30, 180))
run3_space <- NeuroSpace(c(64, 64, 30, 180), spacing = c(3, 3, 3, 2))
run3 <- NeuroVec(run3_data, run3_space)

# Create a NeuroVecSeq to hold all runs
all_runs <- NeuroVecSeq(run1, run2, run3)

# Define metadata for each scan
scan_metadata <- list(
  rest = list(
    subject_id = "sub-01",
    session = "ses-01", 
    task = "rest",
    TR = 2.0,
    date = "2024-01-15"
  ),
  task_motor = list(
    subject_id = "sub-01",
    session = "ses-01",
    task = "motor",
    TR = 2.0,
    n_trials = 40,
    date = "2024-01-15"
  ),
  task_visual = list(
    subject_id = "sub-01", 
    session = "ses-01",
    task = "visual",
    TR = 2.0,
    n_trials = 60,
    date = "2024-01-15"
  )
)

# Convert to HDF5 with time-series optimized chunking
h5_file <- neurovecseq_to_h5(
  all_runs,
  file = "subject01_all_runs.h5",
  scan_names = c("rest", "task_motor", "task_visual"),
  data_type = "FLOAT",
  chunk_dim = c(16, 16, 10, 200),  # Optimize for extracting full time series
  compression = 4,
  scan_metadata = scan_metadata
)

cat("Created HDF5 file:", h5_file, "\n")

# The resulting file structure allows for:
# 1. Lazy loading - data is only read when accessed
# 2. Efficient time-series extraction for individual voxels
# 3. Metadata storage for each scan
# 4. Shared spatial information across all scans

# Example of how you might access the data later:
# h5f <- H5File$new(h5_file, mode = "r")
# rest_data <- h5f[["scans/rest/data"]][]  # Load full rest scan
# voxel_ts <- h5f[["scans/task_motor/data"]][30, 30, 15, ]  # Get time series for one voxel
# h5f$close_all()

# Clean up
unlink(h5_file)