# Data Loader for Spatial Variable Genes Benchmark

This module loads spatial transcriptomics datasets for the omnibenchmark framework.

## Overview

The data loader reads pre-processed h5ad files from the `datasets/` directory and makes them available to the benchmark pipeline. It handles three types of files:

- `dataset.h5ad` - The input dataset (required)
- `solution.h5ad` - Ground truth spatial variability scores (optional)
- `simulated_dataset.h5ad` - Simulated dataset with normalized data (optional)

## Usage

```bash
python load_spatial_data.py \
  --dataset_id spatial_10x_visium/mouse_brain_coronal_section1 \
  --output_dir ./output \
  --data_source_dir ../datasets
```

### Arguments

- `--dataset_id`: Dataset identifier matching the directory structure in `datasets/`
  - Example: `spatial_10x_visium/mouse_brain_coronal_section1`
- `--output_dir`: Directory where the dataset files will be copied
- `--data_source_dir`: (Optional) Base directory containing the datasets. Defaults to `../datasets` relative to the script location
- `--name`: (Optional) Module name for omnibenchmark compatibility
