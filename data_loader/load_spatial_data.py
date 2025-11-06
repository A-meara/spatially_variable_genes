#!/usr/bin/env python3
"""
Data loader for spatially variable genes benchmark.

This module loads spatial transcriptomics datasets and outputs them in the format
expected by the omnibenchmark framework. It can load data from local h5ad files
or download them from a remote source.

Omnibenchmark standard arguments:
  --output_dir: Directory where outputs will be saved
  --name: Dataset name

Component-specific arguments:
  --dataset_id: Dataset identifier (e.g., "spatial_10x_visium/mouse_brain_coronal_section1")
"""

import argparse
import os
import sys
import shutil
from pathlib import Path


def load_dataset(dataset_id, output_dir, name=None, data_source_dir=None):
    """
    Load a spatial transcriptomics dataset and copy/symlink files to output directory.

    Args:
        dataset_id: Dataset identifier (e.g., "spatial_10x_visium/mouse_brain_coronal_section1")
        output_dir: Output directory where dataset files will be saved
        name: Dataset name for output files (default: uses dataset_id basename)
        data_source_dir: Base directory containing the datasets (default: ../datasets)

    Returns:
        Paths to the output files
    """
    # Set default data source directory relative to this script
    if data_source_dir is None:
        script_dir = Path(__file__).parent.resolve()
        data_source_dir = script_dir.parent / "datasets"
    else:
        data_source_dir = Path(data_source_dir)

    # Construct path to dataset
    dataset_path = data_source_dir / dataset_id

    if not dataset_path.exists():
        raise FileNotFoundError(f"Dataset not found: {dataset_path}")

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    # Use name parameter if provided, otherwise use dataset_id basename
    if name is None:
        name = Path(dataset_id).name

    # Expected files in the dataset directory
    expected_files = {
        'dataset': ('dataset.h5ad', f'{name}.dataset.h5ad'),
        'solution': ('solution.h5ad', f'{name}.solution.h5ad'),
        'simulated_dataset': ('simulated_dataset.h5ad', f'{name}.simulated_dataset.h5ad')
    }

    output_files = {}

    # Copy or symlink the files to output directory
    for file_type, (src_filename, dst_filename) in expected_files.items():
        src_file = dataset_path / src_filename
        dst_file = Path(output_dir) / dst_filename

        if src_file.exists():
            # Copy the file to output directory
            print(f"Copying {file_type}: {src_file} -> {dst_file}")
            shutil.copy2(src_file, dst_file)
            output_files[file_type] = str(dst_file)
        else:
            print(f"Warning: {file_type} file not found: {src_file}", file=sys.stderr)

    # Verify at minimum we have dataset file
    if 'dataset' not in output_files:
        raise FileNotFoundError(f"Required dataset.h5ad file not found in {dataset_path}")

    return output_files


def main():
    parser = argparse.ArgumentParser(
        description='Load spatial transcriptomics dataset for omnibenchmark'
    )

    parser.add_argument(
        '--dataset_id',
        type=str,
        required=True,
        help='Dataset identifier (e.g., "spatial_10x_visium/mouse_brain_coronal_section1")'
    )

    parser.add_argument(
        '--dataset_name',
        type=str,
        default=None,
        help='Dataset name for output files (e.g., "mouse_brain_coronal_section1"). If not provided, extracted from dataset_id'
    )

    parser.add_argument(
        '--output_dir',
        type=str,
        required=True,
        help='Output directory where dataset files will be saved'
    )

    parser.add_argument(
        '--data_source_dir',
        type=str,
        default=None,
        help='Base directory containing the datasets (default: ../datasets relative to script)'
    )

    parser.add_argument(
        '--name',
        type=str,
        default=None,
        help='Optional name for this data module (for omnibenchmark compatibility, alias for dataset_name)'
    )

    args = parser.parse_args()

    # Determine dataset name: prioritize --dataset_name, then --name, then extract from dataset_id
    dataset_name = args.dataset_name or args.name
    if dataset_name is None:
        # Extract from dataset_id (last component of path)
        dataset_name = Path(args.dataset_id).name

    try:
        output_files = load_dataset(
            dataset_id=args.dataset_id,
            output_dir=args.output_dir,
            name=dataset_name,
            data_source_dir=args.data_source_dir
        )

        print(f"\nSuccessfully loaded dataset: {args.dataset_id}")
        print(f"Output files:")
        for file_type, path in output_files.items():
            print(f"  {file_type}: {path}")

    except Exception as e:
        print(f"Error loading dataset: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
