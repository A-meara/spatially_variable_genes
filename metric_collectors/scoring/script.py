#!/usr/bin/env python3
"""
Aggregate individual metric scores from multiple h5ad files into a unified table.
This module collects all metric outputs (h5ad files with scores in uns) and
aggregates them into a CSV file for analysis.
"""

import argparse
import glob
import os
import pandas as pd
import anndata as ad

def parse_args():
    parser = argparse.ArgumentParser(description='Aggregate metric scores from multiple h5ad files')
    parser.add_argument('--input_dir', type=str, required=True,
                       help='Directory containing metric h5ad files')
    parser.add_argument('--input_pattern', type=str, default='**/*_scores.h5ad',
                       help='Glob pattern to match metric files (default: **/*_scores.h5ad)')
    parser.add_argument('--output', type=str, required=True,
                       help='Output CSV file path for aggregated scores')
    return parser.parse_args()

def main():
    args = parse_args()

    print(f"Searching for metric files in: {args.input_dir}", flush=True)
    print(f"Using pattern: {args.input_pattern}", flush=True)

    # Find all metric h5ad files
    search_path = os.path.join(args.input_dir, args.input_pattern)
    metric_files = glob.glob(search_path, recursive=True)

    print(f"Found {len(metric_files)} metric files", flush=True)

    if len(metric_files) == 0:
        print("WARNING: No metric files found!", flush=True)
        # Create empty output file
        df = pd.DataFrame(columns=['dataset_id', 'method_id', 'metric_id', 'metric_value'])
        df.to_csv(args.output, index=False)
        return

    # Collect all scores
    records = []

    for metric_file in metric_files:
        print(f"Processing: {metric_file}", flush=True)
        try:
            # Read the h5ad file
            adata = ad.read_h5ad(metric_file)

            # Extract metadata from uns
            dataset_id = adata.uns.get('dataset_id', 'unknown')
            method_id = adata.uns.get('method_id', 'unknown')
            metric_ids = adata.uns.get('metric_ids', [])
            metric_values = adata.uns.get('metric_values', [])

            # Create a record for each metric
            for metric_id, metric_value in zip(metric_ids, metric_values):
                records.append({
                    'dataset_id': dataset_id,
                    'method_id': method_id,
                    'metric_id': metric_id,
                    'metric_value': metric_value
                })

        except Exception as e:
            print(f"ERROR processing {metric_file}: {e}", flush=True)
            continue

    print(f"Aggregated {len(records)} score records", flush=True)

    # Create DataFrame
    df = pd.DataFrame(records)

    # Sort for consistent output
    df = df.sort_values(['dataset_id', 'method_id', 'metric_id']).reset_index(drop=True)

    # Write to CSV
    print(f"Writing aggregated scores to: {args.output}", flush=True)
    df.to_csv(args.output, index=False)

    # Print summary statistics
    print("\n=== Summary ===", flush=True)
    print(f"Total scores: {len(df)}", flush=True)
    print(f"Unique datasets: {df['dataset_id'].nunique()}", flush=True)
    print(f"Unique methods: {df['method_id'].nunique()}", flush=True)
    print(f"Unique metrics: {df['metric_id'].nunique()}", flush=True)

    # Print preview
    print("\n=== Score Preview ===", flush=True)
    print(df.head(10).to_string(), flush=True)

    print("\nMetric aggregation complete!", flush=True)

if __name__ == '__main__':
    main()
