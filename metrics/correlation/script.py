#!/usr/bin/env python3
"""
Correlation metric for spatially variable genes benchmark.

Computes Kendall rank correlation coefficient between predicted and true
spatial variability scores.

Omnibenchmark standard arguments:
  --output_dir: Directory where outputs will be saved
  --name: Dataset name
  --methods.predictions: Method predictions h5ad file
  --data.solution: Ground truth solution h5ad file
"""

import anndata as ad
import pandas as pd
import argparse
import os

# Parse command line arguments
parser = argparse.ArgumentParser(description='Correlation metric for spatially variable genes')
parser.add_argument('--output_dir', type=str, required=True,
                    help='Output directory where files will be saved')
parser.add_argument('--name', type=str, required=True,
                    help='Name of the dataset')
parser.add_argument('--methods.predictions', type=str, required=True,
                    help='Path to method predictions h5ad file')
parser.add_argument('--data.solution', type=str, required=True,
                    help='Path to ground truth solution h5ad file')

args = parser.parse_args()

# Get inputs using getattr for dotted arguments
input_method_path = getattr(args, 'methods.predictions')
input_solution_path = getattr(args, 'data.solution')

# Construct output path
os.makedirs(args.output_dir, exist_ok=True)
output_path = os.path.join(args.output_dir, f"{args.name}.scores.h5ad")

print(f'Reading method predictions from: {input_method_path}', flush=True)
input_method = ad.read_h5ad(input_method_path)

print(f'Reading ground truth solution from: {input_solution_path}', flush=True)
input_solution = ad.read_h5ad(input_solution_path)

print('Computing correlation metric', flush=True)
df = pd.merge(input_method.var, input_solution.var, how='left', on='feature_id')
groupby = df.groupby('orig_feature_name', observed=True)
corr = groupby.apply(lambda x: x['pred_spatial_var_score'].corr(x['true_spatial_var_score'], method='kendall'))

uns_metric_ids = ['correlation']
uns_metric_values = [corr.mean()]

print("Writing output AnnData to file", flush=True)
output = ad.AnnData(
    uns={
        'dataset_id': input_method.uns['dataset_id'],
        'method_id': input_method.uns['method_id'],
        'metric_ids': uns_metric_ids,
        'metric_values': uns_metric_values
    }
)

print(f"Writing output to: {output_path}", flush=True)
output.write_h5ad(output_path, compression='gzip')
print("Done!", flush=True)
