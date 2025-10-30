#!/usr/bin/env python3
"""
Random ranking control method for spatially variable genes.

Assigns random scores to genes as a baseline control method.

Omnibenchmark standard arguments:
  --output_dir: Directory where outputs will be saved
  --name: Dataset name
  --data.dataset: Input dataset h5ad file
  --data.solution: Input solution h5ad file (for control methods)

Method-specific arguments:
  --seed: Random seed for reproducibility
"""

import anndata as ad
import numpy as np
import argparse
import os

# Parse command line arguments
parser = argparse.ArgumentParser(description='Random ranking control method for spatially variable genes')
parser.add_argument('--output_dir', type=str, required=True,
                    help='Output directory where files will be saved')
parser.add_argument('--name', type=str, required=True,
                    help='Name of the dataset')
parser.add_argument('--data.dataset', type=str, required=True,
                    help='Path to input dataset h5ad file')
parser.add_argument('--data.solution', type=str, required=True,
                    help='Path to solution h5ad file')
parser.add_argument('--seed', type=int, default=0,
                    help='Random seed for reproducibility')

args = parser.parse_args()

# Get inputs using getattr for dotted arguments
input_data_path = getattr(args, 'data.dataset')
input_solution_path = getattr(args, 'data.solution')

# Construct output path
os.makedirs(args.output_dir, exist_ok=True)
output_path = os.path.join(args.output_dir, f"{args.name}.h5ad")

print(f'Reading input data from: {input_data_path}', flush=True)
input_data = ad.read_h5ad(input_data_path)

print('Generating random predictions', flush=True)
df = input_data.var[["feature_id"]].copy()

np.random.seed(args.seed)
df['pred_spatial_var_score'] = np.random.rand(len(df['feature_id']))

output = ad.AnnData(
    var=df,
    uns={
        'dataset_id': input_data.uns['dataset_id'],
        'method_id': 'random_ranking'
    }
)

print(f"Writing output to: {output_path}", flush=True)
output.write_h5ad(output_path, compression='gzip')
print("Done!", flush=True)
