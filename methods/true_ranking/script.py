#!/usr/bin/env python3
"""
True ranking control method for spatially variable genes.

Uses ground truth scores as predictions - represents perfect performance upper bound.

Omnibenchmark standard arguments:
  --output_dir: Directory where outputs will be saved
  --name: Dataset name
  --data.dataset: Input dataset h5ad file
  --data.solution: Input solution h5ad file (for control methods)
"""

import anndata as ad
import argparse
import os

# Parse command line arguments
parser = argparse.ArgumentParser(description='True ranking control method for spatially variable genes')
parser.add_argument('--output_dir', type=str, required=True,
                    help='Output directory where files will be saved')
parser.add_argument('--name', type=str, required=True,
                    help='Name of the dataset')
parser.add_argument('--data.dataset', type=str, required=True,
                    help='Path to input dataset h5ad file')
parser.add_argument('--data.solution', type=str, required=True,
                    help='Path to solution h5ad file')

args = parser.parse_args()

# Get inputs using getattr for dotted arguments
input_data_path = getattr(args, 'data.dataset')
input_solution_path = getattr(args, 'data.solution')

# Construct output path
os.makedirs(args.output_dir, exist_ok=True)
output_path = os.path.join(args.output_dir, f"{args.name}.predictions.h5ad")

print(f'Reading solution from: {input_solution_path}', flush=True)
input_solution = ad.read_h5ad(input_solution_path)

print('Copying true scores as predictions', flush=True)
df = input_solution.var[["feature_id", "true_spatial_var_score"]].copy()
df.rename(columns={'true_spatial_var_score': 'pred_spatial_var_score'}, inplace=True)

output = ad.AnnData(
    var=df,
    uns={
        'dataset_id': input_solution.uns['dataset_id'],
        'method_id': 'true_ranking'
    }
)

print(f"Writing output to: {output_path}", flush=True)
output.write_h5ad(output_path, compression='gzip')
print("Done!", flush=True)
