#!/usr/bin/env python3
"""
Spanve method for spatially variable gene detection.

Spanve is a non-parametric statistical approach based on modeling space dependence as a
distance of two distributions for detecting SV genes with high computing efficiency and accuracy.

Omnibenchmark standard arguments:
  --output_dir: Directory where outputs will be saved
  --name: Dataset name
  --data.dataset: Input dataset h5ad file
"""

import anndata as ad
from Spanve import Spanve
import argparse
import os

# Parse command line arguments
parser = argparse.ArgumentParser(description='Spanve method for spatially variable gene detection')
parser.add_argument('--output_dir', type=str, required=True,
                    help='Output directory where files will be saved')
parser.add_argument('--name', type=str, required=True,
                    help='Name of the dataset')
parser.add_argument('--data.dataset', type=str, required=True,
                    help='Path to input dataset h5ad file')

args = parser.parse_args()

# Get input using getattr for dotted arguments
input_data_path = getattr(args, 'data.dataset')

# Construct output path
os.makedirs(args.output_dir, exist_ok=True)
output_path = os.path.join(args.output_dir, f"{args.name}.predictions.h5ad")

print(f'Reading input data from: {input_data_path}', flush=True)
adata = ad.read_h5ad(input_data_path)

print('Running Spanve', flush=True)
adata.X = adata.layers['counts']
spanve = Spanve(adata)
spanve.fit(verbose=False)

# Format output
df = spanve.result_df
df = df.loc[adata.var_names][['ent']]
df = df.reset_index()
df.columns = ['feature_id', 'pred_spatial_var_score']

output = ad.AnnData(
    var=df,
    uns={
        'dataset_id': adata.uns['dataset_id'],
        'method_id': 'spanve'
    }
)

print(f"Writing output to: {output_path}", flush=True)
output.write_h5ad(output_path, compression='gzip')
print("Done!", flush=True)
