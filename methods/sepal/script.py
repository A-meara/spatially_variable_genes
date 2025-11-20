#!/usr/bin/env python3
"""
Sepal method for spatially variable gene detection.

Sepal simulates diffusion of individual transcripts to extract genes with
spatial patterns and assesses the degree of randomness exhibited by each transcript profile.

Omnibenchmark standard arguments:
  --output_dir: Directory where outputs will be saved
  --name: Dataset name
  --data.dataset: Input dataset h5ad file

Method-specific arguments:
  --coord_type_sepal: Type of coordinate system ('grid' or 'generic', default: 'grid')
  --max_neighs_sepal: Maximum number of neighbors (4 or 6, default: 6)
"""

import anndata as ad
import squidpy as sq
import argparse
import os

# Parse command line arguments
parser = argparse.ArgumentParser(description='Sepal method for spatially variable gene detection')
parser.add_argument('--output_dir', type=str, required=True,
                    help='Output directory where files will be saved')
parser.add_argument('--name', type=str, required=True,
                    help='Name of the dataset')
parser.add_argument('--data.dataset', type=str, required=True,
                    help='Path to input dataset h5ad file')
parser.add_argument('--coord_type_sepal', type=str, default='grid',
                    choices=['grid', 'generic'],
                    help='Type of coordinate system')
parser.add_argument('--max_neighs_sepal', type=int, default=6,
                    choices=[4, 6],
                    help='Maximum number of neighbors')

args = parser.parse_args()

# Get input using getattr for dotted arguments
input_data_path = getattr(args, 'data.dataset')

# Construct output path
os.makedirs(args.output_dir, exist_ok=True)
output_path = os.path.join(args.output_dir, f"{args.name}.predictions.h5ad")

print(f'Reading input data from: {input_data_path}', flush=True)
adata = ad.read_h5ad(input_data_path)

print('Running Sepal', flush=True)
sq.gr.spatial_neighbors(adata,
                        coord_type=args.coord_type_sepal,
                        delaunay=False)

sq.gr.sepal(adata,
            layer='normalized',
            max_neighs=args.max_neighs_sepal,
            genes=adata.var_names,
            n_jobs=1)

# Format output
df = adata.uns["sepal_score"]
df = df.loc[adata.var_names][['sepal_score']]
df = df.reset_index()
df.columns = ['feature_id', 'pred_spatial_var_score']

output = ad.AnnData(
    var=df,
    uns={
        'dataset_id': adata.uns['dataset_id'],
        'method_id': 'sepal'
    }
)

print(f"Writing output to: {output_path}", flush=True)
output.write_h5ad(output_path, compression='gzip')
print("Done!", flush=True)
