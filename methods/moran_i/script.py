#!/usr/bin/env python3
"""
Moran's I method for spatially variable gene detection.

The MoranI global spatial auto-correlation statistics evaluates whether features (i.e. genes)
shows a pattern that is clustered, dispersed or random in the tissue are under consideration.

Omnibenchmark standard arguments:
  --output_dir: Directory where outputs will be saved
  --name: Dataset name
  --data.dataset: Input dataset h5ad file

Method-specific arguments:
  --coord_type_moran_i: Type of coordinate system ('grid' or 'generic', default: 'generic')
"""

import warnings
warnings.filterwarnings('ignore')

import anndata as ad
import squidpy as sq
import argparse
import os

# Parse command line arguments
parser = argparse.ArgumentParser(description="Moran's I method for spatially variable gene detection")
parser.add_argument('--output_dir', type=str, required=True,
                    help='Output directory where files will be saved')
parser.add_argument('--name', type=str, required=True,
                    help='Name of the dataset')
parser.add_argument('--data.dataset', type=str, required=True,
                    help='Path to input dataset h5ad file')
parser.add_argument('--coord_type_moran_i', type=str, default='generic',
                    choices=['grid', 'generic'],
                    help='Type of coordinate system')

args = parser.parse_args()

# Get input using getattr for dotted arguments
input_data_path = getattr(args, 'data.dataset')

# Construct output path
os.makedirs(args.output_dir, exist_ok=True)
output_path = os.path.join(args.output_dir, f"{args.name}.predictions.h5ad")

print(f'Reading input data from: {input_data_path}', flush=True)
adata = ad.read_h5ad(input_data_path)

print("Running Moran's I", flush=True)
sq.gr.spatial_neighbors(adata,
                        coord_type=args.coord_type_moran_i,
                        delaunay=True)

sq.gr.spatial_autocorr(adata,
                       mode="moran",
                       layer='normalized',
                       n_perms=100,
                       genes=adata.var_names)

# Format output
df = adata.uns["moranI"]
df = df.loc[adata.var_names][['I']]
df = df.reset_index()
df.columns = ['feature_id', 'pred_spatial_var_score']

output = ad.AnnData(
    var=df,
    uns={
        'dataset_id': adata.uns['dataset_id'],
        'method_id': 'moran_i'
    }
)

print(f"Writing output to: {output_path}", flush=True)
output.write_h5ad(output_path, compression='gzip')
print("Done!", flush=True)
