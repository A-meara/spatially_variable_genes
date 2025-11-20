#!/usr/bin/env python3
"""
SOMDE method for spatially variable gene detection.

SOMDE is a scalable method for identifying spatially variable genes with
self-organizing map, offering 5 to 50 times faster performance than existing methods.

Omnibenchmark standard arguments:
  --output_dir: Directory where outputs will be saved
  --name: Dataset name
  --data.dataset: Input dataset h5ad file
"""

import anndata as ad
import pandas as pd
import numpy as np
import scanpy as sc
from somde import SomNode
import scipy
import argparse
import os

# Parse command line arguments
parser = argparse.ArgumentParser(description='SOMDE method for spatially variable gene detection')
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

print('Running SOMDE', flush=True)
counts = adata.layers["counts"]
if scipy.sparse.issparse(counts):
    counts = counts.todense()

data = pd.DataFrame(
    counts,
    columns=adata.var_names,
    index=adata.obs_names
)

X = pd.DataFrame(adata.obsm["spatial"],
                     index=adata.obs_names,
                     columns=["x", "y"]).values.astype(np.float32)

som = SomNode(X, k=10)
ndf, ninfo = som.mtx(data.transpose())
nres = som.norm()

df, SVnum = som.run()

# Format output
df.set_index("g", inplace=True)
df = df.loc[adata.var_names][['FSV']]
df = df.reset_index()
df.columns = ['feature_id', 'pred_spatial_var_score']

output = ad.AnnData(
    var=df,
    uns={
        'dataset_id': adata.uns['dataset_id'],
        'method_id': 'somde'
    }
)

print(f"Writing output to: {output_path}", flush=True)
output.write_h5ad(output_path, compression='gzip')
print("Done!", flush=True)
