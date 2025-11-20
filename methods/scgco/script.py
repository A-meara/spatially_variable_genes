#!/usr/bin/env python3
"""
scGCO method for spatially variable gene detection.

scGCO uses fast optimization of hidden Markov Random Fields with graph cuts
to identify spatially variable genes with superior performance and scalability.

Omnibenchmark standard arguments:
  --output_dir: Directory where outputs will be saved
  --name: Dataset name
  --data.dataset: Input dataset h5ad file
"""

import warnings
warnings.filterwarnings('ignore')

import pandas as pd
import anndata as ad
import numpy as np
import scipy
import sys
sys.path.append("/opt/scGCO")

from scGCO_simple import *
import argparse
import os

# Parse command line arguments
parser = argparse.ArgumentParser(description='scGCO method for spatially variable gene detection')
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

counts = adata.layers["counts"]
if scipy.sparse.issparse(counts):
    counts = counts.todense()

data = pd.DataFrame(
    counts,
    columns=adata.var_names,
    index=adata.obs_names
)

print('Running scGCO', flush=True)
data_norm = normalize_count_cellranger(data)

exp = data.iloc[:, 0]
locs = adata.obsm['spatial'].copy()

print('Creating graph with weight', flush=True)
cellGraph = create_graph_with_weight(locs, exp)
gmmDict = gmm_model(data_norm)

print('Identifying spatial genes', flush=True)
df = identify_spatial_genes(locs, data_norm, cellGraph, gmmDict)

# Format output
df = df.loc[adata.var_names][['fdr']]
df = df.reset_index()
df.columns = ['feature_id', 'pred_spatial_var_score']

# Transform the values via -log10 to make sure a bigger score represents a
# higher spatial variation
df['pred_spatial_var_score'] = -np.log10(df['pred_spatial_var_score'].tolist())

output = ad.AnnData(
    var=df,
    uns={
        'dataset_id': adata.uns['dataset_id'],
        'method_id': 'scgco'
    }
)

print(f"Writing output to: {output_path}", flush=True)
output.write_h5ad(output_path, compression='gzip')
print("Done!", flush=True)
