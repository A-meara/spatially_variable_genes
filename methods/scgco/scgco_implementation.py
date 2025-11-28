#!/usr/bin/env python3
"""
scGCO implementation - runs inside Python 3.9 container.

This is the actual scGCO method implementation that requires Python 3.9.
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
parser.add_argument('--input_data', type=str, required=True,
                    help='Path to input dataset h5ad file')
parser.add_argument('--output', type=str, required=True,
                    help='Path to output predictions h5ad file')
parser.add_argument('--name', type=str, required=True,
                    help='Name of the dataset')

args = parser.parse_args()

print(f'Reading input data from: {args.input_data}', flush=True)
adata = ad.read_h5ad(args.input_data)

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

print(f"Writing output to: {args.output}", flush=True)
output.write_h5ad(args.output, compression='gzip')
print("Done!", flush=True)
