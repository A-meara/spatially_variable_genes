#!/usr/bin/env python3
"""
SpatialDE method for spatially variable gene detection.

SpatialDE decomposes expression variability into spatial and nonspatial components using
two random effect terms based on a Gaussian Process model.

Omnibenchmark standard arguments:
  --output_dir: Directory where outputs will be saved
  --name: Dataset name
  --data.dataset: Input dataset h5ad file
"""

import warnings
warnings.filterwarnings('ignore')

import scanpy as sc
import anndata as ad
import NaiveDE
import SpatialDE
import argparse
import os

# Parse command line arguments
parser = argparse.ArgumentParser(description='SpatialDE method for spatially variable gene detection')
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

# Run spatialDE
print('Running SpatialDE', flush=True)
sc.pp.calculate_qc_metrics(adata,
                           layer='counts',
                           inplace=True,
                           percent_top=[10])

counts = sc.get.obs_df(adata,
                       keys=list(adata.var_names),
                       use_raw=False,
                       layer='counts')

total_counts = sc.get.obs_df(adata, keys=["total_counts"])
norm_expr = NaiveDE.stabilize(counts.T).T
resid_expr = NaiveDE.regress_out(total_counts,
                                 norm_expr.T,
                                 "np.log(total_counts)").T

df = SpatialDE.run(adata.obsm["spatial"], resid_expr)

# Format output
df.set_index("g", inplace=True)
df = df.loc[adata.var_names][['FSV']]
df = df.reset_index()
df.columns = ['feature_id', 'pred_spatial_var_score']

output = ad.AnnData(
    var=df,
    uns={
        'dataset_id': adata.uns['dataset_id'],
        'method_id': 'spatialde'
    }
)

print(f"Writing output to: {output_path}", flush=True)
output.write_h5ad(output_path, compression='gzip')
print("Done!", flush=True)
