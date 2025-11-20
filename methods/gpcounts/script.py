#!/usr/bin/env python3
"""
GPcounts method for spatially variable gene detection.

GPcounts is non-parametric modelling of temporal and spatial counts data
from RNA-seq experiments using Gaussian Processes with negative binomial likelihood.

Omnibenchmark standard arguments:
  --output_dir: Directory where outputs will be saved
  --name: Dataset name
  --data.dataset: Input dataset h5ad file

Method-specific arguments:
  --n_features: Number of features to include (default: None = all features)
"""

import statsmodels.api as sm
import statsmodels.formula.api as smf
import pandas as pd
import anndata as ad
import scipy
from GPcounts.RNA_seq_GP import rna_seq_gp
import warnings
import argparse
import os

warnings.filterwarnings('ignore')

# Parse command line arguments
parser = argparse.ArgumentParser(description='GPcounts method for spatially variable gene detection')
parser.add_argument('--output_dir', type=str, required=True,
                    help='Output directory where files will be saved')
parser.add_argument('--name', type=str, required=True,
                    help='Name of the dataset')
parser.add_argument('--data.dataset', type=str, required=True,
                    help='Path to input dataset h5ad file')
parser.add_argument('--n_features', type=int, default=None,
                    help='Number of features to include')

args = parser.parse_args()

# Get input using getattr for dotted arguments
input_data_path = getattr(args, 'data.dataset')

# Construct output path
os.makedirs(args.output_dir, exist_ok=True)
output_path = os.path.join(args.output_dir, f"{args.name}.predictions.h5ad")

print(f'Reading input data from: {input_data_path}', flush=True)
adata = ad.read_h5ad(input_data_path)

print('Running GPcounts', flush=True)

# Subset if required
if args.n_features:
    adata = adata[:, :args.n_features]

counts = adata.layers["counts"]
if scipy.sparse.issparse(counts):
    counts = counts.todense()

Y = pd.DataFrame(data=counts,
                index=adata.obs_names,
                columns=adata.var_names)

spatial_locations = pd.DataFrame(data=adata.obsm['spatial'],
                                index=adata.obs_names,
                                columns=['x', 'y'])
spatial_locations['total_counts'] = Y.sum(1)

Y = Y.loc[spatial_locations.index]
X = spatial_locations[['x', 'y']]

# Calculate scales
scales = []
for i in range(0, len(Y.columns)):
    model = smf.glm(formula="Y.iloc[:,i]~0+spatial_locations['total_counts']", data=Y,
                    family=sm.families.NegativeBinomial(sm.families.links.identity())).fit()
    res = model.params[0]*spatial_locations['total_counts']
    scales.append(res)
scalesdf = pd.DataFrame(scales)
scalesdf = scalesdf.T

Y = Y.T
X = X[['x', 'y']]

sparse = True
nb_scaled = True
gene_name = Y.index
likelihood = 'Negative_binomial'
gp_counts = rna_seq_gp(
    X, Y.loc[gene_name], sparse=sparse, M=250, scale=scalesdf, safe_mode=False)

log_likelihood_ratio = gp_counts.One_sample_test(likelihood)

df = gp_counts.calculate_FDR(log_likelihood_ratio)

# Format output
df = df.loc[adata.var_names][['log_likelihood_ratio']]
df = df.reset_index()
df.columns = ['feature_id', 'pred_spatial_var_score']

output = ad.AnnData(
    var=df,
    uns={
        'dataset_id': adata.uns['dataset_id'],
        'method_id': 'gpcounts'
    }
)

print(f"Writing output to: {output_path}", flush=True)
output.write_h5ad(output_path, compression='gzip')
print("Done!", flush=True)
