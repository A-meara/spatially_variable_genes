#!/usr/bin/env python3
"""
SpaGFT method for spatially variable gene detection.

SpaGFT is a graph Fourier transform for tissue module identification from spatially resolved
transcriptomics, transforming complex gene expression patterns into simple but informative signals.

Omnibenchmark standard arguments:
  --output_dir: Directory where outputs will be saved
  --name: Dataset name
  --data.dataset: Input dataset h5ad file
"""

import anndata as ad
import SpaGFT as spg
import argparse
import os

# Parse command line arguments
parser = argparse.ArgumentParser(description='SpaGFT method for spatially variable gene detection')
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

print('Running SpaGFT', flush=True)

adata.X = adata.layers['normalized'].copy()

adata.obs.loc[:, ['array_row', 'array_col']] = adata.obsm['spatial']

(ratio_low, ratio_high) = spg.gft.determine_frequency_ratio(adata, ratio_neighbors=1)

df = spg.detect_svg(adata,
                    spatial_info=['array_row', 'array_col'],
                    ratio_low_freq=ratio_low,
                    ratio_high_freq=ratio_high,
                    ratio_neighbors=1,
                    filter_peaks=True,
                    S=6)

# Format output
df = df.loc[adata.var_names][['gft_score']]
df = df.reset_index()
df.columns = ['feature_id', 'pred_spatial_var_score']

output = ad.AnnData(
    var=df,
    uns={
        'dataset_id': adata.uns['dataset_id'],
        'method_id': 'spagft'
    }
)

print(f"Writing output to: {output_path}", flush=True)
output.write_h5ad(output_path, compression='gzip')
print("Done!", flush=True)
