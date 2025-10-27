import anndata as ad
import numpy as np
import argparse

# Parse command line arguments
parser = argparse.ArgumentParser(description='Random ranking control method for spatially variable genes')
parser.add_argument('--input_data', type=str, required=True,
                    help='Path to input dataset h5ad file')
parser.add_argument('--input_solution', type=str, required=True,
                    help='Path to solution h5ad file')
parser.add_argument('--output', type=str, required=True,
                    help='Path to output h5ad file')
parser.add_argument('--name', type=str, default='random_ranking',
                    help='Name of the method')
parser.add_argument('--seed', type=int, default=0,
                    help='Random seed for reproducibility')

args = parser.parse_args()

# Store parameters in dictionary for compatibility with existing code
par = {
    'input_data': args.input_data,
    'input_solution': args.input_solution,
    'output': args.output
}
meta = {
    'name': args.name
}

print('Generate predictions', flush=True)
input_data = ad.read_h5ad(par['input_data'])

df = input_data.var[["feature_id"]]

np.random.seed(args.seed)
df['pred_spatial_var_score'] = np.random.rand(len(df['feature_id']))

output = ad.AnnData(var=df,
                    uns={'dataset_id': input_data.uns['dataset_id'],
                         'method_id': meta['name']})

print("Write output to file", flush=True)
output.write_h5ad(par['output'])
