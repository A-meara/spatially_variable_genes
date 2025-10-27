import anndata as ad
import argparse

# Parse command line arguments
parser = argparse.ArgumentParser(description='True ranking control method for spatially variable genes')
parser.add_argument('--input_data', type=str, required=True,
                    help='Path to input dataset h5ad file')
parser.add_argument('--input_solution', type=str, required=True,
                    help='Path to solution h5ad file')
parser.add_argument('--output', type=str, required=True,
                    help='Path to output h5ad file')
parser.add_argument('--name', type=str, default='true_ranking',
                    help='Name of the method')

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
input_solution = ad.read_h5ad(par['input_solution'])

df = input_solution.var[["feature_id", "true_spatial_var_score"]]
df.rename(columns={'true_spatial_var_score': 'pred_spatial_var_score'}, inplace=True)

output = ad.AnnData(var=df,
                    uns={'dataset_id': input_solution.uns['dataset_id'],
                         'method_id': meta['name']})

print("Write output to file", flush=True)
output.write_h5ad(par['output'])
