import argparse
import anndata as ad
import numpy as np

def main():
    parser = argparse.ArgumentParser(description='Random ranking method for spatially variable genes')
    parser.add_argument('--input_data', type=str, required=True,
                        help='Path to input dataset (h5ad file)')
    parser.add_argument('--input_solution', type=str, required=False,
                        help='Path to solution dataset (h5ad file) - not used by this method')
    parser.add_argument('--output', type=str, required=True,
                        help='Path to output file (h5ad)')
    parser.add_argument('--method_id', type=str, default='random_ranking',
                        help='Method identifier')
    parser.add_argument('--seed', type=int, default=0,
                        help='Random seed for reproducibility')

    args = parser.parse_args()

    print('Generate predictions', flush=True)
    input_data = ad.read_h5ad(args.input_data)

    df = input_data.var[["feature_id"]].copy()

    np.random.seed(args.seed)
    df['pred_spatial_var_score'] = np.random.rand(len(df['feature_id']))

    output = ad.AnnData(var=df,
                        uns={'dataset_id': input_data.uns['dataset_id'],
                             'method_id': args.method_id})

    print("Write output to file", flush=True)
    output.write_h5ad(args.output)
    print("Done!", flush=True)

if __name__ == '__main__':
    main()
    