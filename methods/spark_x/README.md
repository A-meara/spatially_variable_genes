# SPARK-X Method for Omnibenchmark

## Overview

SPARK-X is a non-parametric method for rapid and effective detection of spatially expressed genes in large spatial transcriptomic studies.

## Files

- **`script.R`**: Omnibenchmark-compatible R script
  - Uses standard arguments: `--output_dir`, `--name`, `--data.dataset`
  - Configurable CPU cores via `--num_cores` (default: 4)

- **`config.cfg`**: Points to script.R

- **`config.vsh.yaml`**: Original Viash config (kept for reference)

## Usage via Dispatcher

```bash
python ob_run_component.py \
  --component methods/spark_x \
  --output_dir ./test_output \
  --name test_dataset \
  --data.dataset ./path/to/dataset.h5ad \
  --num_cores 4
```

## Dependencies

The R environment needs:
- **SPARK** package (from GitHub: xzhoulab/SPARK)
- **anndata** R package
- **optparse** package

## Method Details

**Input:**
- Spatial transcriptomics data with counts in `layers["counts"]`
- Spatial coordinates in `obsm["spatial"]`

**Output:**
- Adjusted p-values transformed via `-log10` as `pred_spatial_var_score`
- Higher scores = higher spatial variation

**Parameters:**
- `num_cores`: Number of CPU cores (default: 4)
- Method uses `option = "mixture"` for SPARK-X

## Testing Locally

Before running the full benchmark:

```bash
# Make sure you're in the repo root
cd /Users/aidan/Robinson_project/spatially_variable_genes

# Test the script
Rscript methods/spark_x/script.R \
  --output_dir ./test_spark_x \
  --name test_dataset \
  --data.dataset ./datasets/spatial_10x_visium/mouse_brain_coronal_section1/dataset.h5ad \
  --num_cores 2
```

## In the YML

SPARK-X is active in `spatially_variable_genes_apptainer.yml`:

```yaml
- id: spark_x
  name: "SPARK-X"
  software_environment: "spark_x"
  repository:
    url: https://github.com/A-meara/spatially_variable_genes
    commit: 40c5186
  parameters:
    - values: ["--component", "methods/spark_x", "--num_cores", "4"]
```
## References

- Paper: https://doi.org/10.1186/s13059-021-02404-0
- Repository: https://github.com/xzhoulab/SPARK
- Documentation: https://xzhoulab.github.io/SPARK/02_SPARK_Example/
