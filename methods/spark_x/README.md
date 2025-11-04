# SPARK-X Method for Omnibenchmark

## Overview

SPARK-X is a non-parametric method for rapid and effective detection of spatially expressed genes in large spatial transcriptomic studies.

## What Was Done

✅ Created omnibenchmark-compatible `script.R`
✅ Created `config.cfg` for omnibenchmark
✅ Updated `spatially_variable_genes_apptainer.yml`

## Files

- **`script.R`**: Omnibenchmark-compatible R script
  - Uses standard arguments: `--output_dir`, `--name`, `--data.dataset`
  - Configurable CPU cores via `--num_cores` (default: 4)

- **`config.cfg`**: Points to script.R

- **`config.vsh.yaml`**: Original Viash config (kept for reference)

## Usage via Dispatcher

```bash
python run_module.py \
  --component methods/spark_x \
  --output_dir ./test_output \
  --name test_dataset \
  --data.dataset ./path/to/dataset.h5ad \
  --num_cores 4
```

## Apptainer Image - ACTION REQUIRED!

⚠️ **You need to build the apptainer image: `envs/spark_x.sif`**

### Option 1: Build from Docker image

The original uses `openproblems/base_r:1.0.0` with SPARK installed:

```bash
# Create a definition file
cat > spark_x.def <<'EOF'
Bootstrap: docker
From: openproblems/base_r:1.0.0

%post
    R -e 'remotes::install_github("xzhoulab/SPARK")'
EOF

# Build the image
apptainer build envs/spark_x.sif spark_x.def
```

### Option 2: Use existing Docker image

If you have a Docker image already:

```bash
apptainer build envs/spark_x.sif docker://your_docker_image:tag
```

### Option 3: Pull from a registry

If the image is in a registry:

```bash
apptainer pull envs/spark_x.sif oras://registry.url/image:tag
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

SPARK-X is now active in `spatially_variable_genes_apptainer.yml`:

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

## Next Steps

1. **Build the apptainer image** (see above)
2. **Test locally** with the test command
3. **Commit and push** changes
4. **Update commit hash** in yml
5. **Run omnibenchmark**:
   ```bash
   ob run validate -b spatially_variable_genes_apptainer.yml
   ob run benchmark --benchmark spatially_variable_genes_apptainer.yml --cores 1 --local-storage --dry
   ```

## References

- Paper: https://doi.org/10.1186/s13059-021-02404-0
- Repository: https://github.com/xzhoulab/SPARK
- Documentation: https://xzhoulab.github.io/SPARK/02_SPARK_Example/
