# Spatially Variable Genes Benchmark

OmniBenchmark implementation of the OpenProblems spatially variable genes task. Evaluates 16 SVG detection methods on 50 spatial transcriptomics datasets using Kendall rank correlation.

## Methods
BOOST-GP, GPcounts, Moran's I, nnSVG, scGCO, SEPAL, SomDE, SpaGCN, SpaGFT, SpaNVE, SPARK, SPARK-X, SpatialDE, SpatialDE2, random_ranking, true_ranking

## Data
50 preprocessed spatial transcriptomics datasets (H5AD format) available from OpenProblems. See [datasets documentation](https://openproblems.bio/benchmarks/spatially_variable_genes?version=v1.0.0) for download instructions.

## Usage
Run the full benchmark:
```bash
ob run benchmark -b spatially_variable_genes_apptainer.yml --local-storage
```

Preview execution without running:
```bash
ob run benchmark -b spatially_variable_genes_apptainer.yml --dry
```

Run in parallel with 8 cores:
```bash
ob run benchmark -b spatially_variable_genes_apptainer.yml --local-storage --cores 8
```

## Requirements
- OmniBenchmark (0.3+)
- Apptainer/Singularity (v1.0+)
- 70 GB disk space

## Related Repositories
- [Metric Collector](https://github.com/A-meara/metric_collector_SVG) - Aggregates scores from all methods and datasets

## References
- [OpenProblems SVG Task](https://github.com/openproblems-bio/task_spatially_variable_genes)
- [OmniBenchmark](https://omnibenchmark.org)
