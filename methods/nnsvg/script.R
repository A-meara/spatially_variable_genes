#!/usr/bin/env Rscript
#
# nnSVG method for spatially variable gene detection.
#
# nnSVG identifies genes that vary in expression continuously across the entire tissue
# or within a priori defined spatial domains using nearest-neighbor Gaussian processes.

suppressMessages(library(SpatialExperiment))
suppressMessages(library(scran))
suppressMessages(library(nnSVG))
suppressMessages(library(anndata))
suppressMessages(library(dplyr))
suppressMessages(library(argparse))

# Define command-line arguments
parser <- ArgumentParser(description='nnSVG method for spatially variable gene detection')
parser$add_argument('--output_dir', type='character', required=TRUE,
                    help='Output directory where files will be saved')
parser$add_argument('--name', type='character', required=TRUE,
                    help='Dataset name')
parser$add_argument('--data.dataset', type='character', required=TRUE,
                    dest='data_dataset',
                    help='Input dataset h5ad file path')
parser$add_argument('--n_threads', type='integer', default=1,
                    help='Number of threads to use')

# Parse arguments
opt <- parser$parse_args()

# Get input path
input_data_path <- opt$data_dataset

# Create output directory and construct output path
dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)
output_file <- file.path(opt$output_dir, paste0(opt$name, ".predictions.h5ad"))

# Load data
cat("Reading data from:", input_data_path, "\n")
adata <- read_h5ad(input_data_path)
counts <- t(as.matrix(adata$layers[['counts']]))

colnames(counts) <- adata$obs_names
rownames(counts) <- adata$var_names

loc <- as.data.frame(adata$obsm[['spatial']])

row_data = adata$var
row_data$gene_id = rownames(row_data)
row_data$feature_type = "Gene Expression"

colnames(loc) <- c("x", "y")
rownames(loc) <- colnames(counts)

spe <- SpatialExperiment(
    assays = list(counts = counts),
    rowData = row_data,
    colData = loc,
    spatialCoordsNames = c("x", "y"))

# Calculate logcounts (log-transformed normalized counts) using scran package
# using library size factors
spe <- computeLibraryFactors(spe)
spe <- logNormCounts(spe)

# Run nnSVG
cat("Running nnSVG with", opt$n_threads, "threads\n")
spe <- nnSVG(spe, n_threads=opt$n_threads)

# Format output
df <- as.data.frame(rowData(spe)) %>%
    subset(select = c('feature_id', 'LR_stat'))

colnames(df) <- c('feature_id', 'pred_spatial_var_score')
rownames(df) <- NULL

# Save output
cat("Writing output to:", output_file, "\n")
output = anndata::AnnData(
    shape = adata$shape,
    var=df,
    uns=list('dataset_id' = adata$uns[['dataset_id']],
             'method_id' = 'nnsvg'))

anndata::write_h5ad(anndata = output, filename = output_file)
cat("Done!\n")
