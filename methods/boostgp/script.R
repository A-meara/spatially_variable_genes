#!/usr/bin/env Rscript
#
# BOOST-GP method for spatially variable gene detection.
#
# Bayesian modeling of spatial molecular profiling data via Gaussian process.

library(RcppDist)
library(anndata)
library(argparse)

# Define command-line arguments
parser <- ArgumentParser(description='BOOST-GP method for spatially variable gene detection')
parser$add_argument('--output_dir', type='character', required=TRUE,
                    help='Output directory where files will be saved')
parser$add_argument('--name', type='character', required=TRUE,
                    help='Dataset name')
parser$add_argument('--data.dataset', type='character', required=TRUE,
                    dest='data_dataset',
                    help='Input dataset h5ad file path')
parser$add_argument('--n_iter', type='integer', default=10,
                    help='Number of iterations')

# Parse arguments
opt <- parser$parse_args()

# Get input path
input_data_path <- opt$data_dataset
n_iter <- opt$n_iter

# Create output directory and construct output path
dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)
output_file <- file.path(opt$output_dir, paste0(opt$name, ".predictions.h5ad"))

cat("Reading data from:", input_data_path, "\n")
adata <- anndata::read_h5ad(input_data_path)

# Change to BOOST-GP directory (if needed)
setwd("/opt/BOOST-GP")
source("./R/boost.gp.R")

# Extract counts and spatial coordinates
counts <- as.matrix(adata$layers[["counts"]])
colnames(counts) <- adata$var_names
rownames(counts) <- adata$obs_names
mode(counts) <- "integer"

loc <- as.data.frame(adata$obsm[["spatial"]])
rownames(loc) <- adata$obs_names
colnames(loc) <- c("x", "y")

cat("Running BOOST-GP with", n_iter, "iterations\n")
df <- as.data.frame(boost.gp(Y = counts, loc = loc, iter = n_iter, burn = 5))

# Format output
df$feature_id <- rownames(df)
df <- subset(df, select = c("feature_id", "PPI"))
colnames(df) <- c("feature_id", "pred_spatial_var_score")

# Save output
cat("Writing output to:", output_file, "\n")
output <- anndata::AnnData(
    shape = adata$shape,
    var = df,
    uns = list(
        "dataset_id" = adata$uns[["dataset_id"]],
        "method_id" = "boostgp"
    )
)

zzz <- output$write_h5ad(output_file, compression = "gzip")
cat("Done!\n")
