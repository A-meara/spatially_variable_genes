#!/usr/bin/env Rscript
#
# BOOST-GP method for spatially variable gene detection.
#
# Omnibenchmark standard arguments:
#   --output_dir: Directory where outputs will be saved
#   --name: Dataset name
#   --data.dataset: Input dataset h5ad file
#
# Method-specific arguments:
#   --n_iter: Number of iterations (default: 10)

library(RcppDist)
library(anndata)
library(optparse)

# Define command-line arguments
option_list <- list(
  make_option(c("--output_dir"), type="character",
              help="Output directory path"),
  make_option(c("--name"), type="character",
              help="Dataset name"),
  make_option(c("--data_dataset"), type="character",
              help="Input dataset h5ad file path"),
  make_option(c("--n_iter"), type="integer", default=10,
              help="Number of iterations [default %default]")
)

# Parse arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Get input path
input_data_path <- opt$data_dataset

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

cat("Running BOOST-GP with", opt$n_iter, "iterations\n")
df <- as.data.frame(boost.gp(Y = counts, loc = loc, iter = opt$n_iter, burn = 5))

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
