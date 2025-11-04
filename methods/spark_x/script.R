#!/usr/bin/env Rscript
#
# SPARK-X method for spatially variable gene detection.
#
# SPARK-X is a non-parametric method for rapid and effective detection of
# spatially expressed genes in large spatial transcriptomic studies.
#
# Omnibenchmark standard arguments:
#   --output_dir: Directory where outputs will be saved
#   --name: Dataset name
#   --data.dataset: Input dataset h5ad file
#
# Method-specific arguments:
#   --num_cores: Number of CPU cores to use (default: 4)

suppressMessages(library(SPARK))
suppressMessages(library(anndata))
suppressMessages(library(optparse))

# Define command-line arguments
option_list <- list(
  make_option(c("--output_dir"), type="character",
              help="Output directory path"),
  make_option(c("--name"), type="character",
              help="Dataset name"),
  make_option(c("--data.dataset"), type="character",
              help="Input dataset h5ad file path"),
  make_option(c("--num_cores"), type="integer", default=4,
              help="Number of CPU cores to use [default %default]")
)

# Parse arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Get input using list accessor for dotted argument
input_data_path <- opt[["data.dataset"]]

# Create output directory and construct output path
dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)
output_file <- file.path(opt$output_dir, paste0(opt$name, ".h5ad"))

# Load data
cat("Reading data from:", input_data_path, "\n")
adata <- anndata::read_h5ad(input_data_path)

# Extract counts and spatial coordinates
cat("Extracting counts and spatial coordinates\n")
counts <- t(as.matrix(adata$layers[["counts"]]))
colnames(counts) <- adata$obs_names
rownames(counts) <- adata$var_names

info <- as.data.frame(adata$obsm[["spatial"]])
rownames(info) <- colnames(counts)
colnames(info) <- c("x", "y")

# Run SPARK-X
cat("Running SPARK-X with", opt$num_cores, "cores\n")
sparkX <- sparkx(counts, info[, 1:2], numCores = opt$num_cores, option = "mixture")

# Format output
cat("Processing results\n")
df <- as.data.frame(sparkX$res_mtest)
df$feature_id <- rownames(df)
df <- subset(df, select = c("feature_id", "adjustedPval"))
colnames(df) <- c("feature_id", "pred_spatial_var_score")
rownames(df) <- NULL

# Transform p-values via -log10 so higher scores = higher spatial variation
df$pred_spatial_var_score <- -log10(df$pred_spatial_var_score)

# Save output
cat("Writing output to:", output_file, "\n")
output <- anndata::AnnData(
    shape = adata$shape,
    var = df,
    uns = list(
        "dataset_id" = adata$uns[["dataset_id"]],
        "method_id" = "spark_x"
    )
)

anndata::write_h5ad(anndata = output, filename = output_file)
cat("Done!\n")
