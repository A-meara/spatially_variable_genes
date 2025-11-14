#!/usr/bin/env Rscript
"SPARK-X method for spatially variable gene detection.

SPARK-X is a non-parametric method for rapid and effective detection of
spatially expressed genes in large spatial transcriptomic studies.

Usage:
  script.R --output_dir=<dir> --name=<name> --data.dataset=<file> [--num_cores=<n>]
  script.R (-h | --help)

Options:
  --output_dir=<dir>     Output directory where files will be saved
  --name=<name>          Dataset name
  --data.dataset=<file>  Input dataset h5ad file path
  --num_cores=<n>        Number of CPU cores to use [default: 4]
  -h --help              Show this help message
" -> doc

suppressMessages(library(SPARK))
suppressMessages(library(anndata))
suppressMessages(library(docopt))

# Parse arguments
opt <- docopt(doc)

# Get input path
input_data_path <- opt$data.dataset
num_cores <- as.integer(opt$num_cores)

# Create output directory and construct output path
dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)
output_file <- file.path(opt$output_dir, paste0(opt$name, ".predictions.h5ad"))

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
cat("Running SPARK-X with", num_cores, "cores\n")
sparkX <- sparkx(counts, info[, 1:2], numCores = num_cores, option = "mixture")

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
