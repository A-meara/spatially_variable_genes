#!/usr/bin/env Rscript
#
# SPARK method for spatially variable gene detection.
#
# SPARK builds upon a generalized linear spatial model (GLSM) with a variety of
# spatial kernels to accommodate count data using a penalized quasi-likelihood algorithm.

suppressMessages(library(SPARK))
suppressMessages(library(anndata))
suppressMessages(library(argparse))

# Define command-line arguments
parser <- ArgumentParser(description='SPARK method for spatially variable gene detection')
parser$add_argument('--output_dir', type='character', required=TRUE,
                    help='Output directory where files will be saved')
parser$add_argument('--name', type='character', required=TRUE,
                    help='Dataset name')
parser$add_argument('--data.dataset', type='character', required=TRUE,
                    dest='data_dataset',
                    help='Input dataset h5ad file path')
parser$add_argument('--num_cores', type='integer', default=1,
                    help='Number of CPU cores to use')

# Parse arguments
opt <- parser$parse_args()

# Get input path
input_data_path <- opt$data_dataset

# Create output directory and construct output path
dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)
output_file <- file.path(opt$output_dir, paste0(opt$name, ".predictions.h5ad"))

# Load data
cat("Reading data from:", input_data_path, "\n")
adata <- anndata::read_h5ad(input_data_path)
counts <- t(as.matrix(adata$layers[["counts"]]))
colnames(counts) <- adata$obs_names
rownames(counts) <- adata$var_names
info <- as.data.frame(adata$obsm[["spatial"]])
rownames(info) <- colnames(counts)
colnames(info) <- c("x", "y")

# Run SPARK
cat("Running SPARK with", opt$num_cores, "cores\n")
spark <- CreateSPARKObject(
    counts = counts, percentage = 0,
    min_total_counts = 0, location = info[, 1:2]
)

spark@lib_size <- apply(spark@counts, 2, sum)
spark <- spark.vc(spark,
    covariates = NULL,
    lib_size = spark@lib_size,
    num_core = opt$num_cores,
    verbose = FALSE
)

## Calculating pval
spark <- spark.test(spark,
    check_positive = T,
    verbose = F
)

df <- as.data.frame(spark@res_mtest)

df$feature_id <- rownames(df)

df <- subset(df, select = c("feature_id", "adjusted_pvalue"))
colnames(df) <- c("feature_id", "pred_spatial_var_score")

# Transform the values via -log10 to make sure a bigger score represents a higher spatial variation
df$pred_spatial_var_score <- -log10(df$pred_spatial_var_score)

# Save output
cat("Writing output to:", output_file, "\n")
output <- anndata::AnnData(
    shape = adata$shape,
    var = df,
    uns = list(
        "dataset_id" = adata$uns[["dataset_id"]],
        "method_id" = "spark"
    )
)

anndata::write_h5ad(anndata = output, filename = output_file)
cat("Done!\n")
