library(RcppDist)
library(anndata)
library(optparse)

# Define command-line arguments
option_list <- list(
  make_option(c("--input_data"), type="character",
              help="Input h5ad file path"),
  make_option(c("--output_dir"), type="character",
              help="Output directory path"),
  make_option(c("--name"), type="character",
              help="Dataset name"),
  make_option(c("--n_iter"), type="integer", default=10,
              help="Number of iterations [default %default]")
)

# Parse arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Set output file path
output_file <- file.path(opt$output_dir, paste0(opt$name, ".h5ad"))

cat("Load data\n")
adata <- anndata::read_h5ad(opt$input_data)

# Change to BOOST-GP directory (if needed)
setwd("/opt/BOOST-GP")
source("./R/boost.gp.R")

counts <- as.matrix(adata$layers[["counts"]])
colnames(counts) <- adata$var_names
rownames(counts) <- adata$obs_names
mode(counts) <- "integer"

loc <- as.data.frame(adata$obsm[["spatial"]])
rownames(loc) <- adata$obs_names
colnames(loc) <- c("x", "y")

cat("Run BOOST-GP\n")
df <- as.data.frame(boost.gp(Y = counts, loc = loc, iter = opt$n_iter, burn = 5))

df$feature_id <- rownames(df)
df <- subset(df, select = c("feature_id", "PPI"))
colnames(df) <- c("feature_id", "pred_spatial_var_score")

# Save output
cat("Write output AnnData to file\n")
output <- anndata::AnnData(
    shape = adata$shape,
    var = df,
    uns = list(
        "dataset_id" = adata$uns[["dataset_id"]],
        "method_id" = "BOOST-GP"
    )
)

zzz <- output$write_h5ad(output_file, compression = "gzip")
cat("Done!\n")x