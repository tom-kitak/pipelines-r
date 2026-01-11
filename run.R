#!/usr/bin/env Rscript

library(argparse)
library(jsonlite)
library(zellkonverter)

parser <- ArgumentParser(description = "Benchmarking entrypoint")

parser$add_argument(
  "--output_dir", "-o",
  dest = "output_dir", type = "character",
  help = "output directory where files will be saved",
  default = getwd(), required = TRUE
)
parser$add_argument(
  "--name", "-n",
  dest = "name", type = "character",
  help = "name of the module",
  required = TRUE
)
parser$add_argument(
  "--data.h5ad",
  dest = "data_path", type = "character",
  help = "input data h5ad path",
  required = TRUE
)
parser$add_argument(
  "--method_name",
  dest = "method_name", type = "character",
  help = "name of the method",
  choices = c("osca", "scrapper", "seurat"), required = TRUE
)

args <- parser$parse_args()

cargs <- commandArgs(trailingOnly = FALSE)
m <- grep("--file=", cargs)
run_dir <- dirname(gsub("--file=", "", cargs[[m]]))

seurat_r_path <- file.path(run_dir, "seurat.R")

timings_path <- file.path(args$output_dir, paste0(args$name, ".timings.json"))
leiden_path <- file.path(args$output_dir, paste0(args$name, ".leiden.tsv"))
louvain_path <- file.path(args$output_dir, paste0(args$name, ".louvain.tsv"))
pca_path <- file.path(args$output_dir, paste0(args$name, ".pca.tsv"))

sce <- readH5AD(args$data_path, use_hdf5 = TRUE, reader = "python")

if (args$method_name == "seurat") {
  source(seurat_r_path)
  seurat_data <- run_seurat(sce)
  write_json(
    seurat_data$time, timings_path,
    auto_unbox = TRUE, pretty = TRUE
  )
  write.table(
    data.frame(seurat_data$cell_ids, seurat_data$leiden), leiden_path,
    sep = "\t", quote = F, row.names = F
  )
  write.table(
    data.frame(seurat_data$cell_ids, seurat_data$louvain), louvain_path,
    sep = "\t", quote = F, row.names = F
  )
  write.table(
    seurat_data$pca, pca_path,
    sep = "\t", quote = F, row.names = F
  )
}
