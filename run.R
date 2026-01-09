#!/usr/bin/env Rscript

library(argparse)

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
