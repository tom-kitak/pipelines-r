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
  "--timings.json",
  dest = "timings_path", type = "character",
  help = "output timings json path",
  required = TRUE
)
parser$add_argument(
  "--leiden.tsv",
  dest = "leiden_path", type = "character",
  help = "leiden tsv path",
  required = TRUE
)
parser$add_argument(
  "--louvain.tsv",
  dest = "louvain_path", type = "character",
  help = "louvain tsv path",
  required = TRUE
)
parser$add_argument(
  "--pca.tsv",
  dest = "pca_path", type = "character",
  help = "pca path",
  required = TRUE
)
parser$add_argument(
  "--method_name",
  dest = "method_name", type = "character",
  help = "name of the method",
  choices = c("osca", "scrapper", "seurat"), required = TRUE
)

args <- parser$parse_args()
sce <- readH5AD(args$data_path, use_hdf5 = TRUE, reader = "python")

if (args$method_name == "seurat") {
  source("seurat.R")
  seurat_data <- run_seurat(sce)
  write_json(
    seurat_data$time, args$timings_path,
    auto_unbox = TRUE, pretty = TRUE
  )
  write.table(
    data.frame(seurat_data$cell_ids, seurat_data$leiden),
    args$leiden_path,
    sep = "\t", quote = F, row.names = F
  )
  write.table(
    data.frame(seurat_data$cell_ids, seurat_data$louvain),
    args$louvain_path,
    sep = "\t", quote = F, row.names = F
  )
  write.table(
    seurat_data$pca,
    args$pca_path,
    sep = "\t", quote = F, row.names = F
  )
}
