library(Seurat)
library(SingleCellExperiment)

run_seurat <- function(sce, time) {
  # data ####
  data <- as.Seurat(sce, counts = "counts", data = NULL, assay = NULL)
  DefaultAssay(data) <- "originalexp"

  # find mitocondrial genes ####
  start_time <- Sys.time()
  data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("Find mitocondrial genes. Time Elapsed:", time_elapsed))
  time$find_mit_gene <- time_elapsed

  # filter data ####
  write(paste0("before: ", dim(data)), stderr())
  start_time <- Sys.time()
  qc <- metadata(sce)$qc_thresholds
  data <- subset(
    data,
    subset = nFeature_originalexp > qc[qc$metric == "nFeature", "min"] &
      nFeature_originalexp < qc[qc$metric == "nFeature", "max"] &
      percent.mt < qc[qc$metric == "percent.mt", "max"] &
      nCount_originalexp < qc[qc$metric == "nCount", "max"]
  )
  end_time <- Sys.time()
  write(paste0("after: ", dim(data)), stderr())
  time_elapsed <- end_time - start_time
  print(paste("Filter data. Time Elapsed:", time_elapsed))
  time$filter <- time_elapsed

  # normalization ####
  start_time <- Sys.time()
  data <- NormalizeData(data, normalization.method = "LogNormalize")
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("Normalization. Time Elapsed:", time_elapsed))
  time$normalization <- time_elapsed

  # Identification of highly variable features (feature selection) ####
  start_time <- Sys.time()
  data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 1000)
  # top10 <- head(VariableFeatures(data), 10)
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("Identification of highly variable features. Time Elapsed:", time_elapsed))
  time$hvg <- time_elapsed

  # Scaling the data ####
  start_time <- Sys.time()
  all.genes <- rownames(data)
  data <- ScaleData(data, features = all.genes)
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("Scaling the data. Time Elapsed:", time_elapsed))
  time$scaling <- time_elapsed

  # PCA ####
  start_time <- Sys.time()
  data <- RunPCA(
    data,
    features = VariableFeatures(object = data), verbose = FALSE
  )
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("PCA. Time Elapsed:", time_elapsed))
  time$pca <- time_elapsed

  # t-sne ####
  start_time <- Sys.time()
  data <- RunTSNE(data,
    reduction = "pca", perplexity = 18
  )
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("t-sne. Time Elapsed:", time_elapsed))
  time$t_sne <- time_elapsed

  # UMAP ####
  start_time <- Sys.time()
  data <- RunUMAP(data, dims = 1:50)
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("UMAP. Time Elapsed:", time_elapsed))
  time$umap <- time_elapsed

  # louvain ####
  start_time <- Sys.time()
  data <- FindNeighbors(data, dims = 1:50, verbose = T)
  data <- FindClusters(
    data,
    resolution = 0.1, algorithm = 1, cluster.name = "louvain"
  )
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("Louvain Clustering. Time Elapsed:", time_elapsed))
  time$louvain <- time_elapsed

  # leiden ####
  start_time <- Sys.time()
  # data <- FindNeighbors(data, dims = 1:50, verbose = T)
  data <- FindClusters(
    data,
    algorithm = 4, resolution = 0.08, cluster.name = "leiden"
  )
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("Leiden Clustering. Time Elapsed:", time_elapsed))
  time$leiden <- time_elapsed

  return(list(
    pca = Embeddings(data, reduction = "pca"),
    hvgs = VariableFeatures(data),
    time = time,
    cell_ids = colnames(data),
    leiden = data$leiden,
    louvain = data$louvain
  ))
}
