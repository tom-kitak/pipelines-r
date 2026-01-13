library(SingleCellExperiment)
library(DelayedArray)
library(scrapper)

run_scrapper <- function(sce, time) {
  nthreads <- 1
  assay(sce) <- DelayedArray(assay(sce))

  #### 1. find mithocondial genes  ####
  start_time <- Sys.time()
  is.mito <- grepl("^MT-", rownames(sce))
  rna.qc.metrics <- computeRnaQcMetrics(assay(sce),
    subsets = list(mt = is.mito),
    num.threads = nthreads
  )
  rna.qc.thresholds <- suggestRnaQcThresholds(rna.qc.metrics)
  rna.qc.filter <- filterRnaQcMetrics(rna.qc.thresholds, rna.qc.metrics)
  filtered <- sce[, rna.qc.filter, drop = FALSE]
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("Find mithocondrial genes. Time Elapsed:", time_elapsed))
  time$find_mit_gene <- time_elapsed

  # filter data ####
  # nothing to do here

  # normalization ####
  start_time <- Sys.time()
  size.factors <- centerSizeFactors(rna.qc.metrics$sum[rna.qc.filter])
  normalized <- normalizeCounts(assay(filtered), size.factors)
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("Normalization. Time Elapsed:", time_elapsed))
  time$normalization <- time_elapsed
  assay(filtered, "normalized") <- normalized

  # Identification of highly variable features (feature selection) ####
  start_time <- Sys.time()
  gene.var <- modelGeneVariances(assay(filtered, "normalized"), num.threads = nthreads)
  hvg.sce.var <- chooseHighlyVariableGenes(gene.var$statistics$residuals, top = 1000)
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("Identification of highly variable features. Time Elapsed:", time_elapsed))
  time$hvg <- time_elapsed

  # PCA ####
  start_time <- Sys.time()
  pca <- runPca((assay(filtered, "normalized")[hvg.sce.var, ]), num.threads = nthreads, number = 50)
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("PCA. Time Elapsed:", time_elapsed))
  time$pca <- time_elapsed
  reducedDim(filtered, "PCA") <- t(pca$components)
  dim(reducedDim(filtered, "PCA")[, 1:2])
  pca_var_pct_scrap <- pca$variance.explained

  # t-sne ####
  start_time <- Sys.time()
  tsne.out <- runTsne(pca$components, num.threads = nthreads)
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("t-sne. Time Elapsed:", time_elapsed))
  time$t_sne <- time_elapsed
  reducedDim(filtered, "TSNE") <- tsne.out

  # umap ####
  start_time <- Sys.time()
  set.seed(1000000)
  umap.out <- runUmap(pca$components, num.threads = nthreads)
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("UMAP. Time Elapsed:", time_elapsed))
  time$umap <- time_elapsed
  reducedDim(filtered, "UMAP") <- umap.out

  # louvain  ####
  start_time <- Sys.time()
  snn.graph <- buildSnnGraph(pca$components, num.threads = nthreads)
  clust.out <- clusterGraph(snn.graph, method = c("multilevel"), multilevel.resolution = 0.18)
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("Louvain clusterings. Time Elapsed:", time_elapsed))
  time$louvain <- time_elapsed
  louvain_clustering <- clust.out$membership

  # leiden ####
  start_time <- Sys.time()
  snn.graph <- buildSnnGraph(pca$components, num.threads = nthreads)
  clust.out <- clusterGraph(snn.graph, method = c("leiden"), leiden.resolution = 0.16)
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("Leiden clusterings. Time Elapsed:", time_elapsed))
  time$leiden <- time_elapsed
  leiden_clustering <- clust.out$membership

  return(list(
    pca = reducedDim(filtered, "PCA"),
    hvgs = hvg.sce.var,
    cell_ids = colnames(filtered),
    time = time,
    leiden = leiden_clustering,
    louvain = louvain_clustering
  ))
}
