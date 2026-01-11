library(SingleCellExperiment)
library(scuttle)
library(AnnotationDbi)
library(scran)
library(scater)
library(bluster)
library(EnsDb.Hsapiens.v75)

run_osca <- function(sce, time) {
  #### 1. find mithocondial genes  ####
  start_time <- Sys.time()
  chr.loc <- mapIds(EnsDb.Hsapiens.v75,
    keys = rownames(sce),
    keytype = "SYMBOL", column = "SEQNAME"
  )
  is.mito <- which(chr.loc == "MT")
  df <- perCellQCMetrics(sce, subsets = list(Mito = is.mito))
  # include them in the object
  colData(sce) <- cbind(colData(sce), df)
  reasons <- perCellQCFilters(df, sub.fields = "subsets_Mito_percent")
  sce$discard <- reasons$discard
  sce <- sce[, !sce$discard]
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("Find mithocondial genes. Time Elapsed:", time_elapsed))
  time$find_mit_gene <- time_elapsed

  # 2. filter data ####
  # nothing to do here

  # normalization ####
  start_time <- Sys.time()
  sce <- logNormCounts(sce)
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("Normalization. Time Elapsed:", time_elapsed))
  time$normalization <- time_elapsed

  # Identification of highly variable features (feature selection) ####
  start_time <- Sys.time()
  dec.sce <- modelGeneVar(sce)
  fit.sce <- metadata(dec.sce)
  hvg.sce.var <- getTopHVGs(dec.sce, n = 1000)
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("Identification of highly variable features. Time Elapsed:", time_elapsed))
  time$hvg <- time_elapsed

  # PCA ####
  counts(sce) <- as.matrix(counts(sce))
  logcounts(sce) <- as.matrix(logcounts(sce))
  start_time <- Sys.time()
  sce <- runPCA(sce, subset_row = hvg.sce.var)
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("PCA.Time Elapsed:", time_elapsed))
  time$pca <- time_elapsed

  # t-sne ####
  start_time <- Sys.time()
  set.seed(100000)
  sce <- runTSNE(sce, dimred = "PCA")
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("t-sne. Time Elapsed:", time_elapsed))
  time$t_sne <- time_elapsed

  # umap ####
  start_time <- Sys.time()
  set.seed(1000000)
  sce <- runUMAP(sce, dimred = "PCA")
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("umap. Time Elapsed:", time_elapsed))
  time$umap <- time_elapsed

  # louvain  ####
  start_time <- Sys.time()
  louvain_clustering <- clusterCells(sce,
    use.dimred = "PCA",
    BLUSPARAM = NNGraphParam(k = 50, cluster.fun = "louvain")
  )
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("Louvain clusterings. Time Elapsed:", time_elapsed))
  time$louvain <- time_elapsed

  # leiden ####
  start_time <- Sys.time()
  leiden_clustering <- clusterCells(sce,
    use.dimred = "PCA",
    BLUSPARAM = NNGraphParam(k = 50, cluster.fun = "leiden")
  )
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("Leiden clusterings. Time Elapsed:", time_elapsed))
  time$leiden <- time_elapsed

  return(list(
    pca = reducedDim(sce, "PCA"),
    cell_ids = colnames(sce),
    time = time,
    leiden = leiden_clustering,
    louvain = louvain_clustering
  ))
}
