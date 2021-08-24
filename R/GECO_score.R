#' Score predefined clusters using GECO.
#'
#' This function takes a set of user-defined clusters and a set of ground truth
#' genes and scores them using the GECO metric. The gene names within the
#' clusters and the gene names within the ground truth gene sets must match -
#' identical naming conventions must be used for both.
#'
#' @param clusters A data.frame with two columns: the first column is the name
#' of a gene and the second is the cluster label/name for that gene.
#' @param gt_dir The absolute path to the directory containing the ground truth
#' gene sets.
#' @return A table with two columns: the first column is the name of the ground
#' truth set and the second is the GECO score for that ground truth set.
#' @examples
#' # Load clusters
#' ex_clusters <- clusters_df
#' # Define ground truth gene set directory
#' GT_dir <- '/path/to/gt_gene_sets'
#' # Score the clusters
#' clusters <- GECO_score(ex_clusters, GT_dir)
#' @export
GECO_score <- function(clusters, gt_dir) {
  gt_set_files <- list.files(path = gt_dir)
  res <- list()

  # Check all gene sets in the given directory
  for(gt_set in gt_set_files) {
    # Gather all genes in the current gene set
    gt_genes <- read.csv(paste0(gt_dir, gt_set), stringsAsFactors = FALSE)[[1]]

    # Find the number of clusters from the data
    num_clusters <- length(unique(clusters[,2]))

    # Containers for scoring information
    gene_name <- vector(mode = "list", num_clusters)
    pos_vec <- vector(mode = "list", num_clusters)
    prob_vec <- vector(mode = "list", num_clusters)

    # Apply the GECO metric
    for(i in seq_len(num_clusters)) {
      clust_size <- sum(clusters$Cluster == i)
      gene_name[[i]] <- clusters[clusters$Cluster == i,]$Genes
      # If cluster is not a singlet
      if(clust_size>1) {
        pos_vec[[i]] <- gene_name[[i]] %in% gt_genes
        gt_bf <- ( (length(gt_genes) - pos_vec[[i]]) / (clust_size - 1) )
        prob_vec[[i]] <- ( (sum(pos_vec[[i]]) - pos_vec[[i]]) / (clust_size - 1) ) / gt_bf
      }
      # If cluster is a singlet
      else {
        pos_vec[[i]] <- gene_name[[i]] %in% gt_genes
        prob_vec[[i]] <- length(gt_genes) / (clust_size - 1)
      }
    }

    # Store the ground truth gene state of each gene (T/F) and GECO metric value
    pos_vec <- unlist(pos_vec)
    prob_vec <- unlist(prob_vec)

    # Get the AUC value based on the GECO metric scores
    pred <- ROCR::prediction(prob_vec, pos_vec)
    GECO_score <- ROCR::performance(pred, "auc")@y.values[[1]]

    # Store the GECO scores for each ground truth set
    res[[tools::file_path_sans_ext(gt_set)]] <- GECO_score
  }
  # Order the GECO results in decreasing order
  res <- res[order(unlist(res), decreasing=TRUE)]

  return(res)
}
