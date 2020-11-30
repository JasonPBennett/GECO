#' Create the clusters.
#'
#' This function takes the normalized data (TPM/FPKM & feature scaled) and uses
#' the k-means function to generate an iterative series of clusters to identify
#' a potentially optimal number of clusters for the dataset.
#'
#' @param df A dataframe containing the normalized reads
#' @param kmin An integer indicating the minimum number of clusters to generate.
#' By default, this is set to 10.
#' @param kmax An integer indicating the maximum number of clusters to generate.
#' By default, this is set to 150.
#' @param ktot An integer indicating how many unique k-values to generate.
#' By default, this is set to 15. This produces 15 values ranging from kmin up
#' to kmax. Increasing this number will significantly impact performance.
#' @param num_iter An integer indicating the number or cluster iterations to
#' generate. By default, this is set to 10. This will perform the same k-means
#' clustering multiple times to account for the stochastic nature of the k-means
#' algorithm, resulting in a mean quality value in the final step that is more
#' reliable than a single iteration would be. Lowering this value will
#' negatively affect the GECO quality assessment, raising it will impact
#' performance.
#' @param km_algo A string indicating which k-means algorithm to use. By
#' default, this is set to 'Hartigan-Wong'.
#' @return A list containing each iteration of the clustering performed. Within
#' each of the iterations are the kmeans objects for use in the second step e.g.
#' score_clusters(clusters).
#' @examples
#' # Create a pseudo RNA-seq counts table
#' df <- data.frame(replicate(10,sample(-1:10,200,rep=TRUE)))
#' rownames(df) <- paste0(rep("Gene.", 200), seq(1:200))
#' # Generate clusters
#' clusters <- generate_clusters(df)
#' @export
generate_clusters <- function(df, kmin, kmax, ktot, num_iter, km_algo) {
  # Define default parameters
  if(missing(kmin)) { kmin = 10 }
  if(missing(kmax)) { kmax = 150 }
  if(missing(ktot)) { ktot = 15 }
  if(missing(num_iter)) { num_iter = 10 }
  if(missing(km_algo)) { km_algo = "Hartigan-Wong" }

  # Check that input was provided in correct format
  if(is.null(df) | !is.data.frame(df)) { stop("Input must be a data frame!") }

  # Generate logarithmically-spaced sequence
  k_vec <- round(exp(seq(log(kmin), log(kmax), length.out = ktot)))

  # Create a list to hold the clustering results
  km_clusters <- vector(mode = "list", num_iter)

  # Run kmeans repeatedly
  for(i in seq_len(num_iter)) {
    # Create list to hold the different k-value clusterings per iteration
    k_clusts <- vector(mode = "list", length(k_vec))

    for(k in k_vec) {
      # Run k-means for each k value
      clust <- stats::kmeans(df, k, iter.max = 30, algorithm = km_algo)

      # Make sure that the algorithm reached convergence
      if(!is.null(clust$ifault)) {
        # If it didn't converge, redo the clustering
        while(clust$ifault != 0) {
          message("Initial k-means clustering could not converge. Attempting to reach convergence by increasing convergence threshold parameter (iter.max = 100).")
          clust <- stats::kmeans(df, k, iter.max = 100, algorithm = km_algo)
        }
      }

      k_clusts[[k]] <- clust
    }
    # Take out the list entries that are null
    k_clusts <- Filter(Negate(is.null), k_clusts)
    names(k_clusts) <- k_vec

    # Store the kmeans iterations
    km_clusters[[i]] <- k_clusts
  }
  names(km_clusters) <- paste0("Iteration ", c(seq_len(num_iter)), sep = "")

  return(km_clusters)
}
