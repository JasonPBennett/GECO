# GECO
GECO: A metric to assess the biological significance of gene clusters

This is still a WIP and primarily intended for in-house testing at the moment.


# Getting Started

The GECO package contains three functions:

1. *make_clusters()*
2. *generate_scores()*
3. *assess_clusters()*

## 1. make_clusters()
This function accepts a data.frame as an arguement. This dataframe should contain your counts table with genes in rows and samples in columns. The reads should already be normalized and ready for clustering. In the case of our data, this meant:

- removing genes with 0 reads across all samples
- TPM normalization of the reads
- filtering out genes with less than 1 TPM expression
- performing Z-score normalization on the reads (or any other method of feature scaling)

This function also has a number of optional parameters that can be reviewed through the man page accessed in the usual way: ?make_clusters

The output will be a set of nested lists, containing each iteration of the clustering performed, and within that will be the kmeans objects generated for each k value selected. The number of iterations, range of k-values, and total number of different k-values to try are all optional parameters as detailed in the man page for the function.

**NOTE: This step is the longest in the process and with default parameters may take about an hour. This time is variable depending on the size of the data as well as well as additional optional parameters selected. Details can be found in the man page. I highly suggest saving the output from this function to allow rapid re-running of subsequent functions.**

## 2. generate_scores()

This function takes the clusters that were just generated as output from make_clusters() in the form of a list as well as the full path to a directory containing your ground truth genes. The files contained within this directory should be csv files formatted with only the name of the ground truth set on the first line and the gene ids (listed one per line).

**NOTE: The naming convention of the genes in the original counts table MUST match the gene ids provided in the ground truth gene csv files.**

The output is a list containing the scores to be used in the final step of assessing the clusters.

## 3. assess_clusters()

This function takes the GECO scores calculated within generate_scores() and performs an ROC/AUC analysis to identify cluster quality.

The output will be a ggplot2 figure that can be examined by running the base plot() function on the returned ggplot2 object.


# Interpreting the Output

The output figure will be a graph containing:

- the k-value (number of clusters) on the x-axis
- the cluster quality on the y-axis (0.5 is random distribution of ground truth genes across the clusters, 1.0 is perfect clustering which means each ground truth gene set was placed into its own cluster with nothing but other members of the same ground truth gene set contained within)
- boxes representing the cluster quality scores for the 'n' iterations of clustering (e.g. default is 10 iterations, so the box contains the 10 quality scores for each k value). The horizontal line within the boxes indicated the mean quality for that value of k.

The take-away is the mean value within the boxes. A trend should be observed, increasing from lower k-values (10, 12, etc.) and increasing to a plateau. The selection of which k-value starts this plateau is heuristic, but the k-value at the beginning of the plateau would represent an optimal number of clusters for your dataset. With this information, you can generate clusters from your dataset using the given k-value and be confident that the clustering has captured biologically significant groups of genes.


# Quick Example

*// Normalized, feature scaled data*
df

*// Create the clusters from the data*
clusters <- GECO::make_clusters(df)

*// Find the GECO scores using the clusters and co-expressed ground truth gene sets*
dir <- "./R/gt_sets/"
GECO_scores <- GECO::generate_scores(clusters, dir)

*// Determine biological significance of clusters from the observed cluster quality*
fig <- GECO::assess_clusters(GECO_scores)
plot(fig)


# Manuscript

IN-PROGRESS
