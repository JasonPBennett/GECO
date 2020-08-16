# GECO
GECO: A metric to assess the biological significance of gene clusters

This is still a WIP and primarily intended for in-house testing at the moment.

The goal of this tool is generate sets of gene clusters that can be evauluated using the GECO metric to determine the highest quality clusters. The level of quality is determined by the distribution of known co-expressed gene sets, with the assumption that genes that are co-expressed should be clustered together if our clustering parameters were correct. The GECO package takes a dataset in the form of a normalized counts table and generates a series of clusters. These clusters are created using the k-means algorithm with wide-ranging k-values (or numbers of clusters). Each permutation is then scored using the GECO metric and ultimately given cluster quality scores with the goal of identifying a certain k-value, or a certain optimal number of clusters, that best captures the delicate biological signals found within the dataset. Following the use of the package, an optimal set of clusters can be generated based on the optimal k-value indicated by the package.


# Getting Started

The GECO package contains three functions:

1. *generate_clusters()*
2. *score_clusters()*
3. *assess_quality()*

## 1. generate_clusters()
This function accepts a dataframe as an argument. This dataframe should contain your counts table with genes in rows and samples in columns. The counts within the dataframe should already be normalized and ready for clustering. In the case of our data, this meant:

- removing genes with 0 reads across all samples
- TPM normalization of the reads
- filtering out genes with less than 1 TPM expression
- filtering out non-protein coding genes
- performing Z-score normalization on the reads (or any other method of feature scaling)

This function also has a number of optional parameters that can be reviewed through the man page accessed in the usual way: ?generate_clusters

The clustering process will be performed iteratively, due to the stochasticity of start conditions used by the k-means algorithm. Multiple iterations allow the package to control for aberrant start conditions. The output will be a set of nested lists, containing each iteration of the clustering performed, and within that will be the kmeans objects generated for each k value selected. The number of iterations, range of k-values, and total number of different k-values to try are all optional parameters as detailed in the man page for the function.

```
A visual example of the output follows, assuming default parameters were used:
- 'clusters' is a list which is the output of the generate_clusters() function
- within 'clusters' will be another list numbered 1 - 10: these refer to the iterations
-- the visual example assumes iteration '1' is being investigated
- within each 'iteration' is another list containing all kmeans clusters, with one entry per k-value
-- the visual example assumes that we are looking at iteration '1' and the kmeans object generated with k = 10

clusters__________________________________________________________________________
|  
|  [1] iterations_________________________________________________________________
|      |
|      |['10']kmeans_objects______________________________________________________
|      |      |
|      |      |  kmeans(data, k=10, ...)
|      |      |

Ultimately, the visual example shows the kmeans object from the first iteration when clustering with a k-value of 10.

```

**NOTE: This step is the longest in the process and with default parameters may take about an hour. This time is variable depending on the size of the data as well as well as additional optional parameters selected. Details can be found in the man page. I highly suggest saving the output from this function to allow rapid re-running of subsequent functions.**

## 2. score_clusters()

This function takes the clusters that were just generated as output from *generate_clusters()* in the form of a list as well as the full path to a directory containing your ground truth genes. The files contained within this directory should be csv files formatted with only the name of the ground truth set on the first line and the gene ids (listed one per line).

**NOTE: The naming convention of the genes in the original counts table MUST match the gene ids provided in the ground truth gene csv files. In the case of ENSEMBL ids, make sure the versions match also (the decimal portion of the id), and as a last resort strip version numbers to still allow matching.**

```
The output is a list containing the scores to be used in the final step of assessing the clusters. Again, a visual representation of the output object is provided:
- 'scores' is a list which is the output of the score_clusters() function
- within 'scores' is another list, each nested list is labeled with the name of the ground truth gene set provided in the first line of the .csv files required by the score_clusters() function
-- the visual example assumes that a ground truth set called 'x' was selected (a terrible label!)
- within 'x' will be another list numbered 1 - 10: these refer to the iterations
-- the visual example assumes iteration '1' is being investigated
- within each 'iteration' is another list containing all kmeans scores, broken down by k-value
-- the visual example assumes a k-value of '10' is being investigated

scores_________________________________________________________________________________
|  
| ['x'] Ground_Truth_Sets______________________________________________________________
|       |
|       | [1]  iterations______________________________________________________________
|       |      |
|       |      |['10']k_value__________________________________________________________
|       |      |      |
|       |      |      |  scores table for all genes in this clustering configuration
|       |      |      |

Ultimate, the example shows the scores from the 'x' ground truth gene set within the first iteration given a k-value of 10.

```

## 3. assess_quality()

This function takes the GECO scores calculated within *score_clusters()* and performs an ROC/AUC analysis to identify cluster quality.

The output will be a ggplot2 object that can be visualized by running the base plot() function on the returned ggplot2 object.


# Interpreting the Output

The output figure will be a graph containing:

- the k-value (number of clusters) on the x-axis
- the cluster quality on the y-axis (0.5 is random distribution of ground truth genes across the clusters, 1.0 is perfect clustering which means each ground truth gene set was placed into its own cluster with nothing but other members of the same ground truth gene set contained within)
- boxes representing the cluster quality scores for the 'n' iterations of clustering (e.g. default is 10 iterations, so the box contains the 10 quality scores for each k value). The horizontal line within the boxes indicates the mean quality for that value of k.

The take-away is the mean value within the boxes. A trend should be observed, increasing from lower k-values (10, 12, etc.) and increasing to a plateau. The selection of which k-value starts this plateau is heuristic, but the k-value at the beginning of the plateau would represent an optimal number of clusters for your dataset. With this information, you can generate clusters from your dataset using the given k-value and be confident that the clustering has captured biologically significant groups of genes.

**NOTE: If a plateau in scores is not observed, consider increasing the 'kmax' parameter and re-running the *generate_clusters()* function with this larger 'kmax' value.**


# Quick Example

*// Normalized, feature scaled data*<br/>
df

*// Create the clusters from the data*<br/>
clusters <- GECO::generate_clusters(df)

*// Find the GECO scores using the clusters and co-expressed ground truth gene sets*<br/>
dir <- "./R/gt_sets/"<br/>
GECO_scores <- GECO::score_clusters(clusters, dir)

*// Determine biological significance of clusters from the observed cluster quality*<br/>
fig <- GECO::assess_quality(GECO_scores)<br/>
plot(fig)


# Manuscript

IN-PROGRESS
