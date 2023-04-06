# NSCLC

# Spatial analysis of human lung cancer reveals organized immune hubs enriched for stem-like CD8 T cells and associated with immunotherapy response

We examined immunity hubs in human pre-immunotherapy lung cancer specimens using multiplex fluorescence imaging.  Two staining panels were analyzed:  

1.  RNA ISH targeting CXCL13, CXCL10/11, CD3E, IFNG, panCK, dapi
2.  IHC IF targeting CD8, panCK, TCF7, Ki67, PD-L1, PD-1, dapi

# Grid windows

grid_windows.R

A computationally fast window approach was used to explore the local tissue microenvironment, where cell phenotypes were enumerated within a 50 x 50μm window.

# Unbiased window identification

kmeans_clustering.R

An unbiased K-means clustering was applied across all windows using the CXCL10/11+ cell fraction (of all cells within a window). Windows that were enriched in CXCL10/11+ cells grouped into cluster 2 (i.e. immunity hub windows), while windows that had little to no CXCL10/11+ cells were assigned to cluster 1. An analogous approach was taken for PanCK-TCF7+ cells. Immunity hub windows were overlaid onto a map of the tissue allowing visual inspection of these regions.

# Aggregation of windows

aggregate_windows.R

Our windowing approach divided large areas enriched in CXCL10/11+ cells into multiple immunity hub windows. Therefore, adjacent CXCL10/11+ windows were combined to facilitate a compositional analysis more reflective of the integrative biology across the entirety of the immunity hub. This allows us to capture a community of interacting cells that span multiple grid windows. An immunity hub is defined as two or more CXCL10/11+ windows that are adjacent and each had greater than 7 cells (15 cells for the two coregistered images). Each hub was assigned a unique identifier. For each immunity hub, cell phenotypes were enumerated and the area was computed. For each patient sample, immunity hubs were mapped back onto a plot of the tissue section for visual inspection. A similar aggregation approach was taken for defining TCF7+ aggregates. An expanded search window (100um) was employed to fully agglomerate large contiguous aggregates. Each aggregate was required to be composed of two or more PanCK-TCF7+ neighboring windows and have a minimum of 35 cells. 

Paired plots were generated comparing phenotype density in immunity hubs compared to all tumor area (both immunity hub and non-immunity hub). Benjamini-Hochberg adjusted p-values for the paired Wilcoxon test were calculated (4 phenotypes compared in RNA panel and 7 in antibody panel).

# T-SNE and Leiden clustering

tsne_Leiden.Rmd

Unsupervised clustering analysis was performed on the hubs using the Rtsne package (approach adapted from Ref58 in manuscript). Cell phenotype counts per hub area were scaled based on the global maximum for each phenotype and used to create tSNE coordinates. 21 cell phenotype features were used in the analysis: CD3E+, CD8+, CD8+Ki67+, CD8+PD-1+, CD8+PD-1+Ki67+, CD8+PD-1+TCF7+, CD8+TCF7+, CD8+TCF7+Ki67+, PanCK+ (IF only panel), PanCK+PD-L1+, PanCK+TCF7+, PanCK-PD-L1+, PanCK-TCF7+, PanCK+ (RNA ISH/IF panel), PanCK+CXCL10/11+, PanCK-CXCL10/11+, CXCL10/11+, IFNG+, PD-1+, PD-L1+, and total cell density. The optimal number of principal components for variance explained was identified as 6 by scree plot. Thereafter a tSNE (t-distributed Stochastic Neighbor Embedding) dimensionality reduction was performed using Rtsne, with perplexity set to 50. To identify subclusters, the tSNE coordinates were used to create an edge list for a fast k-nearest neighbor search algorithm with the maximum number of nearest neighbors equal to 50. Next, distances were converted to weights such that greater weight was given to smaller neighbor distances. A Leiden clustering algorithm was then applied using the modularity method with the resolution parameter set to 0.25, yielding 7 Leiden subclusters. The subclusters were mapped back onto a tissue plot for closer inspection of cell phenotypes in these regions. Plots were also generated showing subclusters present within each patient image and the cell phenotype composition per subcluster. P-values were adjusted for multiple hypothesis testing using the Benjamini-Hochberg method.

