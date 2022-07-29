<!-- badges: start -->
[![](https://www.r-pkg.org/badges/version/CIARA)](https://cran.r-project.org/package=CIARA)
[![](http://cranlogs.r-pkg.org/badges/grand-total/CIARA?color=green)](https://cran.r-project.org/package=CIARA)
[![](http://cranlogs.r-pkg.org/badges/CIARA?color=green)](https://cran.r-project.org/package=CIARA)
[![](http://cranlogs.r-pkg.org/badges/last-day/CIARA?color=green)](https://cran.r-project.org/package=CIARA)


<!-- badges: end -->

- [CIARA](#ciara)
  * [Installation](#installation)
  * [Getting started](#getting-started)
    + [CIARA_gene](#CIARA-gene)
    + [CIARA](#CIARA)
  * [Visualization of rare cell markers](#visualization-of-rare-cell-markers)
  * [Cluster analysis based on CIARA for rare cell types identification](#cluster-analysis-based-on-CIARA-for-rare-cell-types-identification)
  * [Comparison with other published methods for rare cell types identification](#comparison-with-other-published-methods-for-rare-cell-types-identification)
  * [Vignette](#vignette)
  * [Citation](#Citation)
  * [Contributions and Support](#contributions-and-support)
  * 
# CIARA
CIARA (Cluster Independent Algorithm for the identification of markers of RAre cell types) is an R package that identifies potential markers of rare cell types looking at genes whose expression is confined in small regions of the expression space. 
It is possible to use these highly localized genes as features in standard cluster algorithm (i.e. Louvain), for identifying extremely rare population(3/4 cells from thousand of cells)

## Installation

You can install the released version of CIARA from [CRAN](https://CRAN.R-project.org) with:

```r
install.packages("CIARA")
```

And the development version from [GitHub](https://github.com/) with:

```r
devtools::install_github("ScialdoneLab/CIARA",auth_token="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",ref="master")
```

## Getting started 
The main functions of the package are **CIARA_gene** and **CIARA**


### CIARA_gene

```r
CIARA_gene(norm_matrix, knn_matrix, gene_expression, p_value = 0.001, local_region = 1, approximation)
```
requires as input:

1. **norm_matrix**: Norm count matrix (n_genes x n_cells)
2. **knn_matrix**: K-nearest neighbors matrix (n_cells x n_cells)
3. **gene_expression**: numeric vector with the gene expression (length equal to n_cells). The gene expression is binarized (equal to 0/1 in the cells where the value is below/above the median)
4. **p_value**: maximum p value (returned by the R function fisher.test with parameter alternative = "g") for considering a local region enriched 
5. **local_region**: minimum number of local regions (cell with its knn neighbours) where the binarized gene expression is enriched in 1
6. **approximation**.Logical.If TRUE, for a given gene the fisher test is run in the local regions of only the cells where the binarized gene expression is 1

The gene expression is binarized (1/0) if the value in a given cell is above/below the median. Each of cell with its first K nearest neighbors defined a local region. If there are at least **local_region** enriched in 1 according **to fisher.test** (with p value below than **p_value** and odds ratio above or equal to **odds_ratio**) , then the entropy for the gene is computed starting from the probability of having 1/0. The minimum of the entropy across all the enriched local regions is the entropy of mixing. If there are no enriched local regions, then the entropy of mixing  and the p value by default are set to 1
The output of **CIARA_gene**  is a list with one element corresponding to the p value of the gene

### CIARA

```r
CIARA(norm_matrix, knn_matrix, background, cores_number = 1, p_value = 0.001,local_region = 1, approximation)
```
requires as input:

1. **norm_matrix**: Norm count matrix (n_genes x n_cells)
2. **knn_matrix**: K-nearest neighbors matrix (n_cells x n_cells)
3. **background**: Vector of genes for which the function **CIARA_gene** is run
4. **cores_number**: Integer.Number of cores to use.
5. **p_value**: maximum p value (returned by the R function fisher.test with parameter alternative = "g") for considering a local region enriched 
6. **local_region**: minimum number of local regions (cell with its knn neighbours) where the binarized gene expression is enriched in 1
7. **approximation**.Logical. If TRUE, for a given gene the fisher test is run in the local regions of only the cells where the binarized gene expression is 1

Return a dataframe with n_rows equal to the length of **background** . Each row is the output from **CIARA_gene**.

The vector of genes for which the function **CIARA_gene** is run can be obtained with the function **get_background_full**.
This function gives as output a vector with all genes expressed at a level higher than **threshold** in a number of cells between **n_cells_low** and **n_cells_high**

Below an example of input using the development version of **CIARA** from GitHub
```r
load(file = "raw_counts_human_data.Rda")
human_data_seurat <- cluster_analysis_integrate_rare(raw_counts_human_data, "Human_data", 0.1, 5, 30)
norm_matrix <- as.matrix(GetAssayData(human_data_seurat, slot = "data",assay="RNA"))
knn_matrix <- as.matrix(human_data_seurat@graphs$RNA_nn)
background <- get_background_full(norm_matrix, threshold = 1, n_cells_low = 3, n_cells_high = 20)
result <- CIARA(norm_matrix, knn_matrix, background, cores_number = 1, p_value = 0.001, local_region = 1, approximation = FALSE)
```
The input file *raw_counts_human_data.Rda* can be downloaded [here](https://hmgubox2.helmholtz-muenchen.de/index.php/s/x83jDLHobM7Qer6)


In the next two sections (**Visualization of rare cell markers** and **Cluster analysis based on CIARA for rare cell types identification**) it is shown the analysis of the scRNA seq data from human embryo at the gastrulation state from [Tyser *et al.*, 2021](https://www.nature.com/articles/s41586-021-04158-y)

### Visualization of rare cell markers

 We can visualize the highly localized genes identified with CIARA with the functions **plot_gene**, **plot_genes_sum** and even in an interactive way with the function 
**plot_interactive**. For more exhaustive information about the functions offered by CIARA for visualization  see **Tutorials section** below and the help page of the single functions. (*?function_name*).
The pattern expression of the top two genes according to CIARA are shown.
```r
load(system.file("extdata", "result.Rda", package = "CIARA"))
ciara_genes <- row.names(result)[result[, 1] < 1]
ciara_genes_top <- row.names(result)[order(as.numeric(result[, 1]))]
coordinate_umap <- as.data.frame(Embeddings(human_data_seurat, reduction = "umap")[, 1:2])
# In the example below we keep the same umap coordinate used in the original paper
meta_info <- readRDS(system.file("extdata", "annot_umap.rds", package ="CIARA"))
coordinate_umap <- meta_info[,2:3]

p=list()
for(i in ciara_genes_top[1:2]){
  q <- plot_gene(norm_matrix, coordinate_umap, i, i)
  p <- list(p,q)
}
p

```
<img src="https://github.com/ScialdoneLab/CIARA/blob/main/figures/entropy_new_new1.png" width="700" height="500">
<img src="https://github.com/ScialdoneLab/CIARA/blob/main/figures/entropy_new_new_2.png" width="700" height="500">

We can visualize which are the rare cell markers expressed in a particular region of the umap plot with the function **plot_localized_genes** and in an interactive way with **plot_localized_genes_interactive**.
Each cell is coloured according to the number of expressed rare cell markers shared with neighboring cells.
These two functions are currently only present in the development version of **CIARA**.
```r
localized_genes_human <-  detect_localized_genes(knn_matrix,norm_matrix,ciara_genes_top,100)  
list_intersect <-  localized_genes_human[[1]]
rank_intersect <-  localized_genes_human[[2]]
```

```r
genes_name_text <-  names_localized_genes(list_intersect ,ciara_genes_top,max_number = 5)
```

Below the vector of colors is defined.
```r
ramp <- colorRamp(c("white", "blue4"))
ramp.list <- rgb( ramp(seq(0, 1, length = length(unique(rank_intersect)))), max = 255)

index_color=round(length(ramp.list)/2,0)
breaks <- seq(0,max(rank_intersect),length.out=1000)
# library gplots must be installed for executing the following command
gradient1 <- gplots::colorpanel( sum( breaks[-1]<= as.numeric(quantile(breaks,0.15))), "#FFFFFF",ramp.list[index_color])
gradient2 <- gplots::colorpanel( sum( breaks[-1] > as.numeric(quantile(breaks,0.15)) ), ramp.list[index_color], "#00008B" )
hm.colors <- c(gradient1,gradient2)
```

```r
plot_localized_genes(coordinate_umap,norm_matrix,rank_intersect,"Top genes CIARA",hm.colors)

```


<img src="https://github.com/ScialdoneLab/CIARA/blob/main/figures/plot_localized_genes_example.png" width="700" height="500">


```r
plot_localized_genes_interactive(coordinate_umap,norm_matrix,rank_intersect,genes_name_text,hm.colors,min_x=NULL,max_x=NULL,min_y=NULL,max_y=NULL)
```

<img src="https://github.com/ScialdoneLab/CIARA/blob/main/figures/plot_localized_genes_interactive_example.png" width="700" height="500">

### Cluster analysis based on CIARA for rare cell types identification

We can use the genes identified by CIARA as features in standard algorithm (i.e. Louvain) for the identification of extremely rare population of cells (3/4 cells from dataset of thousand of cells).
This approach consists of four steps:
1. **general cluster assignment**: A standard cluster analysis on all the dataset is performed. It is possible to use the function **cluster_analysis_integrate_rare** based on Seurat functions (RunPCA,FindNeighbors,FindClusters) with HVGs as features or an existing cluster assignment available for the dataset.
2. **Identification of rare population**: A standard cluster analysis (with the function **cluster_analysis_integrate_rare**) is performed on all the dataset using as features the highly localized genes provided by the function **CIARA** 
3. **Merging of step 1 and 2**: Any cluster identified at step 2 of size  smaller than **max_number** constitutes an independent cluster.
4. **Identification of markers for small clusters**: It is possible to check if the rare cell types are well defined according to their markers with the function **white_black_markers**. See ?**white_black_markers** for more information.
5. **Identification of cell sub-types**: For each of the cluster identified at step 3, a Fisher test is performed to see if there is a statistically significant enrichment between the HVGs in the cluster and the highly localized genes (function **test_hvg**). For all the clusters where there is an enrichment, a sub cluster analysis starting from the original cluster is performed using as features  HVGs (function **cluster_analysis_sub**)

An example is the following:
```r
# step 1
human_embryo_analysis <- cluster_analysis_integrate_rare(raw_counts = raw_counts_human_data, project_name = "Human_Embryo_data", resolution = 0.1, neighbors=5, max_dimension = 30)
human_embryo_cluster <- as.vector(human_embryo_analysis$seurat_clusters)
# or we can also use the original cluster annotation provided by Tyser et al,2020
original_cluster <- as.vector(meta_info$cluster_id)
human_embryo_cluster <- original_cluster
```

```r
# step 2
human_embryo_analysis_ciara <- cluster_analysis_integrate_rare(raw_counts_human_data, "Human_Embryo_data", 0.01, 5, 30, ciara_genes)
```
```r
# step 3
final_cluster <- merge_cluster(original_cluster, human_embryo_analysis_ciara$seurat_clusters, max_number = 20)
final_cluster[grep("step_2",final_cluster)] <- "PGC"
```
With CIARA is possible to obtain in a completely unsupervised way the primordial germ cells cluster (PGCs), that in the original paper ([Tyser *et al.*, 2021](https://www.nature.com/articles/s41586-021-04158-y)) was detected only with a supervised approach.

```r
# step 4
result_test <- test_hvg(raw_counts_human_data,final_cluster, ciara_genes, background, number_hvg = 100, min_p_value = 0.001)
result_test[[2]]
 "Endoderm"  "Hemogenic Endothelial Progenitors"
# We need to do sub cluster in the above two clusters
raw_endoderm <- raw_counts_human_data[, human_embryo_cluster == "Endoderm"]
raw_emo <- raw_counts_human_data[, human_embryo_cluster == "Hemogenic Endothelial Progenitors"]
combined_endoderm <- cluster_analysis_sub(raw_endoderm, 0.2, 5, 30, "Endoderm")
combined_emo <- cluster_analysis_sub(raw_emo, 0.6, 5, 30, "Hemogenic Endothelial Progenitors")


all_sub_cluster <- c(combined_endoderm$seurat_clusters, combined_emo$seurat_clusters)
names(final_cluster) <- colnames(raw_counts_human_data)
final_cluster_version_sub <- merge_cluster(final_cluster, all_sub_cluster)
```

```r
plot_umap(coordinate_umap, final_cluster_version_sub)
```
<img src="https://github.com/ScialdoneLab/CIARA/blob/main/figures/entropy_cluster.png" width="700" height="500">
With the cluster analysis based on CIARA we are able to detect two clusters (Endoderm_2 and and Hemogenic Endothelial Progenitors_4 highlighted in the plot) that were not reported in the original paper.


For more exhaustive information about the functions offered by CIARA for the identification of rare populations of cells  see **Tutorials section** below and the help page of the single functions. (*?function_name*).

## Comparison with other published methods for rare cell types identification

We use the [Matthews Correlation Coefficient (MCC)](https://en.wikipedia.org/wiki/Phi_coefficient) to evaluate the performance of other methods in the task of identfying the PGCs cluster from human gastrula. 
CIARA is the only able to detect the original PGCs cluster (7 cells).
<img src="https://github.com/ScialdoneLab/CIARA/blob/main/figures/mcc_gap_clust.png" width="700" height="500">


## Vignette

The following vignette is available and completely reproducible. It uses single cell RNA seq from human embryo at the gastrulation state from [Tyser *et al.*, 2021](https://www.nature.com/articles/s41586-021-04158-y). The raw count matrix was downloaded from  [http://human-gastrula.net].
An extremely rare population of primordial germ cells (PGCs-7 cells) is easily identified with CIARA.
It can be accessed within R with:
```r
utils::vignette("CIARA")
```

## Citation
To cite package **CIARA** in publications use:
```
Gabriele Lubatti (2022). CIARA: Cluster Independent Algorithm for Rare Cell Types
Identification. R package version 0.1.0. https://CRAN.R-project.org/package=CIARA

A BibTeX entry for LaTeX users is

  @Manual{,
    title = {CIARA: Cluster Independent Algorithm for Rare Cell Types Identification},
    author = {Gabriele Lubatti},
    year = {2022},
    note = {R package version 0.1.0},
    url = {https://CRAN.R-project.org/package=CIARA},
  }

```

## Contributions and Support
Contributions in the form of feedback, comments, code and bug report are welcome.
* For any contributions, feel free to fork the source code and [submit a pull requests](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request-from-a-fork).
* Please report any issues or bugs here: https://github.com/ScialdoneLab/CIARA/issues.
Any questions and requests for support can also be directed to the package maintainer (gabriele[dot]lubatti[at]helmholtz-muenchen[dot]de).

