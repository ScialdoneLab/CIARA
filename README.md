# CIARA
CIARA (Cluster Independent Algorithm for the identification of RAre cell types) is an R package that identifies potential markers of rare cell types looking at genes whose expression is confined in small regions of the expression space. 
It is possible to use these highly localized genes as features in standard cluster algorithm (i.e. Louvain), for identifying extremely rare population(3/4 cells from thousand of cells)

## Installation

For installation, please use the following command:

```install_github("ScialdoneLab/CIARA",auth_token="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",ref="master")```.
The main function of the package is **CIARA_gene** and **CIARA**


## CIARA_gene

```CIARA_gene(norm_matrix, knn_matrix, gene_expression, p_value = 0.001, odds_ratio=2, local_region = 1, approximation = FALSE)```
requires as input:

1. **norm_matrix**: Norm count matrix (n_genes x n_cells)
2. **knn_matrix**: K-nearest neighbors matrix (n_cells x n_cells)
3. **gene_expression**: numeric vector with the gene expression (length equal to n_cells). The gene expression is binarized (equal to 0/1 in the cells where the value is below/above the median)
4. **p_value**: maximum p value (returned by the R function fisher.test with parameter alternative = "g") for considering a local region enriched 
6. **odds_ratio**: minimum odds ratio (returned by the R function **fisher.test** with parameter alternative = "g") above which a local region is considered enriched
5. **local_region**: minimum number of local regions (cell with its knn neighbours) where the binarized gene expression is enriched in 1
7. **approximation**.Logical.For a given gene, the fisher test is run in the local regions of only the cells where the binarized gene expression is 1

The gene expression is binarized (1/0) if the value in a given cell is above/below the median. Each of cell with its first K nearest neighbors defined a local region. If there are at least **local_region** enriched in 1 according **to fisher.test** (with p value below than **p_value** and odds ratio above or equal to **odds_ratio**) , then the entropy for the gene is computed starting from the probability of having 1/0. The minimum of the entropy across all the enriched local regions is the entropy of mixing. If there are no enriched local regions, then the entropy of mixing  and the p value by default are set to 1
The output of **CIARA_gene**  is a list with one element corresponding to the p value of the gene

## CIARA

```CIARA(norm_matrix, knn_matrix, background, cores_number = 1, p_value = 0.001, odds_ratio = 2,local_region = 1, approximation = FALSE) ```
requires as input:

1. **norm_matrix**: Norm count matrix (n_genes x n_cells)
2. **knn_matrix**: K-nearest neighbors matrix (n_cells x n_cells)
3. **background**: Vector of genes for which the function **CIARA_gene** is run
4. **cores_number**: Integer.Number of cores to use.
5. **p_value**: maximum p value (returned by the R function fisher.test with parameter alternative = "g") for considering a local region enriched 
6. **odds_ratio**: minimum odds ratio (returned by the R function fisher.test with parameter alternative = "g") above which a local region is considered enriched
7. **local_region**: minimum number of local regions (cell with its knn neighbours) where the binarized gene expression is enriched in 1
8. **approximation**.Logical.For a given gene, the fisher test is run in the local regions of only the cells where the binarized gene expression is 1

Return a dataframe with n_rows equal to the length of **background** . Each row is the output from **CIARA_gene**.

The vector of genes for which the function **CIARA_gene** is run can be obtained with the function **Get_background_full**.
This function gives as output a vector with all genes expressed at a level higher than **threshold** in a number of cells between **n_cells_low** and **n_cells_high**

An example of input could be:
```
load(file = "norm_counts.Rda")
load(file = "knn_matrix.Rda")
background <- get_background_full(norm_counts, threshold = 1, n_cells_low = 3, n_cells_high = 20)
result <- CIARA(norm_matrix, knn_matrix, background, cores_number = 1, min_p_value = 0.001, odds_ratio = 2, local_region = 1, approximation = FALSE)
```
The two input files *norm_counts* and *knn_matrix* can be obtained from a seurat object:
```
norm_counts <- as.matrix(GetAssayData(seurat_object, slot = "data",assay="RNA"))
knn_matrix <- as.matrix(seurat_object@graphs$RNA_nn)
```

In the next two sections (**Visualization of highly localized genes** and **Cluster analysis based on CIARA for the identification of extremely rare population of cells**) it is shown the analysis of the scRNA seq data from human embryo at the gastrulation state from [Tyser *et al.*, 2020](https://www.biorxiv.org/content/10.1101/2020.07.21.213512v1)

## Visualization of highly localized genes

 We can visualize the highly localized genes identified with CIARA with the functions **plot_gene**, **plot_genes_sum** and even in an interactive way with the function 
**plot_interactive**. For more exhaustive information about the functions offered by CIARA for visualization  see **Tutorials section** below and the help page of the single functions. (*?function_name*).
The pattern expression of the top two genes according to CIARA are shown.
```
load(system.file("extdata", "result.Rda", package = "CIARA"))
ciara_genes <- row.names(result)[result[, 1] < 1]
geni_top <- row.names(result)[order(as.numeric(result[, 1]))]
coordinate_umap <- as.data.frame(Embeddings(seurat_object, reduction = "umap")[, 1:2])
#In the example below we keep the same umap coordinate used in the original paper
meta_info <- readRDS(system.file("extdata", "annot_umap.rds", package ="CIARA"))
cordinate_umap <- meta_info[,2:3]

p=list()
for(i in geni_top[1:2]){
  q <- Plot_gene(norm_counts, cordinate_umap, i, i)
  p <- list(p,q)
}
p

```
<img src="https://github.com/ScialdoneLab/CIARA/blob/main/figures/entropy_new_new1.png" width="700" height="500">
<img src="https://github.com/ScialdoneLab/CIARA/blob/main/figures/entropy_new_new_2.png" width="700" height="500">

We can visualize in an interactive way which are the genes highly localized in a particular region of the umap plot with the function **plot_interactive**.
Each cell is coloured according to the sum of the normalized expression of the highly localized genes identified with CIARA.
```

norm_counts_small <- apply(norm_counts[ciara_genes, ], 1, function(x) {
  y <- x/sum(x)
  return(y)
  })
gene_sum <- apply(norm_counts_small, 1, sum)
    
    
genes_name_text <- selection_localized_genes(norm_counts, ciara_genes, min_number_cells = 4, max_number_genes = 4)
colnames(coordinate_umap) <- c("UMAP_1","UMAP_2")
plot_interactive(coordinate_umap, gene_sum, genes_name_text, min_x=NULL, max_x=NULL, min_y=NULL, max_y=NULL)
```

<img src="https://github.com/ScialdoneLab/CIARA/blob/main/figures/interactive_plot_new.png" width="700" height="500">

## Cluster analysis based on CIARA for the identification of extremely rare population of cells

We can use the genes identified by CIARA as features in standard algorithm (i.e. Louvain) for the identification of extremely rare population of cells (3/4 cells from dataset of thousand of cells).
This approach consists of four steps:
1. **general cluster assignment**: A standard cluster analysis on all the dataset is performed. It is possible to use the function **cluster_analysis_integrate_rare** based on Seurat functions (RunPCA,FindNeighbors,FindClusters) with HVGs as features or an existing cluster assignment available for the dataset.
2. **Identification of rare population**: A standard cluster analysis (with the function **cluster_analysis_integrate_rare**) is performed on all the dataset using as features the highly localized genes provided by the function **CIARA** 
3. **Merging of step 1 and 2**: Any cluster identified at step 2 of size  smaller than **max_number** constitutes an independent cluster.
4. **Identification of markers for small clusters**: It is possible to check if the rare cell types are well defined according to their markers with the function **white_black_markers**. See ?**white_black_markers** for more information.
5. **Identification of cell sub-types**: For each of the cluster identified at step 3, a Fisher test is performed to see if there is a statistically significant enrichment between the HVGs in the cluster and the highly localized genes (function **test_hvg**). For all the clusters where there is an enrichment, a sub cluster analysis starting from the original cluster is performed using as features  HVGs (function **cluster_analysis_sub**)

An example is the following:
```
# step 1
human_embryo_analysis <- cluster_analysis_integrate_rare(raw_counts = raw_counts, project_name = "Human_Embryo_data", resolution = 0.1, neighbors=5, max_dimension = 30)
human_embryo_cluster <- as.vector(human_embryo_analysis$seurat_clusters)
# or we can also use the original cluster annotation provided by Tyser et al,2020
original_cluster <- as.vector(meta_info$cluster_id)
human_embryo_cluster <- original_cluster
```

```
# step 2
human_embryo_analysis_ciara <- cluster_analysis_integrate_rare(raw_counts,"Human_Embryo_data",0.01,5,30,ciara_genes)
```
```
# step 3
final_cluster <- merge_cluster(original_cluster,human_embryo_analysis_ciara$seurat_clusters,max_number = 20)
final_cluster[grep("step_2",final_cluster)] <- "PGC"
```
```
# step 4
result_test <- test_hvg(raw_counts,final_cluster, ciara_genes, background, number_hvg = 100, min_p_value = 0.05)
result_test[[2]]
 "Endoderm"  "Hemogenic Endothelial Progenitors". "ExE Mesoderm"
# We need to do sub cluster in the above three clusters
raw_endoderm <- raw_counts[, human_embryo_cluster == "Endoderm"]
raw_emo <- raw_elmir[, human_embryo_cluster == "Hemogenic Endothelial Progenitors"]
raw_exe_meso=raw_elmir[, human_embryo_cluster == "ExE Mesoderm"]
combined_endoderm <- Cluster_analysis_sub(raw_endoderm, 0.2, 5, 30, "Endoderm")
combined_exe_meso <- Cluster_analysis_sub(raw_exe_meso, 0.5, 5, 30, "ExE Mesoderm")
combined_emo <- Cluster_analysis_sub(raw_emo, 0.6, 5, 30, "Hemogenic Endothelial Progenitors")


all_sub_cluster <- c(combined_endoderm$seurat_clusters, combined_emo$seurat_clusters, combined_exe_meso$seurat_clusters)
final_cluster_version_sub <- merge_cluster(final_cluster, all_sub_cluster)
```

```
plot_umap(coordinate_umap, final_cluster_version_sub)
```
<img src="https://github.com/ScialdoneLab/CIARA/blob/main/figures/entropy_cluster.png" width="700" height="500">
With the cluster analysis based on CIARA we are able to detect two clusters (Endoderm_2 and and Hemogenic Endothelial Progenitors_4 highlighted in the plot) that were not reported in the original paper.

For more exhaustive information about the functions offered by CIARA for the identification of rare populations of cells  see **Tutorials section** below and the help page of the single functions. (*?function_name*).


## Vignette

The following vignette is available and completely reproducible. It uses single cell RNA seq from human embryo at the gastrulation state from [Tyser *et al.*, 2021](https://www.nature.com/articles/s41586-021-04158-y). The raw count matrix was downloaded from  [http://human-gastrula.net].
An extremely rare population of primordial germ cells (PGCs-7 cells) is easily identified with entropy of mixing.
It can be accessed within R with:
```
utils::vignette("CIARA")
```

