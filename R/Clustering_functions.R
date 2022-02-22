#' cluster_analysis_integrate_rare
#'
#' @param raw_counts Raw count matrix (n_genes X n_cells).
#' @param project_name Character name of the Seurat project.
#' @param resolution Numeric value specifying the parameter \emph{resolution}
#' used in the Seurat function \emph{FindClusters}.
#' @param neighbors Numeric value specifying the parameter \emph{k.param} in
#' the Seurat function \emph{FindNeighbors}
#' @param max_dimension Numeric value specifying the maximum number of the PCA
#' dimensions used in the parameter \emph{dims} for the Seurat function
#' \emph{FindNeighbors}
#' @param feature_genes vector of features specifying the argument
#' \emph{features} in the Seurat function \emph{RunPCA}.
#' @return Seurat object including raw and normalized counts matrices, UMAP coordinates and cluster result.
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#' @seealso
#' \url{https://www.rdocumentation.org/packages/Seurat/versions/4.0.1/topics/FindClusters}
#' \url{https://www.rdocumentation.org/packages/Seurat/versions/4.0.1/topics/FindNeighbors}
#' \url{https://www.rdocumentation.org/packages/Seurat/versions/4.0.1/topics/RunPCA}
#' @export cluster_analysis_integrate_rare
cluster_analysis_integrate_rare <- function(raw_counts, project_name, resolution, neighbors, max_dimension, feature_genes = NULL) {
  if (!(requireNamespace("Seurat", quietly = TRUE))) {
    stop("Package Seurat needed for this function to work. Please install it: install.packages('Seurat')")
  }
  nFeature_RNA <- NULL
  seurat_object <- Seurat::CreateSeuratObject(counts = raw_counts, project = project_name)
  seurat_object <- subset( seurat_object, subset = nFeature_RNA > 0)
  seurat_object <- Seurat::NormalizeData( seurat_object, verbose = FALSE)
  seurat_object <- Seurat::FindVariableFeatures( seurat_object, selection.method = "vst", nfeatures = 2000)
  seurat_object <- Seurat::ScaleData( seurat_object, verbose = FALSE)
  seurat_object <- Seurat::RunPCA( seurat_object, npcs = max_dimension, verbose = FALSE, features = feature_genes)
  seurat_object <- Seurat::RunUMAP( seurat_object, reduction = "pca", dims = 1:20, set.seed = 42)
  seurat_object <- Seurat::FindNeighbors( seurat_object, reduction = "pca", dims = 1:max_dimension, k.param = neighbors)
  seurat_object <- Seurat::FindClusters( seurat_object, resolution = resolution)
  return(seurat_object)}





#' find_resolution
#'
#' @param seurat_object Seurat object as returned by
#' \emph{cluster_analysis_integrate_rare}
#' @param resolution_vector vector with all values of resolution for which the
#' Seurat function \emph{FindClusters} is run
#' @return Clustree object showing the connection between clusters obtained at different level of resolution as specified in \emph{resolution_vector}.
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#' @seealso \url{https://CRAN.R-project.org/package=clustree}
#'
#' @export find_resolution
find_resolution <- function(seurat_object, resolution_vector) {
  if (!(requireNamespace("Seurat", quietly = TRUE) & requireNamespace("clustree", quietly = TRUE))) {
    stop("Package Seurat and clustree needed for this function to work. Please install them: install.packages('Seurat') and install.packages('clustree')")
  }
  for (i in resolution_vector) {
    seurat_object <- Seurat::FindClusters(seurat_object, resolution = i)
  }
  clustree::clustree(seurat_object)
}

#' @importFrom ggraph guide_edge_colourbar
#' @export
ggraph::guide_edge_colourbar



#' build_fisher_test
#' @noRd
build_fisher_test <- function(a, b, c) {
  matrix_test <- matrix(c(sum(c[c%in%a]%in%c[c%in%b]), sum(c[c%in%a]%in%c[!c%in%b]), sum(c[!c%in%a]%in%c[c%in%b]), sum(c[!c%in%a]%in%c[!c%in%b])), nrow = 2, ncol = 2, byrow = T)
  fisher.test(matrix_test)
}






#' de_seurat_cluster
#' @noRd
de_seurat_cluster <- function(seurat_object, cluster, names_cell, max_p_value) {
  level <- levels(as.factor(cluster))
  final_markers <- vector("list", length(level))
  for (i in 1:length(level)) {
    markers <- Seurat::FindMarkers(seurat_object, ident.1 = names_cell[cluster == level[i]], ident.2 = names_cell[cluster!=level[i]], only.pos = T)
    markers_new <- markers[markers$p_val_adj <= max_p_value, ]
    markers_new <- markers_new[order(markers_new$p_val_adj), ]
    markers_final <- row.names(markers_new)
    final_markers[[i]] <- markers_final
  }
  return(final_markers)
}











#' merge_cluster
#'
#' @param old_cluster original cluster assignment that need to be updated
#' @param new_cluster new cluster assignment that need to be integrated with
#' \emph{old_cluster}.
#' @param max_number Threshold in size for clusters in \emph{new_cluster}. Only
#' cluster with number of cells smaller than \emph{max_number} will be
#' integrated in \emph{old cluster}. If \emph{max_number} is NULL, then all the clusters in \emph{new_cluster} are integrated in \emph{old cluster}.
#' @return Numeric vector of length equal to \emph{old_cluster} showing the merged cluster assignment between \emph{old cluster} and \emph{new_cluster}.
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#'
#' @export merge_cluster
merge_cluster <- function(old_cluster, new_cluster, max_number = NULL) {

  if(is.null(max_number)) {
    cluster_final <- old_cluster
    cluster_final[match(names(new_cluster),names(old_cluster))] <- new_cluster
    return(cluster_final)
  }
  else{
    index_small <- which(table(new_cluster)<max_number)
    cluster_small <- names(table(new_cluster))[index_small]
    if (length(cluster_small) > 0) {
      message(paste0("Clusters with number of cells below ", max_number, " in the new partition: ", length(cluster_small)))
      cluster_final <- old_cluster
    for (i in 1:length(cluster_small)) {
      cluster_final[as.vector(new_cluster) == cluster_small[i]] <- paste(cluster_small[i], "step_2", sep = "-")
    }
    return(cluster_final)
    }
    else {
      warning(paste0("There are not clusters with number of cells below ", max_number))
      cluster_final <- old_cluster
      return(cluster_final)
    }
  }
}





#' test_hvg
#'
#' For each cluster in \emph{cluster}, HVGs are defined with
#' Seurat function \emph{FindVariableFeatures}. A Fisher test is performed to
#' see if there is a statistically significant enrichment between the top
#' \emph{number_hvg} and the \emph{localized_genes}
#'
#' @param cluster Vector of length equal to the number of cells, with cluster
#' assignment.
#' @param localized_genes Character vector with localized genes detected by CIARA.
#' @param background Character vector with all the genes names to use as
#' background for the Fisher test.
#' @param number_hvg Integer value. Number of top HVGs provided by the Seurat
#' function \emph{FindVariableFeatures}.
#' @param min_p_value Threshold on p values provided by Fisher test.
#' @inheritParams cluster_analysis_integrate_rare
#' @return A list with two elements. \item{first element}{The first one is a
#' list with length equal to the number of clusters. Each entry is list of
#' three elements. The first two elements contain the p value and the odds
#' ration given by the Fisher test The third is a vector with genes names that
#' are present both in \emph{localized_genes} and in top \emph{number_hvg} HVGs
#' .} \item{second element}{a character vector with the name of the cluster that
#' have a p value smaller than \emph{min_p_value}.}
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#' @seealso
#' \url{https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/fisher.test}
#' @export test_hvg
test_hvg <- function (raw_counts, cluster, localized_genes, background, number_hvg, min_p_value)
{
  if (!(requireNamespace("Seurat", quietly = TRUE))) {
    stop("Package Seurat needed for this function to work. Please install it: install.packages('Seurat')")
  }
  nFeature_RNA <- NULL
  levels_cluster <- levels(as.factor(cluster))
  final_p_value <- rep(list(0), length(levels))
  for (i in 1:length(levels_cluster)) {
    raw_counts_small <- raw_counts[, cluster == levels_cluster[i]]
    seurat_object <- Seurat::CreateSeuratObject(counts = raw_counts_small, project = "sub_cluster")
    seurat_object <- subset(seurat_object, subset = nFeature_RNA > 0)
    seurat_object <- Seurat::NormalizeData(seurat_object, verbose = FALSE)
    seurat_object <- Seurat::FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = number_hvg)
    hvg_genes <- seurat_object@assays$RNA@var.features
    fisher_output <- build_fisher_test(localized_genes, hvg_genes, background)
    final_p_value[[i]] <- list(fisher_output[[1]], fisher_output[[3]], localized_genes[localized_genes %in% hvg_genes])
  }
  names(final_p_value) <- levels_cluster
  select_sub <- lapply(final_p_value, function(x) {if(x[[1]]<min_p_value & x[[2]] > 1) {
    return(TRUE)}
    else{return(FALSE)}
  })
  final_sub <- names(select_sub)[select_sub == TRUE]
  return(list(final_p_value, final_sub))
}
















#' cluster_analysis_sub
#'
#' @param name_cluster Character.Name of the original cluster for which the sub
#' clustering is done.
#' @inheritParams cluster_analysis_integrate_rare
#' @return Seurat object including raw and normalized counts matrices and cluster result.
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#' @seealso
#' \url{https://www.rdocumentation.org/packages/Seurat/versions/4.0.1/topics/RunPCA}
#' \url{https://www.rdocumentation.org/packages/Seurat/versions/4.0.1/topics/FindVariableFeatures}
#' @export cluster_analysis_sub
cluster_analysis_sub <- function (raw_counts, resolution, neighbors, max_dimension, name_cluster) {
  if (!(requireNamespace("Seurat", quietly = TRUE))) {
    stop("Package Seurat needed for this function to work. Please install it: install.packages('Seurat')")
  }
  nFeature_RNA <- NULL
  seurat_object <- Seurat::CreateSeuratObject(counts = raw_counts, project = "sub_cluster")
  seurat_object <- subset(seurat_object, subset = nFeature_RNA > 0)
  seurat_object <- Seurat::NormalizeData(seurat_object, verbose = FALSE)
  seurat_object <- Seurat::FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)
  seurat_object <- Seurat::ScaleData(seurat_object, verbose = FALSE)
  hvg_genes <- seurat_object@assays$RNA@var.features
  genes_cluster <- hvg_genes
  seurat_object <- Seurat::RunPCA(seurat_object, npcs = max_dimension, verbose = FALSE, features = genes_cluster)
  seurat_object <- Seurat::FindNeighbors(seurat_object, reduction = "pca", dims = 1:max_dimension, k.param = neighbors)
  seurat_object <- Seurat::FindClusters(seurat_object, resolution = resolution)
  new_cluster <- as.vector(seurat_object$seurat_clusters)
  for (i in 1:length(levels(seurat_object$seurat_clusters))) {
    new_cluster[as.vector(seurat_object$seurat_clusters) == levels(seurat_object$seurat_clusters)[i]] <- paste(name_cluster,
                                                                                                               new_cluster[as.vector(seurat_object$seurat_clusters) == levels(seurat_object$seurat_clusters)[i]], sep = "_")
  }
  names(new_cluster) <- names(seurat_object$seurat_clusters)
  seurat_object$seurat_clusters <- new_cluster
  return(seurat_object)
}







#' markers_cluster_seurat
#'
#' The Seurat function \emph{FindMarkers} is used to identify general marker
#' for each cluster (specific cluster vs all other cluster). This list of
#' markers is then filtered keeping only the genes that appear as markers in a
#' unique cluster.
#'
#' @param seurat_object Seurat object as returned by
#' \emph{cluster_analysis_sub} or by \emph{cluster_analysis_integrate_rare}.
#' @param cell_names Vector of length equal to the number of cells, with cell
#' names.
#' @param number_top Integer. Number of top marker genes to keep for each
#' cluster.
#' @inheritParams test_hvg
#' @return List of three elements. The first is a vector with \emph{number_top}
#' marker genes for each cluster. The second is a vector with \emph{number_top}
#' marker genes and corresponding cluster. The third element is a vector with
#' all marker genes for each cluster.
#'
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#' @seealso \url{https://www.rdocumentation.org/packages/Seurat/versions/4.0.1/topics/FindMarkers}
#'
#' @export markers_cluster_seurat
markers_cluster_seurat <- function(seurat_object, cluster, cell_names, number_top) {
  if (!(requireNamespace("Seurat", quietly = TRUE))) {
    stop("Package Seurat needed for this function to work. Please install it: install.packages('Seurat')")
  }
  level <- levels(as.factor(cluster))
  marker_cluster <- as.list(rep(0, length(level)))
  marker_all <- de_seurat_cluster(seurat_object, cluster, cell_names, 0.05)

  for (i in 1:length(level)) {
    marker_cluster[[i]] <- marker_all[[i]]
    message(paste0("Cluster ", level[i]))
  }
  marker_all <- unlist(marker_cluster)
  marker_all <- marker_all[Biobase::isUnique(marker_all)]

  for (i in 1:length(level)) {
    marker_cluster[[i]] <- marker_cluster[[i]][marker_cluster[[i]]%in%marker_all]
  }
  marker_to_see <- as.list(rep(0, length(level)))
  marker_complete <- as.list(rep(0, length(level)))
  for (i in 1:length(level)) {
    if (length(marker_cluster[[i]]) >= number_top) {
      marker_to_see[[i]] <- marker_cluster[[i]][1:number_top]
      names(marker_to_see[[i]]) <- rep(level[i], length(marker_to_see[[i]]))
      marker_complete[[i]] <- marker_cluster[[i]]
      names(marker_complete[[i]]) <- rep(level[i], length(marker_complete[[i]]))
    }
    if ((length(marker_cluster[[i]]) > 0)&(length(marker_cluster[[i]])<number_top)) {
      marker_to_see[[i]] <- marker_cluster[[i]][1:length(marker_cluster[[i]])]
      names(marker_to_see[[i]]) <- rep(level[i], length(marker_to_see[[i]]))
      marker_complete[[i]] <- marker_cluster[[i]]
      names(marker_complete[[i]]) <- rep(level[i], length(marker_complete[[i]]))
    }
    if (length(marker_cluster[[i]]) == 0) {
      marker_to_see[[i]] <- NULL
      marker_complete[[i]] <- NULL
    }
  }



  marker_top <- unlist(marker_to_see)
  marker_top <- marker_top[marker_top!=0]
  marker_complete <- unlist(marker_complete)
  marker_complete <- marker_complete[marker_complete!=0]
  marker_all_cluster_sing <- as.list(rep(0, length(level)))


  for (i in 1:length(level)) {
    if(sum(names(marker_top) == level[i]) > 0) {
      marker_all_cluster_sing[[i]] <- paste(marker_top[names(marker_top) == level[i]], level[i], sep="_")
    }
    else{}
  }
  marker_all_cluster <- unlist(marker_all_cluster_sing)
  marker_all_cluster <- marker_all_cluster[marker_all_cluster!=0]

  return(list(marker_top, marker_all_cluster, marker_complete))}







#' white_black_markers
#'
#' A white-marker is a gene whose median expression across cells belong to
#' \emph{single_cluster} is greater than \emph{threshold} and in all the other
#' clusters is equal to zero.
#'
#' @param single_cluster Character. Label of one specify cluster
#' @param marker_list Third element of the output list as returned by the
#' function \emph{markers_cluster_seurat}
#' @param threshold Numeric. The median of the genes across cells belong to
#' \emph{single_cluster} has to be greater than \emph{threshold} in order to be
#' consider as a white-black marker for \emph{single_cluster}
#' @inheritParams test_hvg
#' @inheritParams cluster_analysis_integrate_rare
#' @inheritParams plot_genes_sum
#' @return Logical vector of length equal to \emph{marker_list}, with
#' TRUE/FALSE if the gene is/is not a white-black marker for
#' \emph{single_cluster}.
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#'
#' @export white_black_markers
white_black_markers <- function(cluster, single_cluster, norm_counts, marker_list, threshold = 0) {
  white_black <- apply(norm_counts[marker_list[names(marker_list) == single_cluster], ], 1, function(x) {
  mean_one <- median(x[cluster == single_cluster])
  mean_different <- rep(0, length(levels(factor(cluster[cluster!=single_cluster]))))
  for ( i in 1:length(levels(factor(cluster[cluster!=single_cluster])))) {
    mean_different[i] <- median(x[cluster == levels(factor(cluster[cluster!=single_cluster]))[i]])}
  if (mean_one > threshold&sum(mean_different) == 0) {
    return(TRUE)
  }
  else{
    return(FALSE)
    }
  })
  return(white_black)
}

