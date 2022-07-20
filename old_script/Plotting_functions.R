#' plot_umap
#'
#' @param coordinate_umap Data frame with dimensionality reduction coordinates.
#' Number of rows must be equal to the number of cells
#' @inheritParams test_hvg
#' @return ggplot2 object.
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#' @seealso \url{https://CRAN.R-project.org/package=ggplot2}
#' @export plot_umap
plot_umap <- function(coordinate_umap, cluster){
  umap_plot <- ggplot2::ggplot(coordinate_umap, ggplot2::aes(coordinate_umap[, 1], coordinate_umap[, 2] )) +
    ggplot2::geom_point(ggplot2::aes(colour = as.factor(cluster))) + ggplot2::xlab("UMAP 1") + ggplot2::ylab("UMAP 2") + ggplot2::labs(col = "Cluster")
  return(list(umap_plot))
}



#' gg_color_hue
#' @noRd
gg_color_hue <- function(n) {
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}








#' plot_gene
#'
#' Cells are coloured according to the expression of \emph{gene_id} and plotted
#' according to \emph{coordinate_umap}.
#'
#' @param gene_id Character name of the gene.
#' @param title_name Character name.
#' @inheritParams get_background_full
#' @inheritParams plot_umap
#' @return ggplot2 object.
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#' @seealso \url{https://CRAN.R-project.org/package=ggplot2}
#'
#' @export plot_gene
plot_gene <- function(norm_matrix, coordinate_umap, gene_id, title_name){


  gene_expr <- as.vector(norm_matrix[gene_id, ])
  index_sort <- order(gene_expr)
  row.names(coordinate_umap) <- colnames(norm_matrix)
  coordinate_umap <- coordinate_umap[colnames(norm_matrix)[index_sort], ]
  gene_expr <- sort(gene_expr)
  umap_plot <- ggplot2::ggplot(coordinate_umap, ggplot2::aes(coordinate_umap[, 1], coordinate_umap[, 2] )) +
    ggplot2::geom_point(ggplot2::aes(colour = (gene_expr)))  + ggplot2::xlab("UMAP 1") + ggplot2::ylab("UMAP 2") + ggplot2::labs(col = "Log norm counts")+
    ggplot2::scale_colour_gradient(low = "white", high = "blue") +
    ggplot2::ggtitle(title_name)
  return((umap_plot))
}









#' plot_heatmap_marker
#'
#' @param marker_top First element returned by \emph{markers_cluster_seurat}
#' @param marker_all_cluster Second element returned by
#' \emph{markers_cluster_seurat}
#' @param condition Vector or length equal to the number of cells, specifying
#' the condition of the cells (i.e. batch, dataset of origin..)
#' @param text_size Size of the text in the heatmap plot.
#' @inheritParams test_hvg
#' @inheritParams get_background_full
#' @return Heatmap class object.
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#' @seealso \url{https://www.rdocumentation.org/packages/ComplexHeatmap/versions/1.10.2/topics/Heatmap}
#'
#' @export plot_heatmap_marker
plot_heatmap_marker <- function(marker_top, marker_all_cluster, cluster, condition, norm_matrix, text_size){
  if (!(requireNamespace("ComplexHeatmap", quietly = TRUE) & requireNamespace("circlize", quietly = TRUE))) {
    stop("Package ComplexHeatmap and circlize needed for this function to work. Please install them: BiocManager::install('ComplexHeatmap') and install.packages('circlize')")
  }
  cluster_unique <- unique(cluster)
  cluster_unique <- sort(cluster_unique, decreasing = F)

  index <- as.list(rep(0, length(unique(cluster))))
  for (i in 1:length(cluster_unique)){
    index[[i]] <- which(cluster==cluster_unique[i])
  }

  condition_unique <- as.list(rep(0, length(unique(cluster))))
  for (i in 1:length(cluster_unique)){
    condition_unique[[i]] <- condition[index[[i]]]
  }

  cell_all_day <- colnames(norm_matrix)
  cell_unique <- as.list(rep(0, length(unique(cluster))))
  for (i in 1:length(cluster_unique)){
    cell_unique[[i]] <- cell_all_day[index[[i]]]
  }

  cluster_unique_2 <- as.list(rep(0, length(unique(cluster))))
  for (i in 1:length(cluster_unique)){
    cluster_unique_2[[i]] <- cluster[index[[i]]]
  }

  cluster_finale <- unlist(cluster_unique_2)

  good_finale <- unlist(cell_unique)
  medium_finale <- unlist(condition_unique)

  norm_matrix_plot <- norm_matrix[marker_top, good_finale]

  row.names(norm_matrix_plot)<-marker_all_cluster

  color_cluster <- rep(0, length(unique(cluster)))

  for (i in 1:length(color_cluster) ){
    color_cluster[i] <- gg_color_hue(length(color_cluster))[i]
  }
  names(color_cluster) <- as.character(cluster_unique)


  color_condition <- rep(0, length(unique(condition)))

  for (i in 1:length(color_condition) ){
    color_condition[i] <- gg_color_hue(length(color_condition))[i]
  }
  names(color_condition) <- as.character(unique(condition))

  condition_factor <- factor(condition)

  names(color_condition) <- as.character(levels(condition_factor))



  data_heatmap <- data.frame(Cluster = cluster_finale
                          ,Condition = medium_finale)


  haCol1 <- ComplexHeatmap::HeatmapAnnotation(df = data_heatmap, col = list(Cluster = color_cluster, Condition = color_condition)

                             , show_legend = T)

  ht21 <- ComplexHeatmap::Heatmap(as.matrix(norm_matrix_plot)
                 ,cluster_rows = FALSE
                 , col = circlize::colorRamp2(c(0, round(max(norm_matrix_plot))), c("white", "red"))
                 , name = "log norm counts"
                 , cluster_columns = F
                 , top_annotation = haCol1
                 , row_names_gp = grid::gpar(fontsize = text_size
                 )
                 , show_column_names = F
                 , show_row_names = T)


  ComplexHeatmap::draw(ht21)



}






#' detect_localized_genes
#'
#' @param genes_sort vector of genes names sorted according to the p value given by the output of 
#' \emph{CIARA}.
#' @param top_number Number of \emph{genes_sort} to consider.
#' @inheritParams CIARA_gene
#' @inheritParams get_background_full
#' @return List with two elements. 
#' \item{list_intersect}{The first element is a list with length equal to the number of cells in \emph{norm_matrix}. 
#' Each entry contains the name of \emph{genes_sort} expressed in at least two neighboring cells.}
#' \item{rank_intersect}{The second element is a vector with length equal to the number of cells in \emph{norm_matrix}.Each entry contains the number 
#' of \emph{genes_sort} expressed in at least two neighboring cells.}

#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#'
#' @export detect_localized_genes

detect_localized_genes=function(knn_matrix, norm_matrix, genes_sort, top_number, threshold = 1 ){
  list_intersect <- rep(list(0), length(colnames(knn_matrix)))
  rank_intersect <- rep(0, length(colnames(knn_matrix)))
  for ( i in 1:length(colnames(knn_matrix))){
    knn_name=names(which(knn_matrix[i,] == 1))
    genes_select <- genes_sort[1:top_number]
    sum_knn <- apply(norm_matrix[genes_select, knn_name], 2, function(x){
      somma <- genes_select[which(x > threshold)]
      return(somma)
    })
    list_intersect[[i]] <- unique(as.vector(unlist(sum_knn))[duplicated(as.vector(unlist(sum_knn)))])
    rank_intersect[i] <- length(list_intersect[[i]])
    
  }
  return(list(list_intersect,rank_intersect))
}


#' names_localized_genes
#' 
#' @param list_intersect first element of the list given as output from \emph{detect_localized_genes}.
#' @param max_number maximum number of genes in \emph{list_intersect} to show for each cell.
#' @return A vector with length equal to \emph{list_intersect}. Each entry contains the top \emph{max_number} names of genes from \emph{list_intersect}.
#' @inheritParams detect_localized_genes
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#' 
#' 
#' @export names_localized_genes
names_localized_genes=function (list_intersect, genes_sort, max_number = 5){ 
  text <- rep(0,length(list_intersect))
  
  for (i in 1:length(list_intersect)) {
    if (length(list_intersect[[i]]) == 0){
      localized_genes <- "No specific genes"
      text[i] <- "No specific genes"
      
    }
    else {
      localized_genes <- genes_sort[which(genes_sort %in% list_intersect[[i]])]
      
      max_genes <- min(max_number, length(localized_genes))
      if (max_genes < max_number){
        warning("max_number bigger than length of localized_genes. Set max_number = length(localized_genes) ")
      }
      localized_genes <- localized_genes[1:max_genes]
      text[i] <- paste(localized_genes, collapse = "-")
      
      
    }
  }
  return(text)
}


#' plot_localized_genes_interactive
#'
#' It shows in an interactive plot which are the highly localized genes in each
#' cell. It is based on plotly library
#' @param rank_intersect second element of the list given as output from \emph{detect_localized_genes}.
#' @param ramp_list vector of length equal to n_rows in coordinate_umap.Each cell
#' will be coloured following a gradient according to the corresponding value
#' of this vector.
#' @param text Character vector specifying the highly localized genes in each
#' cell. It is the output from \emph{names_localized_genes}.
#' @param ramp_list \emph{colors} parameter for the funnction \emph{plot_ly} from R package \emph{plotly}.
#' @param min_x Set the min limit on the x axis.
#' @param max_x Set the max limit on the x axis.
#' @param min_y Set the min limit on the y axis.
#' @param max_y Set the min limit on the y axis.
#' @inheritParams plot_umap
#' @inheritParams get_background_full
#' @return plotly object given by \emph{plot_ly function} (from library \emph{plotly}).
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#' @seealso \url{https://plotly.com/r/}
#' @export plot_localized_genes_interactive

plot_localized_genes_interactive = function (coordinate_umap, norm_matrix,rank_intersect, text,ramp_list, min_x = NULL, max_x = NULL, 
                                min_y = NULL, max_y = NULL)
{
  
  
  index_sort <- order(rank_intersect)
  text=text[index_sort]
  row.names(coordinate_umap) <- colnames(norm_matrix)
  coordinate_umap <- coordinate_umap[colnames(norm_matrix)[index_sort], 
                                     ]
  rank_intersect <- sort(rank_intersect)
  
  
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package plotly needed for interactive == TRUE. Please install it: install.packages('plotly')")
  }
  colnames(coordinate_umap) <- c("UMAP_1", "UMAP_2")
  fig <- plotly::plot_ly(data = coordinate_umap, x = ~UMAP_1, 
                         y = ~UMAP_2, color ~ rank_intersect, type = "scatter", mode = "markers", 
                         marker = list(size = 5, width = 2,line = list(color = "black",
                                                                       width = 0.5)), text = ~text, hoverinfo = "text",colors=ramp_list)
  
  
  
  if (!is.null(min_x)) {
    fig <- fig %>% plotly::layout(xaxis = list(range = c(min_x, 
                                                         max_x)), yaxis = list(range = c(min_y, max_y))) %>% 
      plotly::layout(plot_bgcolor='rgb(254, 247, 234)') %>% 
      plotly::layout(paper_bgcolor='rgb(254, 247, 234)')
  }
  return(fig)
}


#' plot_localized_genes
#'
#'
#' @param colors  \emph{colours} parameter for the function scale_colour_gradientn from R package \emph{ggplot2}.
#' the sum in each cell.
#' @param name_title Character value.
#' @inheritParams plot_umap
#' @inheritParams get_background_full
#' @inheritParams plot_localized_genes_interactive
#' @return ggplot2 object.
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#' @seealso \url{ https://CRAN.R-project.org/package=ggplot2}
#'
#' @export plot_localized_genes

plot_localized_genes = function (coordinate_umap, norm_matrix, rank_intersect, name_title, colors) {
  index_sort <- order(rank_intersect)
  row.names(coordinate_umap) <- colnames(norm_matrix)
  coordinate_umap <- coordinate_umap[colnames(norm_matrix)[index_sort],]
  rank_intersect <- sort(rank_intersect)
  umap_plot <- ggplot2::ggplot(coordinate_umap, ggplot2::aes(coordinate_umap[, 1], coordinate_umap[, 2])) + 
    ggplot2::geom_point(ggplot2::aes(colour = rank_intersect),size = 1) +
    ggplot2::xlab("UMAP 1") + ggplot2::ylab("UMAP 2") + ggplot2::labs(col = "Number of localized genes") + 
    ggplot2::scale_colour_gradientn(colours = colors)+
    ggplot2::ggtitle(name_title)
  return((umap_plot))
}



#' plot_balloon_marker
#'
#' @param marker_complete Third element of the output list as
#' returned by the function \emph{markers_cluster_seurat}
#' @param max_number Integer. Maximum number of markers for each cluster for
#' which we want to plot the expression.
#' @param max_size Integer. Size of the dots to be plotted.
#' @inheritParams test_hvg
#' @inheritParams plot_heatmap_marker
#' @inheritParams get_background_full
#' @return ggplot2 object showing balloon plot.
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#' @export plot_balloon_marker
plot_balloon_marker <- function (norm_matrix, cluster, marker_complete,
                          max_number, max_size = 5, text_size = 7) {
  livelli <- levels(factor(cluster))
  all_markers_plot <- rep(list(0), length(livelli))
  all_markers <- rep(list(0), length(livelli))
  for (i in 1:length(livelli)) {
    genes_4_0 <- marker_complete[(names(marker_complete) == livelli[i])]
    if (length(genes_4_0) > 0) {
      genes_4_0 <- genes_4_0[1:max_number]
      all_markers[[i]] <- genes_4_0
      genes_4_0_plot <- paste(livelli[i], genes_4_0, sep = "_")
      all_markers_plot[[i]] <- genes_4_0_plot
    }
  }
  all_markers <- lapply(all_markers, function(x) {
    x[x != 0]
  })
  all_markers_plot <- lapply(all_markers_plot, function(x) {
    x[x != 0]
  })
  all_markers_values <- rep(list(0), length(livelli))
  for (i in 1:length(livelli)) {
    gene_values_0 <- apply(norm_matrix[unlist(all_markers, use.names = FALSE), cluster == livelli[i]], 1, mean)
    all_markers_values[[i]] <- gene_values_0
  }
  values <- rep(list(0), length(livelli))
  for (i in 1:length(livelli)) {
    valore <- apply(norm_matrix[unlist(all_markers, use.names = FALSE), cluster == livelli[i]], 1, function(x) {
      sum(x > 1)/length(x)
      })
    values[[i]] <- valore
  }
  cluster_label_all <- rep(list(0), length(livelli))
  for (i in 1:length(livelli)) {
    cluster_label_all[[i]] <- rep(livelli[i], length(unlist(all_markers)))
  }
  all_markers_plot <- unlist(all_markers_plot)
  cluster_label_all <- unlist(cluster_label_all)
  values <- unlist(values)
  all_markers_values <- unlist(all_markers_values)
  balloon_melted <- data.frame(all_markers_plot, cluster_label_all, values)
  colnames(balloon_melted) <- c("all_markers_plot", "cluster_label_all", "values")
  p <- ggplot2::ggplot(balloon_melted, ggplot2::aes(x = cluster_label_all, y = all_markers_plot))
  p + ggplot2::geom_point(ggplot2::aes(size = values, colour = as.numeric(all_markers_values))) + ggplot2::theme(panel.background = ggplot2::element_blank(), panel.border = ggplot2::element_rect(colour = "blue", fill = NA,  size = 3)) +
    ggplot2::scale_size_area(max_size = max_size) + ggplot2::scale_colour_gradient(low = "white", high = "midnight blue", na.value = NA) +
    ggplot2::labs(colour = "Mean log2 norm expression", size = "Value (fraction cells above 1 log2 norm counts)") +
    ggplot2::xlab("Cluster") + ggplot2::ylab("Markers")  + ggplot2::theme(text = ggplot2::element_text(size = text_size), axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1), axis.title.x = ggplot2::element_blank())
}



