#' get_background_full
#'
#' @param norm_matrix Norm count matrix (n_genes X n_cells).
#' @param threshold threshold in expression for a given gene
#' @param n_cells_low minimum number of cells where a gene is expressed at a
#' level above threshold
#' @param n_cells_high maximum number of cells where a gene is expressed at a
#' level above threshold
#' @return Character vector with all genes expressed at a level higher than
#' \emph{threshold} in a number of cells between \emph{n_cells} and
#' \emph{n_cells_high}.
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#'
#'
#'
#' @export get_background_full
#'
#' @importFrom magrittr %>%
#' @importFrom parallel mclapply
#' @importFrom grDevices hcl
#' @importFrom stats fisher.test
#' @importFrom stats median
get_background_full <- function(norm_matrix, threshold = 1, n_cells_low = 3, n_cells_high = 20) {

  genes_filter <- apply(norm_matrix, 1, function(x) {
    x <- x[x > threshold]
    if (length(x) >= n_cells_low) {
      return(TRUE)
      }
    else{
      return(FALSE)
      }
  })

  genes_filter_2 <- apply(norm_matrix, 1, function(x) {
    x <- x[x > threshold]
    if (length(x) > n_cells_high) {
      return(FALSE)
      }
    else{
      return(TRUE)
      }
  })
  genes_important <- row.names(norm_matrix)[genes_filter]
  genes_important_2 <- row.names(norm_matrix)[genes_filter_2]
  genes_important_final <- intersect(genes_important, genes_important_2)
  if (length(genes_important_final) == 0) {
    stop("There are not genes that pass the filtering steps")
  }
  return(genes_important_final)
}










#' CIARA_gene
#'
#' The gene expression is binarized (1/0) if the value in a given cell is
#' above/below the median. Each of cell with its first K nearest neighbors
#' defined a local region. If there are at least \emph{local_region} enriched
#' in 1 according to \emph{fisher.test}, then the gene is defined as highly
#' localized and a final p value is assigned to it. The final p value is the
#' minimum of the p values from all the enriched local regions. If there are no
#' enriched local regions, then the p value by default is set to 1
#'
#' @param knn_matrix K-nearest neighbors matrix (n_cells X n_cells).
#' @param gene_expression numeric vector with the gene expression (length equal
#' to n_cells). The gene expression is binarized (equal to 0/1 in the cells
#' where the value is below/above the median)
#' @param p_value p value returned by the function \emph{fisher.test} with
#' parameter alternative = "g"
#' @param local_region Integer. Minimum number of local regions (cell with its
#' knn neighbours) where the binarized gene expression is enriched in 1.
#' @param approximation Logical.For a given gene, the fisher test is run in the
#' local regions of only the cells where the binarized gene expression is 1.
#' @inheritParams get_background_full
#' @return List with one element corresponding to the p value of the gene.
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#' @seealso \url{https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/fisher.test}
#' @export CIARA_gene
CIARA_gene <- function(norm_matrix, knn_matrix, gene_expression, p_value = 0.001, local_region = 1, approximation) {




  median_genes <- median(gene_expression)
  binary_expression <- rep(0, length(gene_expression))
  binary_expression[gene_expression > median_genes] <- 1
  binary_expression[gene_expression <= median_genes] <- 0

  if (!(is.logical(approximation))) {
    stop ("approximation must be a logical value. Try with approximation = TRUE or approximation = FALSE")

  }
  if (approximation == FALSE) {
    sub_feature <- colnames(norm_matrix)
    message("approximation == FALSE")
    }
  else {
    sub_feature <- colnames(norm_matrix)[which(binary_expression == 1)]
    message("approximation == TRUE")
  }
  if (!all(colnames(norm_matrix)%in%row.names(knn_matrix))) {
    stop ("row.names in knn_matrix are not equal to colnames in norm_matrix")

  }

  knn_matrix_big <- knn_matrix[sub_feature, ]



  perform_fisher <- function(j) {

    nn_cell <- colnames(norm_matrix)[which(knn_matrix_big[sub_feature[j], ] > 0)]
    nn_gene_expression <- binary_expression[colnames(norm_matrix)%in%nn_cell]
    if(sum(nn_gene_expression) != 0) {
      input_fisher <- matrix(c(sum(nn_gene_expression == 1), sum(binary_expression == 1)-sum(nn_gene_expression == 1), sum(nn_gene_expression == 0), sum(binary_expression == 0)-sum(nn_gene_expression == 0)), nrow = 2, ncol = 2, byrow = T)
      test_output <- fisher.test(input_fisher, alternative = "g")
      return(test_output$p.value)
    }
    else{
      return(1)
    }
  }

  p_value_gene <- unlist(lapply(seq_len(length(sub_feature)), perform_fisher))

  if(sum(p_value_gene < p_value) >= local_region) {
    p_value_final <- min(p_value_gene)
    }
  else {
    p_value_final <- 1
  }
  if (all(p_value_final == 1)){
    warning(paste0("There are not genes enriched in ", local_region, " or more local regions"))
  }

  return((p_value_final))
}




#' CIARA
#'
#' It selects highly localized genes as specified in \emph{CIARA_gene},
#' starting from genes in \emph{background}
#'
#' @param background Vector of genes for which the function \emph{CIARA_gene}
#' is run.
#' @param cores_number Integer.Number of cores to use.
#' @inheritParams CIARA_gene
#' @inheritParams get_background_full
#' @return Dataframe with n_rows equal to the length of
#' \emph{background} . Each row is the output from \emph{CIARA_gene}.
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#'
#'
#'
#' @export CIARA
CIARA <- function(norm_matrix, knn_matrix, background, cores_number = 1, p_value = 0.001, local_region = 1, approximation) {

  if (!(is.logical(approximation))) {
    stop ("approximation must be a logical value. Try with approximation = TRUE or approximation = FALSE")

  }

  run_loop_genes = function(i) {
    if (!all(background %in% row.names(norm_matrix))) {
      stop("Some background genes are not present in norm matrix")
    }
    gene_expression <- as.vector(norm_matrix[background[i],])
    message(paste0("Running CIARA on gene:", background[i]))
    return(CIARA_gene(norm_matrix, knn_matrix, gene_expression, p_value, local_region, approximation))
  }
  result_final <- do.call(rbind, mclapply(seq_len(length(background)),run_loop_genes , mc.cores = cores_number))
  row.names(result_final) <- background
  return(result_final)
}




















