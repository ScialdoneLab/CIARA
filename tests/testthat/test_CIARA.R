# Example with a norm count matrix of 10 genes and 10 samples
# All the genes have the same expression (10) value in all the samples.
norm_test = matrix(10, ncol = 10, nrow = 10)
row.names(norm_test) = paste0("Gene", seq(1:10))
colnames(norm_test) = paste0("Sample", seq(1:10))

test_that("get_background_full returns error when no genes pass the filtering steps", {
  expect_error(get_background_full(norm_test, threshold = 11, n_cells_low = 3, n_cells_high = 20), "There are not genes that pass the filtering steps")
})

test_that("get_background_full returns the full set of genes when the threshold is below the minimum value in norm_test  ", {
  expect_identical(get_background_full(norm_test, threshold = 1, n_cells_low = 3, n_cells_high = 20), row.names(norm_test))
})


# Built a knn_test matrix where all the 10 samples are neighbours of each other. The row.names of knn_test are different from
# colnames of norm_test
knn_test = matrix(1, ncol = 10, nrow = 10)
colnames(knn_test) = paste0("Sample_different_name", seq(1:10))
row.names(knn_test) = paste0("Sample_different_name", seq(1:10))
gene_expression= as.vector(norm_test[1,])

test_that("CIARA_gene returns error when the names of row.names in knn_matrix are not equal to colnames in norm_matrix", {
  expect_error(CIARA_gene(norm_test, knn_test, gene_expression, p_value = 0.001, odds_ratio = 2, local_region = 1, approximation = FALSE), "row.names in knn_matrix are not equal to colnames in norm_matrix")
})

# Built a knn_test matrix where all the 10 samples are neighbours of each other
knn_test = matrix(1, ncol = 10, nrow = 10)
colnames(knn_test) = paste0("Sample", seq(1:10))
row.names(knn_test) = paste0("Sample", seq(1:10))

gene_expression= as.vector(norm_test[1,])

test_that("CIARA_gene returns warning when the names of row.names in knn_matrix are not equal to colnames in norm_matrix", {
  expect_warning(CIARA_gene(norm_test, knn_test, gene_expression, p_value = 0.001, odds_ratio = 2, local_region = 1, approximation = FALSE), "There are not genes enriched in 1 or more local regions")
  expect_message(CIARA_gene(norm_test, knn_test, gene_expression, p_value = 0.001, odds_ratio = 2, local_region = 1, approximation = FALSE), "approximation == FALSE")
  expect_message(CIARA_gene(norm_test, knn_test, gene_expression, p_value = 0.001, odds_ratio = 2, local_region = 1, approximation = TRUE), "approximation == TRUE")
})



# Built a background genes vector with one gene not present in norm_test matrix

background = paste0("Gene", seq(1:11))
test_that("CIARA returns error when some background genes are not present in norm_test matrix", {
  expect_error(CIARA(norm_test, knn_test, background, cores_number = 1,p_value = 0.001,odds_ratio = 2,local_region = 1, approximation = FALSE), "Some background genes are not present in norm matrix")
})


background = paste0("Gene", seq(1:10))
output_ciara_test = matrix(1, ncol = 1, nrow = 10)
row.names(output_ciara_test) = background
test_that("CIARA returns the expected a dataframe with all 1 when there are not genes enriched in 1 or more local regions", {
  expect_identical(CIARA(norm_test, knn_test, background, cores_number = 1,p_value = 0.001,odds_ratio = 2,local_region = 1, approximation = FALSE), output_ciara_test)
})





