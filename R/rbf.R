#' @export
radial_basis_fn <- function(matDR, knn){

  #The Only Differences Are Due To the Slight Difference in KNN algo

  #Median Distance
  median_distances <- rowMedians(knn$nn.idx)
  #n <- floor(ncol(knn$nn.idx)/2)
  #median_distances <- knn[[2]][,n + 1]

  #KNN Graph
  knn_graph <- knn$nn.graph
  knn_graph@x[knn_graph@x > 0] <- 1
  diag(knn_graph) <- 1

  #Create Symmetric Graph
  sym_graph <- (knn_graph + Matrix::t(knn_graph)) > 0
  sym_graph <- drop0(sym_graph)
  sym_graph <- as(sym_graph, "dgCMatrix")

  #Summarize
  df_sym_graph <- Matrix::summary(sym_graph)[, 1:2]

  #RBF
  numerator <- colSums((t(matDR[df_sym_graph$j,]) - t(matDR[df_sym_graph$i,]))^2)
  denominator <- median_distances[df_sym_graph$j] * median_distances[df_sym_graph$i]
  vals <- exp(-numerator / denominator)

  #M
  m <- sparseMatrix(
    i = df_sym_graph$i,
    j = df_sym_graph$j,
    x = vals
  )

  m %c*% Matrix::t(m)

}