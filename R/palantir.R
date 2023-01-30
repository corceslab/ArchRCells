#' @export
diffusion_maps <- function(
  matDR = NULL, 
  knn = NULL,
  n_components = 10
  ){

  #N Cells
  N <- nrow(knn$nn.idx)

  #To DataFrame
  dfKNN <- data.frame(Matrix::summary(knn$nn.graph))

  #Adaptive K
  adaptive_k = floor(ncol(knn$nn.idx) / 3)

  #Adaptive SD
  adaptive_std <- knn$nn.dist[,adaptive_k]
  # adaptive_std <- dfKNN %>% 
  #   {split(.$x, .$i)} %>%
  #     lapply(., function(x){
  #       sort(x)[adaptive_k]
  #     }) %>% unlist(use.names=FALSE)

  #To Symmetric Sparse Matrix
  W <- Matrix::sparseMatrix(
    i = dfKNN$i,
    j = dfKNN$j,
    x = exp(-dfKNN$x / adaptive_std[dfKNN$i]),
    dims = c(N, N)
  )

  #Diffusion Components
  kernel <- W + Matrix::t(W)

  # Markov
  D <- Matrix::colSums(kernel)

  #Norm
  D[D != 0] = 1 / D[D != 0]

  #CSR
  T <- Matrix::sparseMatrix(
    i = seq_len(N),
    j = seq_len(N),
    x = D,
    dims = c(N, N)
  )
  TK <- T %*% kernel

  #Eigen Value Decomp
  e <- RSpectra::eigs(TK, k = n_components)
  eD <- as.numeric(e$values)
  eV <- apply(e$vectors, 2, as.numeric)
  idx <- order(eD,decreasing=TRUE)
  eD <- eD[idx]
  eV <- eV[,idx]

  #Normalize
  eV <- apply(eV, 2, function(x){
    x / sqrt(sum(x^2))
  })

  #Create Object
  out <- SimpleList(
    T = TK,
    EigenVectors = eV,
    EigenValues = eD,
    Kernel = kernel  
  )

  out

}

# diffusion_maps <- function(
#   matDR = NULL, 
#   knn = NULL,
#   n_components = 10
#   ){

#   #N Cells
#   N <- nrow(knn$nn.idx)

#   #To DataFrame
#   dfKNN <- data.frame(Matrix::summary(knn$nn.graph))

#   #Adaptive K
#   adaptive_k = floor(ncol(knn$nn.idx) / 3)

#   #Adaptive SD
#   adaptive_std <- dfKNN %>% 
#     {split(.$x, .$i)} %>%
#       lapply(., function(x){
#         sort(x)[adaptive_k]
#       }) %>% unlist(use.names=FALSE)

#   #To Symmetric Sparse Matrix
#   W <- Matrix::sparseMatrix(
#     i = dfKNN$i,
#     j = dfKNN$j,
#     x = exp(-dfKNN$x / adaptive_std[dfKNN$i]),
#     dims = c(N, N)
#   )

#   #Diffusion Components
#   kernel = W + Matrix::t(W)

#   # Markov
#   D <- Matrix::colSums(kernel)

#   #Norm
#   D[D != 0] = 1 / D[D != 0]

#   #CSR
#   T <- Matrix::sparseMatrix(
#     i = seq_len(N),
#     j = seq_len(N),
#     x = D,
#     dims = c(N, N)
#   )
#   TK <- T %c*% kernel

#   #Eigen Value Decomp
#   e <- RSpectra::eigs(TK, k = n_components)
#   eD <- as.numeric(e$values)
#   eV <- apply(e$vectors, 2, as.numeric)
#   idx <- order(eD,decreasing=TRUE)
#   eD <- eD[idx]
#   eV <- eV[,idx]

#   #Normalize
#   eV <- apply(eV, 2, function(x){
#     x / sqrt(sum(x^2))
#   })

#   #Create Object
#   out <- SimpleList(
#     T = TK,
#     EigenVectors = eV,
#     EigenValues = eD,
#     Kernel = kernel  
#   )

#   out

# }

#' @export
determine_multiscale_space <- function(
  dm_res = NULL, 
  n_eigs = NULL
  ){
  idx <- seq(2, n_eigs)
  eig_vals <- dm_res[["EigenValues"]][idx]
  data <- t(t(dm_res[["EigenVectors"]][, idx]) * (eig_vals / (1 - eig_vals)))
  data  
}

