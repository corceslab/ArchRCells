#' @export
np_linalg_norm <- function(x){
  sqrt(sum(x^2))
}

#' @export
compute_reconstruction <- function(matDR, A, B){
  #(Matrix::t(matDR) %*% B) %*% A
  (Matrix::t(matDR) %c*% B) %c*% A
}

#' @export
compute_RSS <- function(matDR, A, B){
  reconstruction <- compute_reconstruction(matDR, A, B)
  np_linalg_norm(Matrix::t(matDR) - reconstruction)
}

#' @export
binarize_matrix_rows <- function(T){
  idx <- apply(T, 2, which.max)
  sparseMatrix(
    i = idx,
    j = seq_along(idx),
    x = 1,
    dims = c(nrow(T), ncol(T))
  ) 
}

#' @export
get_assignments <- function(A, B){
  
    bin_A <- binarize_matrix_rows(A)
    bin_B <- binarize_matrix_rows(B)

    labels <- (Matrix::t(bin_A) %*% (seq_len(nrow(bin_A))))[,1]

    df <- data.frame(
      ArchRCell = labels,
      MainCell = Matrix::rowSums(bin_B) > 0
    )

    df

}

#' @export
fast_knn <- function(
  X = NULL, 
  k = 15, 
  includeSelf = TRUE,
  method = "uwot",
  metric = "euclidean",
  returnGraph = FALSE,
  binarizeGraph = FALSE,
  threads = 8,
  verbose = TRUE,
  ...
  ){

  if(tolower(method) == "uwot"){

    if(verbose) message("Finding KNN with uwot annoy_nn...")
    
    knn <- uwot:::annoy_nn(
      X = X, 
      k = k, 
      n_threads = threads, 
      metric = metric, 
      verbose = verbose,
      ...
    )
    knn <- SimpleList(
      nn.idx = knn$idx,
      nn.dist = knn$dist,
      nn.recall = knn$recall
    )
    if(!includeSelf){
      if(all(knn$nn.idx[,1]==seq_len(nrow(knn$nn.idx)))){
        knn$nn.idx <- knn$nn.idx[,-1,drop=FALSE]
        knn$nn.dist <- knn$nn.dist[,-1,drop=FALSE]
      }else{
        stop("Not all first indices are self cant remove!")
      }
    }

  }else if(tolower(method) == "seurat"){
    
    if(verbose) message("Finding KNN with Seurat...")

    rownames(X) <- paste0(seq_len(nrow(X)))
    knn <- Seurat::FindNeighbors(
      object = X, 
      k.param = k, 
      annoy.metric = metric,
      return.neighbor = TRUE, 
      verbose = verbose, 
      ...
    )
    knn <- SimpleList(
      nn.idx = knn@nn.idx,
      nn.dist = knn@nn.dists,
    )
    if(!includeSelf){
      if(all(knn$nn.idx[,1]==seq_len(nrow(knn$nn.idx)))){
        knn$nn.idx <- knn$nn.idx[,-1,drop=FALSE]
        knn$nn.dist <- knn$nn.dist[,-1,drop=FALSE]
      }else{
        stop("Not all first indices are self cant remove!")
      }
    }

  }else{
    
    if(verbose) message("Finding KNN with Nabor...")

    knn <- nabor::knn(data = X, k = k, ...)
    knn <- SimpleList(
      nn.idx = knn$nn.idx,
      nn.dist = knn$nn.dists
    )
    if(!includeSelf){
      if(all(knn$nn.idx[,1]==seq_len(nrow(knn$nn.idx)))){
        knn$nn.idx <- knn$nn.idx[,-1,drop=FALSE]
        knn$nn.dist <- knn$nn.dist[,-1,drop=FALSE]
      }else{
        stop("Not all first indices are self cant remove!")
      }
    }

  }

  if(returnGraph){

    if(binarizeGraph){
      
      knn$nn.graph <- sparseMatrix(
        i = rep(seq_len(nrow(X)), times = ncol(knn$nn.idx)),
        j = as.vector(knn$nn.idx),
        x = rep(1, length(as.vector(knn$nn.dist)))
      )

    }else{

      knn$nn.graph <- sparseMatrix(
        i = rep(seq_len(nrow(X)), times = ncol(knn$nn.idx)),
        j = as.vector(knn$nn.idx),
        x = as.vector(knn$nn.dist)
      ) 

    }

  }

  knn

}

#' @export
`%c*%` <- function(X,Y){

    #determine Class
    zzzX <- -1
    zzzY <- -1
    if(is(X, "dgCMatrix")) zzzX <- 1
    if(is(Y, "dgCMatrix")) zzzY <- 1
    if(is(X, "matrix")) zzzX <- 0
    if(is(Y, "matrix")) zzzY <- 0

    if(zzzX < 0) stop()
    if(zzzY < 0) stop()

    if(zzzX == 0 & zzzY == 0){
        mat_mult_dens_dens_cpp(X, Y)
    }else if(zzzX == 1 & zzzY == 0){
        mat_mult_sp_dens_cpp(X, Y)
    }else if(zzzX == 0 & zzzY == 1){
        mat_mult_dens_sp_cpp(X, Y)
    }else if(zzzX == 1 & zzzY == 1){
        mat_mult_sp_sp_cpp(X, Y)
    }else{
        stop()
    }

}
