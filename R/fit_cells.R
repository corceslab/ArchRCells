#' @export
fit_cells <- function(
  K = NULL, 
  matDR = NULL, 
  archetypes = NULL,
  ni = 10,
  gi = 30,
  epsilon = 10^-5,
  ai = 5,
  version = 2,
  tstart = Sys.time(),
  tmpA = tempfile(pattern = "matrix_A_", tmpdir = "tmpArchRCells", fileext = ".rds"),
  tmpB = tempfile(pattern = "matrix_B_", tmpdir = "tmpArchRCells", fileext = ".rds"),
  verbose = TRUE
  ){

  #Create Tmp Dir
  dir.create(dirname(tmpA), showWarnings = FALSE)

  # initialize B (update this to allow initialization from RRQR)
  n <- ncol(K)
  k <- length(archetypes)

  # Create B Matrix
  B <- Matrix::sparseMatrix(
    i = archetypes,
    j = seq_along(archetypes),
    x = rep(1, length(archetypes)),
    dims = c(n, k)
  )
  dt <- round(as.numeric(difftime(Sys.time(), tstart, units = "min")), 3)
  if(verbose) message("Initialized B : ", dt, " min...")

  #Iterate A
  t2 <- Matrix::t(K %c*% B)
  t1 <- (t2 %c*% B)
  for(i in seq_len(ai)){

    dt <- round(as.numeric(difftime(Sys.time(), tstart, units = "min")), 3)
    if(verbose) message("Randomly Initializing A ( ", i , " of ", ai, " ) : ",  dt, " min...")

    A_i <- sparseMatrix(
      i = initialize_A_idx(t1 = t1, t2 = t2, verbose = verbose) + 1, 
      j = seq_len(nrow(B)), 
      x = rep(1, nrow(B)), 
      dims = c(ncol(B), nrow(B))
    )

    RSS_i <- compute_RSS(A=A_i, B=B, matDR=matDR)

    if(i == 1){
      A <- A_i
      RSS_min <- RSS_i
    }else{
      if(RSS_i < RSS_min){
        A <- A_i
        RSS_min <- RSS_i        
      }
    }
    rm(A_i)

  }
  rm(t2, t1)
  gc()
  dt <- round(as.numeric(difftime(Sys.time(), tstart, units = "min")), 3)
  if(verbose) message("Initialized A : ", dt, " min...")

  #Shallow Update A Since its Dense
  A <- update_A(A=A, B=B, K=K, it=gi, t=1, version=2, verbose=verbose)

  #Print
  dt <- round(as.numeric(difftime(Sys.time(), tstart, units = "min")), 3)
  if(verbose) message("Update A Initial : ", dt, " min...")

  # Create convergence threshold
  RSS_i <- compute_RSS(A=A, B=Matrix(B, sparse=TRUE), matDR=matDR)
  RSS_iters <- c(RSS_i)

  #Treshold
  convergence_threshold <- epsilon * RSS_i

  #Save Initial
  saveRDS(A, tmpA, compress = FALSE)
  saveRDS(B, tmpB, compress = FALSE)

  #Iterate
  converged <- FALSE
  iter <- 0
  while((!converged & iter < ni) | (iter < floor(ni / 2))){
    
    iter <- iter + 1
    dt <- round(as.numeric(difftime(Sys.time(), tstart, units = "min")), 3)
    if(verbose) message("Fitting iteration ", iter, " of ", ni, " : ", dt, " min...")

    #1 Update A     
    A <- update_A(A=A, B=B, K=K, it=gi, t=0, version=version, verbose=verbose)

    #Compute RSS
    RSS_i <- compute_RSS(A=A, B=B, matDR=matDR)

    #If Lowest Save A+B Combo
    if(RSS_i < min(RSS_iters)){
      saveRDS(A, tmpA, compress = FALSE)
      saveRDS(B, tmpB, compress = FALSE)
    }

    #Append
    RSS_iters <- c(RSS_iters, RSS_i)      
  
    dt <- round(as.numeric(difftime(Sys.time(), tstart, units = "min")), 3)
    if(verbose) message("\tUpdate A (RSS = ", round(RSS_i,1), ") ", iter, " of ", ni, " : ", dt, " min...")

    #2 Update B
    B <- update_B(A=A, B=B, K=K, it=gi, t=0, version=version, verbose=verbose)

    #Compute RSS
    RSS_i <- compute_RSS(A=A, B=B, matDR=matDR)

    #If Lowest Save A+B Combo
    if(RSS_i < min(RSS_iters)){
      saveRDS(A, tmpA, compress = FALSE)
      saveRDS(B, tmpB, compress = FALSE)
    }

    #Append
    RSS_iters <- c(RSS_iters, RSS_i)

    dt <- round(as.numeric(difftime(Sys.time(), tstart, units = "min")), 3)
    if(verbose) message("\tUpdate B (RSS = ", round(RSS_i,1), ") ", iter, " of ", ni, " : ", dt, " min...")

    #4 Converged?
    if(abs(tail(RSS_iters,2)[1] - tail(RSS_iters,1)) < convergence_threshold){
      converged <- TRUE
    }
    gc()

    message("\n")

  }

  #Load Minima
  A <- readRDS(tmpA)
  B <- readRDS(tmpB)
  RSS_final <- compute_RSS(A=A, B=B, matDR=matDR)
  if(verbose) message("Final RSS = ", round(RSS_final,3))

  #File Remove
  if(file.exists(tmpA)) file.remove(tmpA)
  if(file.exists(tmpB)) file.remove(tmpB)

  #Assignments
  dfLabel <- get_assignments(A, B)

  #Return
  out <- SimpleList(
    dfLabel = dfLabel,
    A = A,
    B = B,
    K = K,
    RSS = RSS_iters
  )

  out

}

#' @export
update_A <- function(
  A = NULL,
  B = NULL,
  K = NULL,
  it = 20,
  t = 0,
  version = 2,
  verbose = TRUE
){

  if(version == 1){

    A <- update_A_cpp(
      A = Matrix(A, sparse=TRUE), 
      B = Matrix(B, sparse=TRUE), 
      K = Matrix(K, sparse=TRUE), 
      it = it,
      t = t,
      verbose = verbose
    ) 

  }else if(version == 2){

    A <- update_A_cpp_sp(
      A = Matrix(A, sparse=TRUE), 
      B = Matrix(B, sparse=TRUE), 
      K = Matrix(K, sparse=TRUE), 
      it = it,
      t = t,
      verbose = verbose
    ) 

  }else{

    stop("Error 2 versions available!")

  }

  A

}

#' @export
update_B <- function(
  A = NULL,
  B = NULL,
  K = NULL,
  it = 20,
  t = 0,
  version = 2,
  verbose = TRUE
){

  if(version == 1){

    A <- update_B_cpp(
      A = Matrix(A, sparse=TRUE), 
      B = Matrix(B, sparse=TRUE), 
      K = Matrix(K, sparse=TRUE), 
      it = it,
      t = t,
      verbose = verbose
    ) 

  }else if(version == 2){

    A <- update_B_cpp_sp(
      A = Matrix(A, sparse=TRUE), 
      B = Matrix(B, sparse=TRUE), 
      K = Matrix(K, sparse=TRUE), 
      it = it,
      t = t,
      verbose = verbose
    ) 

  }else{

    stop("Error 2 versions available!")

  }

  A

}


# #' @export
# update_A_dense <- function(A, B, K, it = 5){

#   #Update A
#   t2 <- Matrix::t(K %c*% B)
#   t1 <- (t2 %c*% B)
#   t2s <- data.frame(Matrix::summary(t2))

#   #Compute Gradient
#   G <- as.matrix(t1 %c*% A)
#   G[cbind(t2s[,1], t2s[,2])] <- 2*(G[cbind(t2s[,1], t2s[,2])] - t2s[,3])

#   #Get argmins
#   amins <- col_mins_id_cpp(G)[,1] + 1
#   rm(G)
#   gc()

#   #loop free implementation
#   A <- sparseMatrix(
#     i = amins,
#     j = seq_along(amins),
#     x = rep(1, length(amins))
#   )

#   #Sparse A Update
#   for(t in seq_len(it - 1)){

#     # compute gradient 
#     G <- 2 * (Matrix(t1 %c*% A, sparse=TRUE) - t2)

#     # get argmins
#     amins <- apply(G, 2, which.min)

#     #loop free implementation
#     e <- sparseMatrix(
#       i = amins,
#       j = seq_along(amins),
#       x = rep(1, length(amins)),
#       dims = c(nrow(A), ncol(A))
#     )

#     #Update
#     A <- A + 2. / (t + 2.) * (e - A)

#   }

#   A

# }




# #' @export
# fit_cells <- function(
#   K = NULL, 
#   matDR = NULL, 
#   archetypes = NULL,
#   ni = 10,
#   gi = 30,
#   epsilon = 10^-5,
#   ai = 5,
#   tstart = Sys.time(),
#   tmpA = tempfile(pattern = "matrix_A_", tmpdir = "tmpArchRCells", fileext = ".rds"),
#   tmpB = tempfile(pattern = "matrix_B_", tmpdir = "tmpArchRCells", fileext = ".rds")
#   ){

#   #Create Tmp Dir
#   dir.create(dirname(tmpA), showWarnings = FALSE)

#   # initialize B (update this to allow initialization from RRQR)
#   n <- ncol(K)
#   k <- length(archetypes)

#   # Create B Matrix
#   B <- Matrix::sparseMatrix(
#     i = archetypes,
#     j = seq_along(archetypes),
#     x = rep(1, length(archetypes)),
#     dims = c(n, k)
#   )
#   dt <- round(as.numeric(difftime(Sys.time(), tstart, units = "min")), 3)
#   message("Initialized B : ", dt, " min...")

#   #Create A Matrix
#   A <- initialize_A_cpp(k, n)

#   #RSS
#   RSS_min <- compute_RSS(A=A, B=B, matDR=matDR)

#   #Iterate
#   for(i in seq_len(ai)){
#     A_i <- initialize_A_cpp(k, n)
#     RSS_i <- compute_RSS(A=A_i, B=B, matDR=matDR)
#     if(RSS_i < RSS_min){
#       RSS_min <- RSS_i
#       A <- A_i
#     }
#     rm(A_i)
#   }
#   gc()
#   dt <- round(as.numeric(difftime(Sys.time(), tstart, units = "min")), 3)
#   message("Initialized A : ", dt, " min...")

#   #Shallow Update A Since its Dense
#   A <- update_A_dense(
#     A = A, 
#     B = Matrix(B, sparse=TRUE), 
#     K = Matrix(K, sparse=TRUE), 
#     it = 5
#   )

#   dt <- round(as.numeric(difftime(Sys.time(), tstart, units = "min")), 3)
#   message("Update A Shallow : ", dt, " min...")

#   # Create convergence threshold
#   RSS_i <- compute_RSS(A=A, B=Matrix(B, sparse=TRUE), matDR=matDR)
#   RSS_iters <- c(RSS_i)

#   #Treshold
#   convergence_threshold <- epsilon * RSS_i

#   #Save Initial
#   saveRDS(A, tmpA, compress = FALSE)
#   saveRDS(B, tmpB, compress = FALSE)

#   #Iterate
#   converged <- FALSE
#   iter <- 0
#   while((!converged & iter < ni) | (iter < floor(ni / 2))){
    
#     iter <- iter + 1
#     dt <- round(as.numeric(difftime(Sys.time(), tstart, units = "min")), 3)
#     message("Fitting iteration ", iter, " of ", ni, " : ", dt, " min...")

#     #1 Update A     
#     A <- update_A_cpp(
#       A = Matrix(A, sparse=TRUE), 
#       B = Matrix(B, sparse=TRUE), 
#       K = Matrix(K, sparse=TRUE), 
#       it = gi
#     )

#     #Compute RSS
#     RSS_i <- compute_RSS(A=A, B=B, matDR=matDR)

#     #If Lowest Save A+B Combo
#     if(RSS_i < min(RSS_iters)){
#       saveRDS(A, tmpA, compress = FALSE)
#       saveRDS(B, tmpB, compress = FALSE)
#     }

#     #Append
#     RSS_iters <- c(RSS_iters, RSS_i)      
  
#     dt <- round(as.numeric(difftime(Sys.time(), tstart, units = "min")), 3)
#     message("\tUpdate A (RSS = ", round(RSS_i,1), ") ", iter, " of ", ni, " : ", dt, " min...")

#     #2 Update B
#     B <- update_B_cpp(
#       A = Matrix(A, sparse=TRUE), 
#       B = Matrix(B, sparse=TRUE), 
#       K = Matrix(K, sparse=TRUE), 
#       it = gi
#     )

#     #Compute RSS
#     RSS_i <- compute_RSS(A=A, B=B, matDR=matDR)

#     #If Lowest Save A+B Combo
#     if(RSS_i < min(RSS_iters)){
#       saveRDS(A, tmpA, compress = FALSE)
#       saveRDS(B, tmpB, compress = FALSE)
#     }

#     #Append
#     RSS_iters <- c(RSS_iters, RSS_i)

#     dt <- round(as.numeric(difftime(Sys.time(), tstart, units = "min")), 3)
#     message("\tUpdate B (RSS = ", round(RSS_i,1), ") ", iter, " of ", ni, " : ", dt, " min...")

#     #4 Converged?
#     if(abs(tail(RSS_iters,2)[1] - tail(RSS_iters,1)) < convergence_threshold){
#       converged <- TRUE
#     }
#     gc()

#     message("\n")

#   }

#   #Load Minima
#   A <- readRDS(tmpA)
#   B <- readRDS(tmpB)
#   RSS_final <- compute_RSS(A=A, B=B, matDR=matDR)
#   message("Final RSS = ", round(RSS_final,3))

#   #File Remove
#   if(file.exists(tmpA)) file.remove(tmpA)
#   if(file.exists(tmpB)) file.remove(tmpB)

#   #Assignments
#   dfLabel <- get_assignments(A, B)

#   #Return
#   out <- SimpleList(
#     dfLabel = dfLabel,
#     A = A,
#     B = B,
#     K = K,
#     RSS = RSS_iters
#   )

# }
