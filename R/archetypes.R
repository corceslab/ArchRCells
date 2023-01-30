#' @export
initialize_archetypes <- function(
  matDR = NULL,
  K = NULL,
  k = 15,
  knn = NULL,
  nCells = NULL,
  nEigs = 10,
  pw = NULL,
  threads = 8,
  verbose = TRUE
  ){

  #1. Diffustion Maps
  dm_res <- diffusion_maps(matDR = matDR, knn = knn, n_components = k)

  #2. Multiscale Space
  dc_components <- determine_multiscale_space(dm_res = dm_res, n_eigs = nEigs)

  #3. Initialize SEACells via waypoint
  if(pw > 0){
    waypoint_ix <- max_min_sampling(data = dc_components, num_waypoints = nCells, threads = threads, verbose = verbose)
    if(pw < 1){
      waypoint_ix <- sort(sample(waypoint_ix, floor(length(waypoint_ix) * pw)))
    }
  }else{
    waypoint_ix <- c()
  }

  #4. Greedy Cells
  if(pw < 1){
    nGreedy <- nCells - length(waypoint_ix) + 10
    greedy_ix <- get_greedy_centers(K, n_mcs = nGreedy) 
  }else{
    greedy_ix <- c()
  }

  #5. Combine
  all_ix <- c(waypoint_ix, greedy_ix)
  all_ix <- unique(all_ix)
  all_ix <- head(all_ix, nCells)
  all_ix

}

#' @export
max_min_sampling <- function(
  data = NULL, 
  num_waypoints = NULL,
  threads = 8,
  verbose = TRUE
  ){

  waypoint_set <- c()
  no_iterations <- round(((num_waypoints) / ncol(data)))

  # Sample along each component
  N <- nrow(data)

  #Parallel Sampling
  waypoint_set <- parallel::mclapply(seq_len(ncol(data)), function(ind){

      # Data vector
      vec <- data[, ind]

      # Random initialzlation
      iter_set <- sample(seq_len(N), 1)

      #Waypoints
      determine_waypoints_cpp(vec, iter_set, N, no_iterations, ind, verbose)[,1] + 1

  }, mc.cores = threads) %>% unlist

}

determine_waypoints_R <- function(v, idx, N, it){

  # Random initialzlation
  iter_set <- idx

  # Distances along the component
  dists <- matrix(0, nrow = N, ncol = it)
  dists[,1] <- abs(vec - vec[iter_set])

  #Determine
  for(k in seq(1, it-1)){

    if(k %% 100 == 0) message(k)

      # Minimum distances across the current set
      min_dists <- rowMins(dists[, 1:k,drop=FALSE])

      # Point with the maximum of the minimum distances is the new waypoint
      new_wp <- which.max(min_dists)

      # Append
      iter_set <- c(iter_set, new_wp)

      # Update distances
      dists[, k+1] <- abs(vec - vec[new_wp])

  }

  iter_set

}

#' @export
get_greedy_centers <- function(
  ATA = NULL, 
  n_mcs = NULL,
  verbose = TRUE
  ){

  n <- ncol(ATA)

  # initialize
  # X = K
  # ATA = K
  f <- sparseMatrixStats::rowSums2(ATA^2)
  g <- Matrix::diag(ATA)
  # d = as(matrix(0, nrow = n_mcs, ncol = n), "dgCMatrix")
  omega <- as(matrix(0, nrow = n_mcs, ncol = n), "dgCMatrix")

  # keep track of selected indices
  centers <- rep(0, n_mcs)

  # sampling
  for(j in seq_len(n_mcs)){

    if(j %% 50 == 0 & verbose) message("Getting Greedy Center : ", j , " of ", length(n_mcs))

    score <- as.vector(f / g)
    p <- which.max(score)

    residual <- sum(f)

    delta_term1 <- ATA[, p]
    delta_term2 <- Matrix::rowSums(as.matrix(t(omega) * omega[, p, drop = FALSE]))
    delta <- as(delta_term1 - delta_term2, "sparseVector")

    # some weird rounding errors
    delta[p] <- pmax(delta[p], 0)

    o <- as(delta / max(sqrt(delta[p]), 10^-6), "sparseVector")
    omega_square_norm <- np_linalg_norm(o)^2
    omega_hadamard <- o * o
    term1 <- omega_square_norm * omega_hadamard

    pl <- colSums(omega * as.vector(omega %*% o))

    ATAo <- ATA %*% o
    term2 <- o * (ATAo - pl)

    # update f
    f <- f + (-2 * term2 + term1)

    # update g
    g <- g + omega_hadamard

    # store omega and delta
    # d[j, ] = delta
    omega[j, , drop = FALSE] <- o

    # add index
    centers[j] <- p

  }

  centers
  
}
