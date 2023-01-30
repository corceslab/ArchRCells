#' @export
RunArchRCells <- function(
  matDR = NULL,
  k = 15,
  nCells = ceiling(nrow(matDR) / 100),
  ni = 20,
  gi = 25,
  nEigs = 10,
  pw = 1,
  ai = 2,
  epsilon = 10^-5,
  threads = nEigs - 1,
  version = 2,
  verbose = TRUE
  ){

  tstart <- Sys.time()

  #1. KNN Euclidean For RBF
  dt <- round(as.numeric(difftime(Sys.time(), tstart, units = "min")), 3)
  if(verbose) message("1. Running KNN : ", dt, " min...")
  knn <- fast_knn(
    X = matDR, 
    k = k, 
    returnGraph = TRUE, 
    metric = "euclidean", 
    threads = threads, 
    verbose = verbose
  )

  #2. Create RBF Kernel
  dt <- round(as.numeric(difftime(Sys.time(), tstart, units = "min")), 3)
  if(verbose) message("2. Creating RBF Kernel : ", dt, " min...")
  K <- radial_basis_fn(
    matDR = matDR, 
    knn = knn
  )

  #3. Initialize Archetypes
  dt <- round(as.numeric(difftime(Sys.time(), tstart, units = "min")), 3)
  if(verbose) message("3. Initializing Archetypes : ", dt, " min...")
  archetypes <- initialize_archetypes(
    matDR = matDR,
    K = K,
    k = k,
    knn = knn,
    nCells = nCells,
    nEigs = nEigs,
    pw = pw, 
    threads = threads,
    verbose = verbose
  )

  #Clear
  rm(knn)
  gc()

  #4. Fit Cells
  dt <- round(as.numeric(difftime(Sys.time(), tstart, units = "min")), 3)
  if(verbose) message("4. Fitting SEACells : ", dt, " min...")
  obj <- fit_cells(
    K = K,
    matDR = matDR,
    archetypes = archetypes,
    ni = ni,
    gi = gi,
    ai = ai,
    tstart = tstart,
    version = version,
    verbose = verbose
  )

  obj

}
