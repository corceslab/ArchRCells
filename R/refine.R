#' @export
refineArchRCells <- function(
  matDR = NULL, 
  ids = NULL, 
  meanCells = 100,
  minCells = 10, 
  maxCells = 150,
  maxIter = 10
  ){

  idList <- list(as.vector(ids))

  it <- 0

  #Refine Small Groups
  ids <- refine_cells_small(
    matDR = matDR, 
    ids = ids, 
    minCells = minCells
  )

  while(max(table(ids)) > maxCells & it < maxIter){

    it <- it + 1
    message(it)

    #1. Refine Large Groups
    ids <- refine_cells_large(
      matDR = matDR, 
      ids = ids, 
      meanCells = meanCells, 
      maxCells = maxCells
    )

    #2. Refine Small Groups
    ids <- refine_cells_small(
      matDR = matDR, 
      ids = ids, 
      minCells = minCells
    )

    idList[[length(idList) + 1]] <- paste0(match(ids, unique(ids)))

  }

  dfLabels <- DataFrame(Reduce("cbind", idList))
  colnames(dfLabels) <- paste0("ArchRCells_", (1:ncol(dfLabels))-1)
  dfLabels <- dfLabels[,rev(1:ncol(dfLabels)),drop=FALSE]
  colnames(dfLabels)[1] <- "ArchRCells"
  dfLabels

}

refine_cells_small <- function(matDR, ids, minCells = 10){

  message("Refining Small ArchRCells...")

  #1. Count Occurences
  tabId <- sort(table(ids))

  #2. Determine Ids To Re-Assign
  fail <- names(which(tabId < minCells))
  idx1 <- ids %in% fail

  #3. Mean Clusters
  meanDR <- t(group_means(t(matDR), ids))

  #4. KNN
  knn_idx <- nabor::knn(meanDR[!(rownames(meanDR) %in% fail),,drop=FALSE], matDR[idx1,,drop=FALSE], k = 1)[[1]]

  #5. Get New Ids
  new_ids <- rownames(meanDR[!(rownames(meanDR) %in% fail),,drop=FALSE])[knn_idx]

  #6. Assign
  ids[idx1] <- new_ids

  ids

}

refine_cells_large <- function(matDR, ids, meanCells = 100, maxCells = 200){

  message("Refining Large ArchRCells...")

  #1. Count Occurences
  tabId <- sample(table(ids))

  #2. Determine Ids To Re-Assign
  fail <- names(which(tabId > maxCells))

  #3. K-Means Assign To New Cells
  fail_new <- parallel::mclapply(seq_along(fail), function(x){

    if(x %% 10 == 0) message(x , " of ", length(fail), " ArchRCells")

    #Idx
    i <- fail[x]

    #N
    n <- ceiling(sum(ids == i) / meanCells)

    #Kmeans
    km <- kmeans(
      x = matDR[ids == i,,drop=FALSE],
      centers = n,
      iter.max = 100L,
      nstart = 10L
    )

    vals <- paste0(i,".",km$cluster)
    names(vals) <- rownames(matDR)[ids == i]
    vals

  }, mc.cores = 8, mc.preschedule = FALSE) %>% unlist 

  #4. Assign
  names(ids) <- rownames(matDR)
  ids[names(fail_new)] <- as.vector(fail_new)

  ids

}