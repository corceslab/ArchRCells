#' @useDynLib ArchRCells

logo <- "
    ___              __    ____  ______     ____    
   /   |  __________/ /_  / __ \\/ ____/__  / / /____
  / /| | / ___/ ___/ __ \\/ /_/ / /   / _ \\/ / / ___/
 / ___ |/ /  / /__/ / / / _, _/ /___/  __/ / (__  ) 
/_/  |_/_/   \\___/_/ /_/_/ |_|\\____/\\___/_/_/____/                                                                              
"

dependency <- c(
  "Rcpp",
  "RcppArmadillo",
  "RcppEigen",
  "uwot",
  "Matrix",
  "S4Vectors",
  "matrixStats",
  "magrittr"
)

.onAttach <- function(libname, pkgname){
  
  cat(logo)
  packageStartupMessage("Logo created from patorjk.com")

  #Package Startup
  v <- packageVersion("ArchRCells")
  packageStartupMessage("ArchRCells : Version ", v)
  
  #Load Packages
  packageStartupMessage("Loading Required Packages...")
  pkgs <- dependency
  for(i in seq_along(pkgs)){
    tryCatch({
      packageStartupMessage("\tLoading Package : ", pkgs[i], " v", packageVersion(pkgs[i]))
      suppressPackageStartupMessages(require(pkgs[i], character.only=TRUE))
    }, error = function(e){
      packageStartupMessage("\tFailed To Load Package : ", pkgs[i])
    })
  }

}
