
#' @export
print.pathsInfo <- function(x){
  cat("\n")
  cat("Paths Parameters")
  cat("\n\n")
  cat(x$type,":",x$A,"; outcome:", x$Y,"; covariates:", x$cov_x)
  cat("\n\n")
  cat("mediators:")
  cat("\n\n")
  print(x$M.list)
}


#' @export
print.pathsEffect <- function(x, ...) {
  cat("\n")
  cat("Causal Paths Analysis: ")
  cat(x$decomposition,"; scale:", x$scale, "; estimation:", x$estimation)
  cat("\n\n")

  # Print function call
  cat("Call: ")
  print(x$call)
  cat("\n\n")

  # print results
  if(is.null( x$boot_results)){
  cat("Causal Paths Estimates: \n\n")
  print(x$results)
  cat("\n\n")}

  if(!is.null( x$boot_results)){
    cat("boot strap results: \n\n")
    print(x$boot_results)
    cat("\n\n")
  }

  #
}



#' @export
print.flexPotential <- function(x, ...) {
  cat("\n")
  # Print function call
  cat("Call: ")
  print(x$call)
  cat("\n")



  cat("Causal Paths Analysis: \n")
  cat("mean counterfactual outcome for flexible path \n")
  cat("estimation:", x$estimation,"\n\n")


  print(x$results, ...)
}



#' @export
print.mflexPotential <- function(x, ...) {
  cat("\n")
  # Print function call
  cat("Call: ")
  print(x$call)
  cat("\n")



  cat("Causal Paths Analysis: \n")
  cat("mean counterfactual outcome for flexible path \n")
  cat("estimation:", x$estimation,"\n\n")


    cat("mutiple treatment matrix: \n")
    print(x$active)
    cat("\n\n")



  print(x$results, ...)
}


#' @export
print.flexEffect <- function(x, ...) {
  cat("\n")
  cat("Flexible Paths Analysis ")
  cat("\n\n")

  # Print function call
  cat("Call: ")
  print(x$call)
  cat("\n\n")

  # print results
  if(is.null( x$boot_results)){
    cat("Causal Paths Estimates: \n\n")
    print(x$results)
    cat("\n\n")}

  if(!is.null( x$boot_results)){
    cat("boot strap results: \n\n")
    print(x$boot_results)
    cat("\n\n")
  }

  #
}





