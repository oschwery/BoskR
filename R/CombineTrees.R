#' Rearrange Input Trees
#'
#' Rearranges different input trees to correct format for downstream analyses.
#'
#' This function accepts different kinds of input phylogenies (see below) and rearranges them into a format that will work for the remaining functions of the package.
#'
#' @param trees Vector of tree objects: individual trees, lists of trees, and/or objects of class multiPhylo. Also accepts the output of `GetMetricTrees` from specified parameters.
#' @return List of trees in correct format to be used by downstream functions.

CombineTrees <- function(trees) {
  outtrees <- list()
  outnames <- list()
  for (i in 1:length(trees)) {

  }
  outtrees
}
