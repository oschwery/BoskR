#' Rearrange Input Trees
#'
#' Rearranges different input trees to correct format for downstream analyses.
#'
#' This function accepts different kinds of input phylogenies (see below) and rearranges them into a format that will work for the remaining functions of the package.
#'
#' @param trees Vector of tree objects: individual trees, lists of trees, and/or objects of class multiPhylo. Also accepts the output of `GetMetricTrees` from specified parameters.
#' @return List of trees in correct format to be used by downstream functions.
#'
#' @export

CombineTrees <- function(trees) {
  outtrees <- list()
  if (class(trees) == "phylo") {
    temptree <- list(trees)
    class(temptree) <- "multiPhylo"
    outtrees <- temptree
  } else if (class(trees) == "multiPhylo") {
    outtrees <- trees
  } else if (class(trees) == "list") {
    for (i in 1:length(trees)) {
     if (class(trees[[i]]) == "phylo") {
       temptree <- list(trees[[i]])
       class(temptree) <- "multiPhylo"
       outtrees <- c(outtrees, temptree)
       temptree <- list()
     } else if (class(trees[[i]]) == "multiPhylo") {
       outtrees <- c(outtrees, trees[[i]])
     } else if (is.list(trees[[i]]) & class(trees[[i]][[1]]) == "phylo") {
       class(trees[[i]]) <- "multiPhylo"
       outtrees <- c(outtrees, trees[[i]])
     }
   }
  }
  class(outtrees) <- "multiPhylo"
  outtrees
}
