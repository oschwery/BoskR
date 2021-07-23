#' Rearrange Input Trees
#'
#' Rearranges different input trees to correct format for downstream analyses.
#'
#' This function accepts different kinds of input phylogenies (see below) and rearranges them into a format that will work for the remaining functions of the package.
#'
#' @param trees Vector of tree objects: individual trees, lists of trees, and/or objects of class multiPhylo. Also accepts the output of `GetMetricTrees` from specified parameters.
#' @param sims Logical, `FALSE` if combining empirical trees, `TRUE` if combining simulations based on an empirical tree. This setting is mainly for the case where trees were simulated under a model that is not implemented, so they can be supplied to `GetTreeMetrics`. Default is `FALSE`.
#' @return List of trees in correct format to be used by downstream functions.
#'
#' @export
#'
#' @import ape

CombineTrees <- function(trees, sims=FALSE) {
  outtrees <- list()
  if (class(trees) == "phylo") {  # turn single tree into a list containing that tree, class multiPhylo
    temptree <- list(trees)
    class(temptree) <- "multiPhylo"
    outtrees <- temptree
  } else if (class(trees) == "multiPhylo") {  # keep multiPhylo as they are
    if (sims==TRUE) {
      temptrees <- list()
      temptrees$metricTreeSet <- list(trees)
      outtrees <- temptrees
    } else {
      outtrees <- trees
    }
  } else if (class(trees) == "list") {  # turn list of trees into...
    for (i in 1:length(trees)) {
     if (class(trees[[i]]) == "phylo") {  # if list of trees: turn all into list of multiPhylo
       temptree <- list(trees[[i]])
       class(temptree) <- "multiPhylo"
       outtrees <- c(outtrees, temptree)
       temptree <- list()
     } else if (class(trees[[i]]) == "multiPhylo") {  # if list of multiPhylo objects, just append them to list
       outtrees <- c(outtrees, trees[[i]])
     } else if (is.list(trees[[i]]) & class(trees[[i]][[1]]) == "phylo") {  # if list of phylo, turn multiphylo and append
       class(trees[[i]]) <- "multiPhylo"
       outtrees <- c(outtrees, trees[[i]])
     }
   }
  }
  if (sims==FALSE) {
      class(outtrees) <- "multiPhylo"
  }
  outtrees
}
