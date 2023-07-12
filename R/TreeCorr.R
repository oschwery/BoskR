#' Run Tests and Corrections for Trees
#'
#' Tests input treeset for branch length rounding errors, zero length branches, and order.
#'
#' The function is a wrapper around the internals `CorrUltramet`, `CorrZerobranch`, and `ReorderCladewise`. Trees which are not ultrametric due to rounding errors are being corrected using `nnls.tree` as [discribed on the phytools blog](http://blog.phytools.org/2016/08/fixing-ultrametric-tree-whose-edges-are.html), polytomies are randomly resolved and all trees are reordered to 'cladewise' using the ape functions `multi2di` and `reorder.phylo` respectively.
#'
#' @param emptrees Tree or list of trees.
#' @return Same tree set as input, but corrected if necessary.
#'
#' @export
#'
#' @import ape

TreeCorr <- function(emptrees) {
  CheckRootedness(emptrees)
  emptrees <- CorrUltramet(emptrees)
  emptrees <- CorrZerobranch(emptrees)
  emptrees <- ReorderCladewise(emptrees)
  emptrees
}

#' Test whether trees are rooted
#'
#' Tests input treeset for whether they are rooted trees.
#'
#' The internal function tests whether all trees are rooted, which they have to be for subsequent analyses. If any are unrooted, the function returns an error and indicates the indices of the unrooted trees.
#'
#' @param emptrees Tree or list of trees
#' @return Either message that all trees are rooted, or error and indices of unrootes trees (to be removed/rooted).
#'
#' @noRd

CheckRootedness <- function(emptrees) {
  rootstatus <- c()
    for (i in 1:length(emptrees)) {
      rootstatus <- c(rootstatus, is.rooted(emptrees[[i]]))
    }
    if (FALSE %in% rootstatus) {
      stop(paste("The following trees are not rooted: ", which(rootstatus==FALSE), sep=""))
    } else {
      print("All trees are rooted.")
    }
}



#' Correct Non-Ultrameitric Trees
#'
#' Tests input treeset for branch length rounding errors.
#'
#' The internal function tests whether trees are not ultrametric due to rounding errors and corrects them using `nnls.tree` as [discribed on the phytools blog](http://blog.phytools.org/2016/08/fixing-ultrametric-tree-whose-edges-are.html).
#'
#' @param emptrees Tree or list of trees.
#' @return Same tree set as input, but corrected if necessary.
#'
#' @noRd
#'
#' @import ape
#' @importFrom phangorn nnls.tree

CorrUltramet <- function(emptrees) {
  for (i in 1:length(emptrees)) {
    tree <- emptrees[[i]]
    nnls <- c()
    if (is.ultrametric(tree) == FALSE) {
      print(paste("not ultrametric", names(emptrees[i]), sep=" "))
      nnls<-nnls.tree(cophenetic(tree), tree, rooted=TRUE)
      if (is.ultrametric(nnls) == TRUE) {
        emptrees[[i]] <- nnls
        print(paste("fixed", names(emptrees[i]), sep=" "))
      } else if (is.ultrametric(nnls) == FALSE) {
        print(paste("Still not ultrametric", names(emptrees[i]), sep=" "))
      }
    }
    tree <- c()
  }
  emptrees
}

#' Correct Zero Length Branches/Polytomies
#'
#' Tests input treeset for branches of lenght zero/politomies and randomly resolves them.
#'
#' The internal function tests whether trees have polytomies and randomly resolves them using the ape functions `multi2di`.
#'
#' @param emptrees Tree or list of trees.
#' @return Same tree set as input, but corrected if necessary.
#'
#' @noRd
#'
#' @import ape

CorrZerobranch <- function(emptrees) {
  for (i in 1:length(emptrees)) {
    tree <- emptrees[[i]]
    if (is.binary(tree) == FALSE) {
      print(paste("not binary", names(emptrees[i]), sep=" "))
      newtree <- multi2di(tree, random=TRUE)
      if (is.binary(newtree) == TRUE) {
        emptrees[[i]] <- newtree
        print(paste("fixed", names(emptrees[i]), sep=" "))
      } else if (is.binary(newtree) == FALSE) {
        print(paste("Still not binary", names(emptrees[i]), sep=" "))
      }
    }
    tree <- c()
  }
  emptrees
}

#' Reorder phylogenies to 'cladewise'
#'
#' Rerders input tree set.
#'
#' The internal function reorders the trees to 'cladewise' using the ape function `reorder.phylo`.
#'
#' @param emptrees Tree or list of trees.
#' @return Same tree set as input, but reordered.
#'
#' @noRd
#'
#' @import ape

ReorderCladewise <- function(emptrees) {
  for (i in 1:length(emptrees)) {
    emptrees[[i]] <- reorder.phylo(emptrees[[i]],"cladewise")
  }
  emptrees
}
