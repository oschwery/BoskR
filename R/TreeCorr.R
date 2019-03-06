#' Run Tests and Corrections for Trees
#'
#' Tests input treeset for branch length rounding errors, zero length branches, and order.
#'
#' The function is a wrapper around the internals `CorrUltramet`, `CorrZerobranch`, and `ReorderCladewise`. Trees which are not ultrametric due to rounding errors are being corrected using `nnls.tree` as [discribed on the phytools blog](http://blog.phytools.org/2016/08/fixing-ultrametric-tree-whose-edges-are.html), polytomies are randomly resolved and all trees are reordered to 'cladewise' using the ape functions `multi2di` and `reorder.phylo` respectively.
#'
#' @param emptrees Tree or list of trees.
#' @return Same tree set as input, but corrected if necessary.

TreeCorr <- function(emptrees) {
  emptrees <- CorrUltramet(emptrees)
  emptrees <- CorrZerobranch(emptrees)
  emptrees <- ReorderCladewise(emptrees)
  emptrees
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

CorrUltramet <- function(emptrees) {
  for (i in 1:length(emptrees)) {
    tree <- emptrees[[i]]
    nnls <- c()
    if (is.ultrametric(tree) == FALSE) {
      print(paste("not ultrametric", names(emptrees[i]), sep=" "))
      nnls<-nnls.tree(cophenetic(tree),tree,rooted=TRUE)
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

#' Correct Zelo Length Branches/Polytomies
#'
#' Tests input treeset for branches of lenght zero/politomies and randomly resolves them.
#'
#' The internal function tests whether trees have polytomies and randomly resolves them using the ape functions `multi2di`.
#'
#' @param emptrees Tree or list of trees.
#' @return Same tree set as input, but corrected if necessary.
#'
#' @noRd

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

ReorderCladewise <- function(emptrees) {
  for (i in 1:length(emptrees)) {
    emptrees[[i]] <- reorder.phylo(emptrees[[i]],"cladewise")
  }
  emptrees
}
