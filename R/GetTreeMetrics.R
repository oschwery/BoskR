#' Get metrics describing tree shape
#'
#' `GetTreeMetrics` calculates a number of metrics describing tree shape for a tree or a set of trees.
#'
#' The function wraps around the internal 'GetMetrics', which will calculate five 'traditional' tree metrics (Colless, Sackin, number of cherries, number of pitchforks, ladder sizes), the node age and branch length-based metrics (number of taxa, minimum, median, and maximum node age and branch length respectively), as well as standard and normalised graph Laplacian spectra and the associated summary metrics (principal eigenvalue, asymmetry, peakedness, eigengap), as implemented in `RPANDA`.
#'
#' @param trees Tree or set of trees, list or multiPhylo-object, or list of tree sets
#' @param empirical_start `TRUE` if started out from empirical trees, `FALSE` if started from user-specified parameters (currently just affects status updates).
#' @param return_spectra `TRUE` to return the full Laplacian spectra of each tree, `FALSE` (the default) to only return the list of metrics (traditional and based on spectra).
#' @param quick_run If `TRUE`, ommits the time and memory intensive step of calculating the Laplacian spectra and estimates the remaining metrics only (default is `FALSE`).
#' @return A list with two elements: `metrics`: a matrix with the values for all tree metrics for each tree, and `spectra`: a list of raw values for the standard and normalised graph Laplacian spectra for each tree. If applied to the simulated trees based on a tree set, it will be one such two-element list for each tree set provided in a nested list. If `return_spectra` is `FALSE`, only `metrics` will be returned.
#'
#' @export
#'
#' @import ape

GetTreeMetrics <- function(trees, empirical_start=FALSE, return_spectra=FALSE, quick_run=FALSE) {
  outlist <- c()
  if (!is.list(trees[[1]][[1]][[1]]) & !is.na(trees[[1]][[1]][[1]])[1]) {
    outlist <- try(GetMetrics(trees, empirical_start, return_spectra, quick_run), FALSE)
  } else if (is.list(trees[[1]][[1]][[1]]) | is.na(trees[[1]][[1]][[1]])[1]) {
    for (k in 1:length(trees$metricTreeSet)) {
      metricname <- c()
      metrics <- c()
      print(names(trees$metricTreeSet[k]))
      if (is.na(trees$metricTreeSet[[k]])[1]) {
        metrics <- NA
      } else {
        metrics <- try(GetMetrics(trees$metricTreeSet[[k]], empirical_start, return_spectra, quick_run), FALSE)
      }
      metricname <- paste("metrics", names(trees$metricTreeSet[k]), sep="_")
      if (class(metrics) == "try-error") {
        metrics <- NA
      }
      outlist[[metricname]] <- metrics
    }
  }
  outlist
}

#' Internal to get metrics describing tree shape
#'
#' `GetMetrics` calculates a number of metrics describing tree shape for a tree or a set of trees, within the wrapper of 'GetTreeMetrics', which allows to use it for lists of tree sets (e.g. simulated trees based on a set of trees)
#'
#' The function will calculate five 'traditional' tree metrics (Colless, Sackin, number of cherries, number of pitchforks, ladder sizes), the node age and branch length-based metrics (number of taxa, minimum, median, and maximum node age and branch length respectively), as well as standard and normalised graph Laplacian spectra and the associated summary metrics (principal eigenvalue, asymmetry, peakedness, eigengap), as implemented in `RPANDA`.
#'
#' @param trees Tree or set of trees, list or multiPhylo-object
#' @param empirical_start `TRUE` if started out from empirical trees, `FALSE` if started from user-specified parameters (currently just affects status updates).
#' @param return_spectra `TRUE` to return the full Laplacian spectra of each tree, `FALSE` to only return the list of metrics (traditional and based on spectra). This defaults to `FALSE` if `quick_run` is `TRUE`.
#' @param quick_run If `TRUE`, ommits the time and memory intensive step of calculating the Laplacian spectra and estimates the remaining metrics only.
#' @return A list with two elements: `metrics`: a matrix with the values for all tree metrics for each tree, and `spectra`: a list of raw values for the standard and normalised graph Laplacian spectra for each tree. If `return_spectra` is `FALSE`, only `metrics` will be returned.
#'
#' @noRd
#'
#' @import ape
#' @import RPANDA
#' @import phyloTop

GetMetrics <- function(trees, empirical_start=FALSE, return_spectra=FALSE, quick_run=FALSE) {
  metricsnames <- c("Colless", "Sackin", "Cherries", "pitchforks", "AvgLadder", "Gamma", "N_tax", "Min_NodeAge", "Median_NodeAge", "Max_NodeAge", "Min_BranchLength", "Median_BranchLength", "Max_BranchLength")
  if (quick_run != TRUE) {
    metricsnames <- c(metricsnames, "Princ_Eigenv_St", "Asymmetry_St", "Peakedness_St", "Eigengap_St", "Princ_Eigenv_Nor", "Asymmetry_Nor", "Peakedness_Nor")
  }
  metricsmatrix <- matrix(nrow=length(trees), ncol=length(metricsnames))
  colnames(metricsmatrix) <- metricsnames
  rownames(metricsmatrix) <- names(trees)
  spectrallist <- list()
  for (i in 1:length(trees)) {
    tree <- c()
    tree <- trees[[i]]
    #Colless
    metricsmatrix[i, 1] <- colless.phylo(tree, normalise = TRUE)
    #Sackin
    metricsmatrix[i, 2] <- sackin.phylo(tree, normalise = FALSE)
    ######################################total cophenetic indices
    #cherries
    metricsmatrix[i, 3] <- cherries(tree, normalise = FALSE)
    #pitchforks
    metricsmatrix[i, 4] <- pitchforks(tree, normalise = FALSE)
    #ladder sizes
    metricsmatrix[i, 5] <- avgLadder(tree, normalise = FALSE)
    # gamma statistic
    metricsmatrix[i, 6] <- gammaStat(tree)
    # Number of tips
    metricsmatrix[i, 7] <- length(tree$tip.label)
    # node ages
    ages <- branching.times(tree)
    metricsmatrix[i, 8] <- min(ages)
    metricsmatrix[i, 9] <- median(ages)
    metricsmatrix[i, 10] <- max(ages)
    # branch lengths
    metricsmatrix[i, 11] <- min(tree$edge.length)
    metricsmatrix[i, 12] <- median(tree$edge.length)
    metricsmatrix[i, 13] <- max(tree$edge.length)
    # RPANDA metrics
    if(quick_run != TRUE) {
      # standard spectral
      if (return_spectra == TRUE) {
        standardspec <- spectR(tree, meth="standard")
      } else if (return_spectra == FALSE) {
        standardspec <- spectR(tree, meth="standard")[2:5]
      }
      metricsmatrix[i, 14] <- standardspec$principal_eigenvalue
      metricsmatrix[i, 15] <- standardspec$asymmetry
      metricsmatrix[i, 16] <- standardspec$peakedness
      metricsmatrix[i, 17] <- standardspec$eigengap
      if (return_spectra == FALSE) {
        standardspec <- NULL  #empty this before next round
        gc()
      }
      #normalised spectral (disregard eigengap)
      if (return_spectra == TRUE) {
        normalspec <- spectR(tree, meth="normal")
      } else if (return_spectra == FALSE) {
        normalspec <- spectR(tree, meth="normal")[2:4]
      }
      metricsmatrix[i, 18] <- normalspec$principal_eigenvalue
      metricsmatrix[i, 19] <- normalspec$asymmetry
      metricsmatrix[i, 20] <- normalspec$peakedness
      if (return_spectra == FALSE) {
        normalspec <- NULL  #empty this before next round
        gc()
      }
      # drop to list
      if (return_spectra == TRUE) {
        spectra <- list(standardspec=standardspec, normalspec=normalspec)
        spectrallist[[i]] <- spectra
        spectra <- NULL  # empty this before next round
        gc()
        }
      standardspec <- NULL  #empty this before next round
      normalspec <- NULL  #empty this before next round
      gc()
    }
    if (empirical_start == TRUE) {
      print(paste("tree", i, "of", length(trees), sep=" "))
    } else if (empirical_start == FALSE) {
      print(paste("set", k, "tree", i, sep=" "))
    }
  }
  if (return_spectra == TRUE) {
    names(spectrallist) <- names(trees)
    outlist <- list(metrics=metricsmatrix, spectra=spectrallist)
    spectrallist <- NULL
    metricsmatrix <- NULL
    gc()
  } else {
    outlist <- list(metrics=metricsmatrix)
    metricsmatrix <- NULL
    gc()
  }
  outlist
}
