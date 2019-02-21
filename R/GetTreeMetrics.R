#' Get metrics describing tree shape
#'
#' `GetTreeMetrics` calculates a number of metrics describing tree shape for a tree or a set of trees.
#'
#' The function will calculate five 'traditional' tree metrics (Colless, Sackin, number of cherries, number of pitchforks, ladder sizes), as well as standard and normalised graph Laplacian spectra and the associated summary metrics (principal eigenvalue, asymmetry, peakedness, eigengap), as implemented in `RPANDA`.
#'
#' @param trees Tree or set of trees, list or multiPhylo-object
#' @param empirical_start `TRUE` if started out from empirical trees, `FALSE` if started from user-specified parameters
#' @return A list with two elements: `metrics`: a matrix with the values for all tree metrics for each tree, and `spectra`: a list of raw values for the standard and normalised graph Laplacian spectra for each tree.
GetTreeMetrics <- function(trees, empirical_start=FALSE) {
  metricsmatrix <- matrix(nrow=length(trees), ncol=12)
  colnames(metricsmatrix) <- c("Colless", "Sackin", "Cherries", "pitchforks", "AvgLadder", "Princ_Eigenv_St", "Asymmetry_St", "Peakedness_St", "Eigengap_St", "Princ_Eigenv_Nor", "Asymmetry_Nor", "Peakedness_Nor")
  rownames(metricsmatrix) <- names(trees)
  spectrallist <- list()
  for (i in 1:length(trees)) {
    tree <- c()
    if (empirical_start == TRUE) {
      tree <- trees[[i]]
    } else if (empirical_start == FALSE) {
      tree <- trees[[i]]
    }
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
    # RPANDA metrics
    # standard spectral
    standardspec <- spectR(tree, method="standard")
    metricsmatrix[i, 6] <- standardspec$principal_eigenvalue
    metricsmatrix[i, 7] <- standardspec$asymmetry
    metricsmatrix[i, 8] <- standardspec$peakedness
    metricsmatrix[i, 9] <- standardspec$eigengap
    #normalised spectral (disregard eigengap)
    normalspec <- spectR(tree, method="normal") #disabled until fixed
    metricsmatrix[i, 10] <- normalspec$principal_eigenvalue
    metricsmatrix[i, 11] <- normalspec$asymmetry
    metricsmatrix[i, 12] <- normalspec$peakedness
    # drop to list
    spectra <- list(standardspec=standardspec, normalspec=normalspec)
    spectrallist[[i]] <- spectra
    standardspec <- NULL  #empty this before next round
    normalspec <- NULL  #empty this before next round
    spectra <- NULL  # empty this before next round
    if (empirical_start == TRUE) {
      print(paste("tree", i, "of", length(emptrees), sep=" "))
    } else if (empirical_start == FALSE) {
      print(paste("set", k, "tree", i, sep=" "))
    }
  }
  names(spectrallist) <- names(trees)
  outlist <- list(metrics=metricsmatrix, spectra=spectrallist)
  outlist
}
