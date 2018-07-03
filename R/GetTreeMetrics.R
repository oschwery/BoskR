##########################################
# measure the tree metrics
##########################################
GetTreeMetrics <- function(trees, rownum, empirical_start=FALSE) {
  metricsmatrix <- matrix(nrow=rownum, ncol=14)
  colnames(metricsmatrix) <- c("Colless", "Sackin", "Cherries", "pitchforks", "AvgLadder", "Princ_Eigenv_St", "Asymmetry_St", "Peakedness1_St", "Peakedness2_St", "Eigengap_St", "Princ_Eigenv_Nor", "Asymmetry_Nor", "Peakedness1_Nor", "Peakedness2_Nor")
#  rownam <- c()
#  for (i in 1:rownum) {
#    rownam <- c(rownam, paste("tree", i, sep=""))
#  }
  rownames(metricsmatrix) <- names(trees)
  spectrallist <- list()
  for (i in 1:length(trees)) {
    tree <- c()
    if (empirical_start == TRUE) {
#      tree <- trees[[i]]$tree[[1]]
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
    #widths __ paused bcz bigger output
    #widths(tree)
    #node-imbalance
    #node-depth fractions

    # RPANDA
    # standard spectral
    standardspec <- spectR(tree, method="standard")
    metricsmatrix[i, 6] <- standardspec$principal_eigenvalue
    metricsmatrix[i, 7] <- standardspec$asymmetry
    metricsmatrix[i, 8] <- standardspec$peakedness1
    metricsmatrix[i, 9] <- standardspec$peakedness2
    metricsmatrix[i, 10] <- standardspec$eigengap
    #normalised spectral (forget eigengap)
    #normalspec <- spectR(tree, method="normal") #disabled until fixed
    # metricsmatrix[i, 11] <- normalspec$principal_eigenvalue
    # metricsmatrix[i, 12] <- normalspec$asymmetry
    # metricsmatrix[i, 13] <- normalspec$peakedness1
    # metricsmatrix[i, 14] <- normalspec$peakedness2
    # drop to list
#    spectra <- list(standardspec=standardspec, normalspec=normalspec)
    spectra <- list(standardspec=standardspec, normalspec=NA)
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
  # measure the tree metrics
  #LTT
  #branch length distributions
  #AdequacyLTT(tree, trees, PDF=TRUE)
  #AccuracyBoxplot(simulated_solution, PDF=TRUE)
  #BLdistPlot(tree, trees, PDF=TRUE)
  #IntTermBLdistPlot(tree, trees, PDF=TRUE)


  #RPANDA for several trees or for later(BICOMPARE)
  # BICompare  # BIC test for modalities (clusters) in tree, need to be set apriori --> CHECK WHAT THAT MEANS
  # JSDtree  # Jensen-Shannon distance between phylogenies (based on spectral densities)  --> CHECK WHAT THAT MEANS
  # JSDtree_cluster  # clustering & plotting based on JSDtree
  outlist <- list(metrics=metricsmatrix, spectra=spectrallist)
  outlist
}
