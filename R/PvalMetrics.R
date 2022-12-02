#' Get p-values for tree metrics
#'
#' Estimates p-values based on simulated and empirical distributions of tree metrics
#'
#' The function uses an Empirical Cumulative Distribution function to determine the area under the curve of the metric values of the simulated trees, to get to a p-value for the position of the metrics of the empirical tree on that distribution. The argument `metricset` allows to chose between: `"spectR"`- the standard (i.e. unnormalised) spectral densities, `"spectRnorm"` - the normalised spectral densities, `"classic"` - a couple of more 'conventional' measures of tree shape, being Colless index, Sackin index, number of cherries, number of pitchforks, average ladder size, and gamma statistic; finally `"nodibranch" - includes minimum, maximum, and median for both node ages and branch lengths respectively; option `"All"` returns all sets, and `"AllQuick"` returns all except the metrics based on the Laplacian spectra. For more information on the spectral densities, i.e. the Eigenvalues of the tree's modified graph Laplacian, see R package RPANDA and associated papers.
#'
#' @param empMetrics Metrics of empirical tree or set of trees; output of `GetTreeMetrics` or formatted the same way.
#' @param simMetrics Metrics of sets of simulated trees; output of `GetTreeMetrics` or formatted the same way.
#' @param empirical_start Indicator whether empMetrics is based on empirical or simulated initial trees, default is `TRUE` (=empirical); mainly important for data format reasons.
#' @param methodnr Integral specifying which method is used: 1: BD, 2: TimeD-BD, 3: DD; is only used if `empirical_start` is `TRUE`
#' @param metricset String specifying which tree metrics to use; default is "spectR", other options are "spectrRnorm", "classic", and "nodibranch", or "ALL" and "AllQuick"; for more information on the options see Details.
#' @return A list with two entries: `ECDs` is a list of Empirical Cumulative Distributions; `pValues` is a matrix with p-values for the targeted metrics
#'
#' @export

PvalMetrics <- function(empMetrics, simMetrics, empirical_start=TRUE, methodnr, metricset="spectR") {
  # loop getting distributions for all sim trees
  dists <- list()  # primer for ECDs
  targetmetrics <-c()
  if (metricset == "spectR") {
    targetmetrics <- c("Princ_Eigenv_St", "Asymmetry_St", "Peakedness_St", "Eigengap_St")  # CHANGE once you have likelihood implemented!
  } else if (metricset == "spectRnorm") {
    targetmetrics <- c("Princ_Eigenv_Nor", "Asymmetry_Nor", "Peakedness_Nor")
  } else if (metricset == "classic") {
    targetmetrics <- c("Colless", "Sackin", "Cherries", "pitchforks", "AvgLadder", "Gamma")
  } else if (metricset == "nodibranch") {
    targetmetrics <- c("N_tax", "Min_NodeAge", "Median_NodeAge", "Max_NodeAge", "Min_BranchLength", "Median_BranchLength", "Max_BranchLength")
  } else if (metricset == "All") {
    targetmetrics <- c("Colless", "Sackin", "Cherries", "pitchforks", "AvgLadder", "Gamma", "N_tax", "Min_NodeAge", "Median_NodeAge", "Max_NodeAge", "Min_BranchLength", "Median_BranchLength", "Max_BranchLength", "Princ_Eigenv_St", "Asymmetry_St", "Peakedness_St", "Eigengap_St", "Princ_Eigenv_Nor", "Asymmetry_Nor", "Peakedness_Nor")
  } else if (metricset == "AllQuick") {
    targetmetrics <- c("Colless", "Sackin", "Cherries", "pitchforks", "AvgLadder", "Gamma", "N_tax", "Min_NodeAge", "Median_NodeAge", "Max_NodeAge", "Min_BranchLength", "Median_BranchLength", "Max_BranchLength")
  }
  
  Methods <- c("BD", "TimeD-BD", "DD")
  current_method <- Methods[methodnr]
  pval <- matrix(nrow=nrow(empMetrics$metrics), ncol=2*length(targetmetrics))
  if (empirical_start == FALSE) {
    k <- 1
    for (j in 1:nrow(empMetrics$metrics)) {
      dist <- list()  # primer for ECD's
      if (paste("metrics", rownames(empMetrics$metrics)[j], current_method, sep="_") == names(simMetrics[k])) {
        for (i in 1:length(targetmetrics)) {
          deP <- c()  # primer for actual p-Value
          dist[[i]] <- ecdf(simMetrics[[k]]$metrics[, targetmetrics[i]])  # create ecd function for a metric
          pval[j, i] <- empMetrics[1]$metrics[j, targetmetrics[i]]  # add empirical value to be evaluated into output matrix
          deP <- dist[[i]](c(empMetrics[1]$metrics[j, targetmetrics[i]]))  # evaluate ecdf for this value
          pval[j, i+length(targetmetrics)] <- (min(c(deP, 1-deP)))*2  # take smaller value (i.e. low or high end) and multiply by two for two-tailedness
        }
        names(dist) <- targetmetrics
        dists[[j]] <- dist
        k <- k+1
      } else {
        pval[j, ] <- NA
        dists[[j]] <- NA
      }
    }
  } else {
    for (j in 1:nrow(empMetrics$metrics)) {
      dist <- list()  # primer for ECD's
      if (is.na(simMetrics[[j]])) {
        pval[j, ] <- NA
        dists[[j]] <- NA
      } else {
        for (i in 1:length(targetmetrics)) {
          deP <- c()  # primer for actual p-Value
          dist[[i]] <- ecdf(simMetrics[[j]]$metrics[, targetmetrics[i]])  # create ecd function for a metric
          pval[j, i] <- empMetrics[1]$metrics[j, targetmetrics[i]]  # add empirical value to be evaluated into output matrix
          deP <- dist[[i]](c(empMetrics[1]$metrics[j, targetmetrics[i]]))  # evaluate ecdf for this value
          pval[j, i+length(targetmetrics)] <- (min(c(deP, 1-deP)))*2  # take smaller value (i.e. low or high end) and multiply by two for two-tailedness
        }
        names(dist) <- targetmetrics
        dists[[j]] <- dist
      }
    }
  }
  colnames(pval)  <- c(targetmetrics, paste("p-Value ", targetmetrics, sep=""))
  pvals <- list(ECDs=dists, pValues=pval)
  pvals
}
