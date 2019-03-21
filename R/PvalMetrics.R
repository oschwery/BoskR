#' Get p-values for tree metrics
#'
#' Estimates p-values based on simulated and empirical distributions of tree metrics
#'
#' The function uses an Empirical Cumulative Distribution function to determine the area under the curve of the metric values of the simulated trees, to get to a p-value for the position of the metrics of the empirical tree on that distribution.
#'
#' @param empMetrics Metrics of empirical tree or set of trees; output of `GetTreeMetrics` or formatted the same way.
#' @param simMetrics Metrics of sets of simulated trees; output of `GetTreeMetrics` or formatted the same way.
#' @param empirical_start Indicator whether empMetrics is based on empirical or simulated initial trees, default is TRUE (=empirical); mainly important for data format reasons.
#' @param methodnr Integral specifying which method is used: 1: BD, 2: TimeD-BD, 3: DD.
#' @return A list with two entries: `ECDs` is a list of Empirical Cumulative Distributions; `pValues` is a matrix with p-values for the targeted metrics
#'
#' @export

PvalMetrics <- function(empMetrics, simMetrics, empirical_start=TRUE, methodnr) {
  # loop getting distributions for all sim trees
  dists <- list()  # primer for ECDs
  targetmetrics <- c("Princ_Eigenv_St", "Asymmetry_St", "Peakedness_St")  # CHANGE once you have likelihood implemented!
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
  colnames(pval)  <- c(targetmetrics, paste("p-Value ", targetmetrics, sep=""))
  pvals <- list(ECDs=dists, pValues=pval)
  pvals
}
