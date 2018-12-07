#' Get p-values for tree metrics
#'
#' Estimates p-values based on simulated and empirical distributions of tree metrics
#'
#' The function uses an Empirical Cumulative Distribution function to determine the area under the curve of the metric values of the simulated trees, to get to a p-value for the position of the metrics of the empirical tree on that distribution.
#'
#' @param empMetrics Metrics of empirical tree or set of trees; output of `GetTreeMetrics` or formatted the same way.
#' @param simMetrics Metrics of sets of simulated trees; output of `GetTreeMetrics` or formatted the same way.
#' @param current_case Indicator whether empMetrics is based of simulated initial trees, in which case the value should be "simed"; mainly important for data format reasons.
#' @param methodnr Integral specifying which method is used: 1: BD, 2: TimeD-BD, 3: DD.
#' @return A list with two entries: `ECDs` is a list of Empirical Cumulative Distributions; `pValues` is a matrix with p-values for the targeted metrics
PvalMetrics <- function(empMetrics, simMetrics, current_case, methodnr) {
  # loop getting distributions for all sim trees
  dists <- list()
  targetmetrics <- c("Princ_Eigenv_St", "Asymmetry_St", "Peakedness1_St")  # CHANGE once you have likelihood implemented!
  Methods <- c("BD", "TimeD-BD", "DD")
  current_method <- Methods[methodnr]
  pval <- matrix(nrow=nrow(empMetrics$metrics), ncol=2*length(targetmetrics))
  if (current_case == "simed") {
    k <- 1
    for (j in 1:nrow(empMetrics$metrics)) {
      dist <- list()
      if (paste("metrics", rownames(empMetrics$metrics)[j], current_method, sep="_") == names(simMetrics[k])) {
        for (i in 1:length(targetmetrics)) {
          deP <- c()
          dist[[i]] <- ecdf(simMetrics[[k]]$metrics[, targetmetrics[i]])
          pval[j, i] <- empMetrics[1]$metrics[j, targetmetrics[i]]
          deP <- dist[[i]](c(empMetrics[1]$metrics[j, targetmetrics[i]]))
          pval[j, i+length(targetmetrics)] <- (min(c(deP, 1-deP)))/2
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
      dist <- list()
      for (i in 1:length(targetmetrics)) {
        deP <- c()
        dist[[i]] <- ecdf(simMetrics[[j]]$metrics[, targetmetrics[i]])
        pval[j, i] <- empMetrics[1]$metrics[j, targetmetrics[i]]
        deP <- dist[[i]](c(empMetrics[1]$metrics[j, targetmetrics[i]]))
        pval[j, i+length(targetmetrics)] <- (min(c(deP, 1-deP)))/2
      }
      names(dist) <- targetmetrics
      dists[[j]] <- dist
    }
  }
  colnames(pval)  <- c(targetmetrics, paste("p-Value ", targetmetrics, sep=""))
  pvals <- list(ECDs=dists, pValues=pval)
  pvals
}
