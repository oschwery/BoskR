#' Plot p-values on PDF for sets of tree metrics
#'
#' Creates plots of p-values on their corresponding probability density function, based on sets of simulated and empirical distributions of tree metrics,
#'
#' @param empMetrics Metrics of empirical tree or set of trees; output of `GetTreeMetrics` or formatted the same way.
#' @param simMetrics Metrics of sets of simulated trees; output of `GetTreeMetrics` or formatted the same way.
#' @param set Numerical index for which of the sets of pairs of empirical and simulated metrics to be plotted; default NULL will plot all sets.
#' @return An array of plots.
plotPvalMetricsPDF <- function(empMetrics, simMetrics, set=NULL) {
  if (is.null(set)) {
    plotcounter <- 0
    for (k in 1:length(empMetrics$metrics[, 1])) {
      if (plotcounter %% 4 == 0) {  # open new plot window after each 4 plots
        quartz()
        par(mfrow=c(4, 3))
      }
      if (is.na(empMetrics$metrics[k])) {
        plot(1, type="n", axes=F, xlab="", ylab="")
        plot(1, type="n", axes=F, xlab="", ylab="")
        plot(1, type="n", axes=F, xlab="", ylab="")
        plotcounter <- plotcounter+1
      } else {
        plotPvalsPDF(empMetrics, simMetrics, set=k, inloop=TRUE)
        plotcounter <- plotcounter+1
      }
    }
  } else {
    plotPvalsPDF(empMetrics, simMetrics, set=set, inloop=FALSE)
  }
}


#' Plot p-values for tree metrics
#'
#' Internal function to create plots of p-values based on simulated and empirical distributions of tree metrics
#'
#' @param empMetrics Metrics of empirical tree or set of trees; output of `GetTreeMetrics` or formatted the same way.
#' @param simMetrics Metrics of sets of simulated trees; output of `GetTreeMetrics` or formatted the same way.
#' @param set Numerical index for which of the sets of pairs of empirical and simulated metrics to be plotted.
#' @param inloop Logical indicating whether the function is called from within a loop (TRUE) or not (FALSE).
#' @return An array of plots.
plotPvalsPDF <- function(empMetrics, simMetrics, set, inloop=FALSE) {
  targetmetrics <- c("Princ_Eigenv_St", "Asymmetry_St", "Peakedness_St")
  if (inloop == FALSE) {
    par(mfrow=c(length(set), length(targetmetrics)))
  }
  for (j in 1:length(targetmetrics)) {
    plot(density(as.numeric(simMetrics[[set]]$metrics[, targetmetrics[j]])), xlim=c(min(c(simMetrics[[set]]$metrics[, targetmetrics[j]], empMetrics[1]$metrics[set, targetmetrics[j]])), max((c(simMetrics[[set]]$metrics[, targetmetrics[j]], empMetrics[1]$metrics[set, targetmetrics[j]])))), xlab=targetmetrics[j], ylab="Density", main=paste("PDF Tree", set, sep=" "))
    polygon(density(as.numeric(simMetrics[[set]]$metrics[, targetmetrics[j]])), col="lightgray", border="darkgray")
    abline(v=empMetrics[1]$metrics[set, targetmetrics[j]], col="black", lwd=1.5)
  }
}
