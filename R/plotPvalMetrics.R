#' Plot p-values for sets of tree metrics
#'
#' Creates plots of p-values based on sets of simulated and empirical distributions of tree metrics
#'
#' @param pmetrics Object with ECDs and p-values of empirical and simulated tree shapes, output of `PvalMetrics` or formatted the same way.
#' @param set Numerical index for which of the sets of pairs of empirical and simulated metrics to be plotted; default NULL will plot all sets.
#' @return An array of plots.
plotPvalMetrics <- function(pmetrics, set=NULL) {
  if (is.null(set)) {
    plotcounter <- 0
    for (j in 1:length(pmetrics$ECDs)) {
      if (plotcounter %% 4 == 0) {  # open new plot window after each 4 plots
        quartz()
        par(mfrow=c(4, 3))
      }
      if (is.na(pmetrics$ECDs[j])) {
        plot(1, type="n", axes=F, xlab="", ylab="")
        plot(1, type="n", axes=F, xlab="", ylab="")
        plot(1, type="n", axes=F, xlab="", ylab="")
        plotcounter <- plotcounter+1
      } else {
        plotPvals(pmetrics, set=j, inloop=TRUE)
        plotcounter <- plotcounter+1
      }
    }
  } else {
    plotPvals(pmetrics, set=set, inloop=FALSE)
  }
}


#' Plot p-values for tree metrics
#'
#' Internal function to create plots of p-values based on simulated and empirical distributions of tree metrics
#'
#' @param pmetrics Object with ECDs and p-values of empirical and simulated tree shapes, output of `PvalMetrics` or formatted the same way.
#' @param set Numerical index for which of the sets of pairs of empirical and simulated metrics to be plotted.
#' @param inloop Logical indicating whether the function is called from within a loop (TRUE) or not (FALSE).
#' @return An array of plots.
plotPvals <- function(pmetrics, set, inloop=FALSE) {
  if (inloop == FALSE) {
    par(mfrow=c(length(set), ncol(pmetrics$pValues)/2))
  }
  for (i in 1:(ncol(pmetrics$pValues)/2)) {
    plot(pmetrics$ECDs[[set]][[i]], xlab= colnames(pmetrics$pValues)[i])
    abline(v=pmetrics$pValues[set, i])
    abline(h=pmetrics$pValues[set, i+(ncol(pmetrics$pValues)/2)])
  }
}
