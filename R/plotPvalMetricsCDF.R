#' Plot p-values on CDF for sets of tree metrics
#'
#' Creates plots of p-values on their corresponding cumulative distribution function, based on sets of simulated and empirical distributions of tree metrics,
#'
#' @param pmetrics Object with ECDs and p-values of empirical and simulated tree shapes, output of `PvalMetrics` or formatted the same way.
#' @param set Numerical index for which of the sets of pairs of empirical and simulated metrics to be plotted; default NULL will plot all sets.
#' @param metricset String specifying which tree metrics to use; default is "spectR", other options are "spectrRnorm", "classic", and "nodibranch" (but not "ALL" and "AllQuick", to reduce the number of plots crammed in one); for more information on the options see Details of `PvalMetrics()`.
#' @return An array of plots.
#'
#' @export

plotPvalMetricsCDF <- function(pmetrics, set=NULL, metricset="spectR") {
  if (is.null(set)) {
    plotcounter <- 0
       empno <- c()
    if (metricset == "spectR") {
      empno = 4
    } else if (metricset == "spectRnorm") {
      empno = 3
    } else if (metricset == "classic") {
      empno = 6
    } else if (metricset == "nodibranch") {
      empno = 7
    }
    for (k in 1:length(pmetrics$ECDs)) {
      if (plotcounter %% 4 == 0) {  # open new plot window after each 4 plots
        if(.Platform$OS.type=="windows") {  # make windows usable
          quartz<-function() windows()
        }
        quartz()
        par(mfrow=c(4, empno))
      }
      if (is.na(pmetrics$ECDs[k])) {
        for (mets in 1:empno) {
        plot(1, type="n", axes=F, xlab="", ylab="")
        }
        plotcounter <- plotcounter+1
      } else {
        plotPvalsCDF(pmetrics, set=k, metricset, inloop=TRUE)
        plotcounter <- plotcounter+1
      }
    }
  } else {
    plotPvalsCDF(pmetrics, set=set, metricset, inloop=FALSE)
  }
}


#' Plot p-values for tree metrics
#'
#' Internal function to create plots of p-values based on simulated and empirical distributions of tree metrics
#'
#' @param pmetrics Object with ECDs and p-values of empirical and simulated tree shapes, output of `PvalMetrics` or formatted the same way.
#' @param set Numerical index for which of the sets of pairs of empirical and simulated metrics to be plotted.
#' @param metricset String specifying which tree metrics to use; default is "spectR", other options are "spectrRnorm", "classic", and "nodibranch" (but not "ALL" and "AllQuick", to reduce the number of plots crammed in one); for more information on the options see Details of `PvalMetrics()`.
#' @param inloop Logical indicating whether the function is called from within a loop (TRUE) or not (FALSE).
#' @return An array of plots.
#'
#' @noRd

plotPvalsCDF <- function(pmetrics, set, metricset, inloop=FALSE) {
  targetmetrics <-c()
  if (metricset == "spectR") {
    targetmetrics <- c("Princ_Eigenv_St", "Asymmetry_St", "Peakedness_St", "Eigengap_St")  # CHANGE once you have likelihood implemented!
  } else if (metricset == "spectRnorm") {
    targetmetrics <- c("Princ_Eigenv_Nor", "Asymmetry_Nor", "Peakedness_Nor")
  } else if (metricset == "classic") {
    targetmetrics <- c("Colless", "Sackin", "Cherries", "pitchforks", "AvgLadder", "Gamma")
  } else if (metricset == "nodibranch") {
    targetmetrics <- c("N_tax", "Min_NodeAge", "Median_NodeAge", "Max_NodeAge", "Min_BranchLength", "Median_BranchLength", "Max_BranchLength")
  }
  if (inloop == FALSE) {
    par(mfrow=c(length(set), length(targetmetrics)))
  }
  for (i in 1:length(targetmetrics)) {
    plot(pmetrics$ECDs[[set]][[targetmetrics[i]]], xlab=targetmetrics[i], main=paste("CDF Tree", set, sep=" "))
    abline(v=pmetrics$pValues[set, targetmetrics[i]])  # vertical line for the values of the empirical metric
    I <- which(colnames(pvalues$pValues)==targetmetrics[i])
    # Horizontal lines for the two tails of pvalue (i.e. where the vertical line meets the CDF)
    abline(h=pmetrics$pValues[set, I+(ncol(pmetrics$pValues)/2)]/2)
    abline(h=1-(pmetrics$pValues[set, I+(ncol(pmetrics$pValues)/2)]/2))
  }
}
