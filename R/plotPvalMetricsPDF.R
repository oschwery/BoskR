#' Plot p-values on PDF for sets of tree metrics
#'
#' Creates plots of p-values on their corresponding probability density functio (PDF), based on sets of simulated and empirical distributions of tree metrics.
#'
#' @param empMetrics Metrics of empirical tree or set of trees; output of `GetTreeMetrics` or formatted the same way.
#' @param simMetrics Metrics of sets of simulated trees; output of `GetTreeMetrics` or formatted the same way.
#' @param set Numerical index for which of the sets of pairs of empirical and simulated metrics to be plotted; default NULL will plot all sets.
#' @param metricset String specifying which tree metrics to use; either "classic" (default) or "nodibranch" (or if Laplacian spectra are used also either "spectR" or "spectrRnorm"; but not "ALL" and "AllQuick", to reduce the number of plots crammed in one); for more information on the options see Details of `PvalMetrics()`.
#' @return An array of plots.
#'
#' @export

plotPvalMetricsPDF <- function(empMetrics, simMetrics, set=NULL, metricset="classic") {
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
    for (k in 1:length(empMetrics$metrics[, 1])) {
      if (plotcounter %% 4 == 0) {  # open new plot window after each 4 plots
        if(.Platform$OS.type=="windows") {  # make windows usable
          quartz<-function() windows()
        }
        quartz()
        par(mfrow=c(4, empno))
      }
      if (is.na(empMetrics$metrics[k])) {
        for (emp in 1:empno) {
          plot(1, type="n", axes=F, xlab="", ylab="")
        }
        plotcounter <- plotcounter+1
      } else {
        plotPvalsPDF(empMetrics, simMetrics, set=k, metricset, inloop=TRUE)
        plotcounter <- plotcounter+1
      }
    }
  } else {
    plotPvalsPDF(empMetrics, simMetrics, set=set, metricset, inloop=FALSE)
  }
}


#' Plot p-values for tree metrics
#'
#' Internal function to create plots of p-values based on simulated and empirical distributions of tree metrics
#'
#' @param empMetrics Metrics of empirical tree or set of trees; output of `GetTreeMetrics` or formatted the same way.
#' @param simMetrics Metrics of sets of simulated trees; output of `GetTreeMetrics` or formatted the same way.
#' @param set Numerical index for which of the sets of pairs of empirical and simulated metrics to be plotted.
#' @param metricset String specifying which tree metrics to use; either "classic" (default) or "nodibranch" (or if Laplacian spectra are used also either "spectR" or "spectrRnorm"; but not "ALL" and "AllQuick", to reduce the number of plots crammed in one); for more information on the options see Details of `PvalMetrics()`.
#' @param inloop Logical indicating whether the function is called from within a loop (TRUE) or not (FALSE).
#' @return An array of plots.
#'
#' @noRd

plotPvalsPDF <- function(empMetrics, simMetrics, set, metricset, inloop=FALSE) {
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
  for (j in 1:length(targetmetrics)) {
    dat <- na.omit(simMetrics[[set]]$metrics[, targetmetrics[j]])
    plot(density(as.numeric(dat)), xlim=c(min(c(dat, empMetrics[1]$metrics[set, targetmetrics[j]])), max((c(dat, empMetrics[1]$metrics[set, targetmetrics[j]])))), xlab=targetmetrics[j], ylab="Density", main=paste("PDF Tree", set, sep=" "))
    polygon(density(as.numeric(dat)), col="lightgray", border="darkgray")
    abline(v=empMetrics[1]$metrics[set, targetmetrics[j]], col="black", lwd=1.5)
  }
}
