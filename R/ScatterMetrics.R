#' 3D Metrics Scatterplot
#'
#' Plots empirical trees and their simulations in tree metric space using a 3D scatterplot.
#'
#'
#'
#' @param
#' @return 3D scatterplot of trees in metric space.
#'
#' @export
#'
#' @import scatterplot3d scatterplot3d
#' @import scales alpha


ScatterMetrics <- function(empMetrics, simMetrics, pair=1, skim=FALSE, combine=FALSE, colours=c("black", "red"), transparencyEmp=0.8 ,transparencySim=0.2, pttype=16, ptsize=1.5, plottitle=paste("Empirical vs. Simulated Metrics Set", pair, sep=" "), perspective=-230) {
  if (combine == TRUE) {
    simComb <- c()
    for (i in 1:length(simMetrics)) {
      simComb <- rbind(simComb, simMetrics[[i]]$metrics)
    }
    ScatterMetricsCombo <-
  }
  if (skim == TRUE) {
    for (i in 1:nrow(empMetrics$metrics)) {
      ScatterMetricsPair(empMetrics, simMetrics, pair=i, colours=colours, transparencyEmp=transparencyEmp, transparencySim=transparencySim, pttype=pttype, ptsize=ptsize, plottitle=plottitle, perspective=perspective)
      invisible(readline(prompt="Press [enter] to continue"))
    }
  } else {
    ScatterMetricsPair(empMetrics, simMetrics, pair=pair, colours=colours, transparencyEmp=transparencyEmp, transparencySim=transparencySim, pttype=pttype, ptsize=ptsize, plottitle=plottitle, perspective=perspective)
  }
}


ScatterMetricsPair <- function(empMetrics, simMetrics, pair=1, colours=c("black", "red"), transparencyEmp=0.8, transparencySim=0.2, pttype=16, ptsize=0.8, plottitle=paste("Empirical vs. Simulated Metrics Set", pair, sep=" "), perspective=-230) {
  scatterplot3d(c(empMetrics$metrics[pair, "Asymmetry_St"], simMetrics[[pair]]$metrics[, "Asymmetry_St"]), c(empMetrics$metrics[pair, "Peakedness_St"], simMetrics[[pair]]$metrics[, "Peakedness_St"]), c(empMetrics$metrics[pair, "Princ_Eigenv_St"], simMetrics[[pair]]$metrics[, "Princ_Eigenv_St"]), xlab="Skewness", ylab="Kurtosis", zlab="Princ. Eigenvalue", pch=pttype, cex.symbols=ptsize, highlight.3d=FALSE, color=c(alpha(colours[1], transparencyEmp), alpha(c(rep(colours[2], times=nrow(simMetrics[[pair]]$metrics))), transparencySim)), type="p", main=plottitle, angle=perspective)
}

ScatterMetricsCombo <- function(empMetrics, simComb, colours=c("black", "red"), transparencyEmp=0.8, transparencySim=0.2, pttype=16, ptsize=0.8, plottitle="Empirical vs. Simulated Metrics Sets Combined", perspective=-230) {
  scatterplot3d(c(empMetrics$metrics[, "Asymmetry_St"], simComb[, "Asymmetry_St"]), c(empMetrics$metrics[, "Peakedness_St"], simComb[, "Peakedness_St"]), c(empMetrics$metrics[, "Princ_Eigenv_St"], simComb[, "Princ_Eigenv_St"]), xlab="Skewness", ylab="Kurtosis", zlab="Princ. Eigenvalue", pch=pttype, cex.symbols=ptsize, highlight.3d=FALSE, color=c(alpha(rep(colours[1], times=nrow(empMetrics$metrics)), transparencyEmp), alpha(c(rep(colours[2], times=nrow(simComb$metrics))), transparencySim)), type="p", main=plottitle, angle=perspective)
}


library(scatterplot3d)
library(scales)
# flip axes
#scatterplot3d(fulldata[,2], fulldata[,3], fulldata[,1], xlab="Skewness", ylab="Kurtosis", zlab="Princ. Eigenvalue", pch=16, cex.symbols=0.5, highlight.3d=FALSE, color=fulltype, type="p", main="All Models", angle=-230)  # type="h" for vertical drop lines




quartz()
scatterplot3d(c(BDMetricsMerg[, "Asymmetry_St"], EmpMetrics$metrics[, "Asymmetry_St"]), c(BDMetricsMerg[, "Peakedness1_St"], EmpMetrics$metrics[, "Peakedness1_St"]), c(BDMetricsMerg[,"Princ_Eigenv_St"], EmpMetrics$metrics[, "Princ_Eigenv_St"]), xlab="Skewness", ylab="Kurtosis", zlab="Princ. Eigenvalue", pch=16, cex.symbols=0.8, highlight.3d=FALSE, color=c(alpha(BDMetricsMerg[, "Colour"], 0.2), EmpMetrics$metrics[, "Colour"]), type="p", main="Birth-Death", angle=-230)  # type="h" for vertical drop lines


#Emp alone


quartz()
scatterplot3d(c(EmpMetrics$metrics[, "Asymmetry_St"]), c(EmpMetrics$metrics[, "Peakedness1_St"]), c(EmpMetrics$metrics[, "Princ_Eigenv_St"]), xlab="Skewness", ylab="Kurtosis", zlab="Princ. Eigenvalue", pch=16, cex.symbols=0.8, highlight.3d=FALSE, color=c(EmpMetrics$metrics[, "Colour"]), type="p", main="Empirical", angle=-230)  # type="h" for vertical drop lines
