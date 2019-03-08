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

ScatterMetrics <- function(empMetrics, simMetrics) {
  scatterplot3d(c())
}


library(scatterplot3d)
library(scales)
# flip axes
#scatterplot3d(fulldata[,2], fulldata[,3], fulldata[,1], xlab="Skewness", ylab="Kurtosis", zlab="Princ. Eigenvalue", pch=16, cex.symbols=0.5, highlight.3d=FALSE, color=fulltype, type="p", main="All Models", angle=-230)  # type="h" for vertical drop lines

#Emp alone


quartz()
scatterplot3d(c(EmpMetrics$metrics[, "Asymmetry_St"]), c(EmpMetrics$metrics[, "Peakedness1_St"]), c(EmpMetrics$metrics[, "Princ_Eigenv_St"]), xlab="Skewness", ylab="Kurtosis", zlab="Princ. Eigenvalue", pch=16, cex.symbols=0.8, highlight.3d=FALSE, color=c(EmpMetrics$metrics[, "Colour"]), type="p", main="Empirical", angle=-230)  # type="h" for vertical drop lines


quartz()
scatterplot3d(c(BDMetricsMerg[, "Asymmetry_St"], EmpMetrics$metrics[, "Asymmetry_St"]), c(BDMetricsMerg[, "Peakedness1_St"], EmpMetrics$metrics[, "Peakedness1_St"]), c(BDMetricsMerg[,"Princ_Eigenv_St"], EmpMetrics$metrics[, "Princ_Eigenv_St"]), xlab="Skewness", ylab="Kurtosis", zlab="Princ. Eigenvalue", pch=16, cex.symbols=0.8, highlight.3d=FALSE, color=c(alpha(BDMetricsMerg[, "Colour"], 0.2), EmpMetrics$metrics[, "Colour"]), type="p", main="Birth-Death", angle=-230)  # type="h" for vertical drop lines
