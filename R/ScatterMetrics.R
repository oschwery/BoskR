#' 3D Metrics Scatterplot
#'
#' Plots empirical trees and their simulations in tree metric space using a 3D scatterplot.
#'
#' The function uses the internals `ScatterMetricsPair` and `ScatterMetricsCombo` and plots the empirical input-trees and their corresponding simulations in the metric space (asymmetry x peakedness x principal Eigenvalue) as a 3D scatterplot. It allows to either plot them all combined, or pairwise. The latter meaning each empirical tree is plotted with its corresponding simulations only, either one at a time or all together interactively (one advances through the plots by pressing enter). The basic function used is `scatterplot3d`, from the package with the same name.
#'
#' @param empMetrics Metrics of empirical tree or set of trees; output of `GetTreeMetrics` or formatted the same way.
#' @param simMetrics Metrics of sets of simulated trees; output of `GetTreeMetrics` or formatted the same way.
#' @param pair Numerical index for which of the sets of pairs of empirical and simulated metrics to be plotted. Value is ignored if `skim` or `combine` are `TRUE`, or `sets` is larger than 1 (the default).
#' @param sets Integer for number of different subsets within `empMetrics`, e.g., sets originating from different models; is ignored if 1 (default).
#' @param setsize Size of each of the sets. If `NA` (default), set size will be determined automatically by dividing the metrics samples by the number of sets; single numbers get recycled for each set, and a vector of sizes can be provided if the sets all differ.
#' @param skim Logical, creates interactive plot of all pairs of empirical trees and their simulations if `TRUE`; one can advance through the plots by hitting enter.
#' @param separate Boolean for whether sets should be plotted in same plot or separate plots (default is `FALSE`); does not automatically open new plot devices or panels.
#' @param setnames Names for each set, to be used to label plots when `separate=TRUE`
#' @param combine Logical, combines all empirical and simulated trees into one plot if `TRUE`.
#' @param colours Vector of length two, indicating the desired colours for empirical trees and simulated treees, in that order (defaults are "black" and "red", respectively); or vector ndicating the desired colours for each set, in same order as sets are in metrics object.
#' @param transparencyEmp Value determining the transparency of the empirical tree plot points (0: completely transparent, 1: completely opaque; corresponding to `alpha` from package `scales`).
#' @param transparencySim Value determining the transparency of the simulated tree plot points (0: completely transparent, 1: completely opaque; corresponding to `alpha` from package `scales`).
#' @param pch Shape of plot symbols; default 16.
#' @param cex.symbols Size of plot symbols; default 1.5.
#' @param main String for plot title; default "Empirical vs. Simulated Metrics Set", followed by pair number plotted, or "Combined", or "Metrics Sets" preceeded by setname.
#' @param angle Rotation of the plot, determined by angle between x and y axis (corresponding to `scatterplot3d`); default -230.
#' @param xlim Min and max values of x-axis, default is `NULL` (fits to data automatically), can be manually adjusted as `c(min, max)`; same for `ylim` and `zlim`.
#' @return 3D scatterplot of trees in metric space, or a series of such plots to skip through.
#'
#' @export
#'
#' @importFrom scatterplot3d scatterplot3d
#' @importFrom scales alpha

ScatterMetrics <- function(empMetrics, simMetrics=NULL, pair=1, sets=1, setsize=NA, skim=FALSE, separate=FALSE, setnames=NA, combine=FALSE, colours=c("black", "red"), transparencyEmp=0.8, transparencySim=0.2, pch=16, cex.symbols=1.5, main=paste("Empirical vs. Simulated Metrics Set", pair, sep=" "), angle=-230, xlim=NULL, ylim=NULL, zlim=NULL) {
  if (combine == TRUE) {
    simComb <- c()
    for (i in 1:length(simMetrics)) {
      simComb <- rbind(simComb, simMetrics[[i]]$metrics)
    }
    ScatterMetricsCombo(empMetrics, simComb, colours=colours, transparencyEmp=transparencyEmp, transparencySim=transparencySim, pch=pch, cex.symbols=cex.symbols, main=main, angle=angle)
  } else if (sets > 1) {
    ScatterMetricsSets(empMetrics, sets=sets, setsize=setsize, separate=separate, setnames=setnames, colours=colours, transparencyEmp=transparencyEmp, pch=pch, cex.symbols=cex.symbols, main="Metrics Set", angle=-230, xlim=NULL, ylim=NULL, zlim=NULL) {

  } else if (skim == TRUE) {
    for (i in 1:nrow(empMetrics$metrics)) {
      ScatterMetricsPair(empMetrics, simMetrics, pair=i, colours=colours, transparencyEmp=transparencyEmp, transparencySim=transparencySim, pch=pch, cex.symbols=cex.symbols, main=main, angle=angle, xlim=xlim, ylim=ylim, zlim=zlim)
      invisible(readline(prompt="Press [enter] to continue"))
    }
  } else {
    ScatterMetricsPair(empMetrics, simMetrics, pair=pair, colours=colours, transparencyEmp=transparencyEmp, transparencySim=transparencySim, pch=pch, cex.symbols=cex.symbols, main=main, angle=angle, xlim=xlim, ylim=ylim, zlim=zlim)
  }
}

#' 3D Metrics Scatterplot of Tree Set
#'
#' Plots an empirical tree and its simulations in tree metric space using a 3D scatterplot.
#'
#' Internal function for `ScatterMetrics`, plots an empirical input-tree and its corresponding simulations in the metric space (asymmetry x peakedness x principal Eigenvalue) as a 3D scatterplot. Each empirical tree is plotted with its corresponding simulations only. The basic function used is `scatterplot3d`, from the package with the same name.
#'
#' @param empMetrics Metrics of empirical tree or set of trees; output of `GetTreeMetrics` or formatted the same way.
#' @param simMetrics Metrics of sets of simulated trees; output of `GetTreeMetrics` or formatted the same way.
#' @param pair Numerical index for which of the sets of pairs of empirical and simulated metrics to be plotted.
#' @param colours Vector of length two, indicating the desired colours for empirical trees and simulated treees, in that order (defaults are "black" and "red", respectively).
#' @param transparencyEmp Value determining the transparency of the empirical tree plot points (0: completely transparent, 1: completely opaque; corresponding to `alpha` from package `scales`).
#' @param transparencySim Value determining the transparency of the simulated tree plot points (0: completely transparent, 1: completely opaque; corresponding to `alpha` from package `scales`).
#' @param pch Shape of plot symbols; default 16.
#' @param cex.symbols Size of plot symbols; default 1.5.
#' @param main String for plot title; default "Empirical vs. Simulated Metrics Set", followed by pair number plotted.
#' @param angle Rotation of the plot, determined by angle between x and y axis (corresponding to `scatterplot3d`); default -230.
#' @param xlim Min and max values of x-axis, default is `NULL` (fits to data automatically), can be manually adjusted as `c(min, max)`; same for `ylim` and `zlim`.
#' @return 3D scatterplot of trees in metric space.
#'
#' @noRd
#'
#' @importFrom scatterplot3d scatterplot3d
#' @importFrom scales alpha

ScatterMetricsPair <- function(empMetrics, simMetrics, pair=1, colours=c("black", "red"), transparencyEmp=0.8, transparencySim=0.2, pch=16, cex.symbols=0.8, main=paste("Empirical vs. Simulated Metrics Set", pair, sep=" "), angle=-230, xlim=NULL, ylim=NULL, zlim=NULL) {
  scatterplot3d(c(empMetrics$metrics[pair, "Asymmetry_St"], simMetrics[[pair]]$metrics[, "Asymmetry_St"]), c(empMetrics$metrics[pair, "Peakedness_St"], simMetrics[[pair]]$metrics[, "Peakedness_St"]), c(empMetrics$metrics[pair, "Princ_Eigenv_St"], simMetrics[[pair]]$metrics[, "Princ_Eigenv_St"]), xlab="Skewness", ylab="Kurtosis", zlab="Princ. Eigenvalue", pch=pch, cex.symbols=cex.symbols, highlight.3d=FALSE, color=c(alpha(colours[1], transparencyEmp), alpha(c(rep(colours[2], times=nrow(simMetrics[[pair]]$metrics))), transparencySim)), type="p", main=main, angle=angle, xlim=xlim, ylim=ylim, zlim=zlim)
}

#' 3D Metrics Scatterplot of Combined Tree Sets
#'
#' Plots the combined set of empirical trees and their simulations in tree metric space using a 3D scatterplot.
#'
#' Internal function for `ScatterMetrics`, plots all empirical input-trees and their corresponding simulations in the metric space (asymmetry x peakedness x principal Eigenvalue) as a single 3D scatterplot. The basic function used is `scatterplot3d`, from the package with the same name.
#'
#' @param empMetrics Metrics of empirical tree or set of trees; output of `GetTreeMetrics` or formatted the same way.
#' @param simComb Matrix of metrics of sets of simulated trees; `ScatterMetrics` will create this by combining the output of `GetTreeMetrics` or data formatted the same way.
#' @param colours Vector of length two, indicating the desired colours for empirical trees and simulated treees, in that order (defaults are "black" and "red", respectively).
#' @param transparencyEmp Value determining the transparency of the empirical tree plot points (0: completely transparent, 1: completely opaque; corresponding to `alpha` from package `scales`).
#' @param transparencySim Value determining the transparency of the simulated tree plot points (0: completely transparent, 1: completely opaque; corresponding to `alpha` from package `scales`).
#' @param pch Shape of plot symbols; default 16.
#' @param cex.symbols Size of plot symbols; default 1.5.
#' @param main String for plot title; default "Empirical vs. Simulated Metrics Sets Combined".
#' @param angle Rotation of the plot, determined by angle between x and y axis (corresponding to `scatterplot3d`); default -230.
#' @param xlim Min and max values of x-axis, default is `NULL` (fits to data automatically), can be manually adjusted as `c(min, max)`; same for `ylim` and `zlim`.
#' @return 3D scatterplot of all sets of trees in metric space.
#'
#' @noRd
#'
#' @importFrom scatterplot3d scatterplot3d
#' @importFrom scales alpha

ScatterMetricsCombo <- function(empMetrics, simComb, colours=c("black", "red"), transparencyEmp=0.8, transparencySim=0.2, pch=16, cex.symbols=0.8, main="Empirical vs. Simulated Metrics Sets Combined", angle=-230, xlim=NULL, ylim=NULL, zlim=NULL) {
  scatterplot3d(c(empMetrics$metrics[, "Asymmetry_St"], simComb[, "Asymmetry_St"]), c(empMetrics$metrics[, "Peakedness_St"], simComb[, "Peakedness_St"]), c(empMetrics$metrics[, "Princ_Eigenv_St"], simComb[, "Princ_Eigenv_St"]), xlab="Skewness", ylab="Kurtosis", zlab="Princ. Eigenvalue", pch=pch, cex.symbols=cex.symbols, highlight.3d=FALSE, color=c(alpha(rep(colours[1], times=nrow(empMetrics$metrics)), transparencyEmp), alpha(c(rep(colours[2], times=nrow(simComb))), transparencySim)), type="p", main=main, angle=angle, xlim=xlim, ylim=ylim, zlim=zlim)
}

#' 3D Metrics Scatterplot of Combined Tree Sets
#'
#' Plots the combined set of empirical trees and their simulations in tree metric space using a 3D scatterplot.
#'
#' Internal function for `ScatterMetrics`, plots all empirical input-trees and their corresponding simulations in the metric space (asymmetry x peakedness x principal Eigenvalue) as a single 3D scatterplot. The basic function used is `scatterplot3d`, from the package with the same name.
#'
#' @param empMetrics Metrics of empirical tree or set of trees; output of `GetTreeMetrics` or formatted the same way.
#' @param sets Integer for number of different subsets within `empMetrics`, e.g., sets originating from different models.
#' @param setsize Size of each of the sets. If `NA` (default), set size will be determined automatically by dividing the metrics samples by the number of sets; single numbers get recycled for each set, and a vector of sizes can be provided if the sets all differ.
#' @param separate Boolean for whether sets should be plotted in same plot or separate plots (default is `FALSE`); does not automatically open new plot devices or panels.
#' @param setnames Names for each set, to be used to label plots when `separate=TRUE`
#' @param colours Vector indicating the desired colours for each set, in same order as sets are in metrics object.
#' @param transparencyEmp Value determining the transparency of the tree plot points (0: completely transparent, 1: completely opaque; corresponding to `alpha` from package `scales`).
#' @param pch Shape of plot symbols; default 16.
#' @param cex.symbols Size of plot symbols; default 1.5.
#' @param main String for plot title; default "Metrics Sets".
#' @param angle Rotation of the plot, determined by angle between x and y axis (corresponding to `scatterplot3d`); default -230.
#' @param xlim Min and max values of x-axis, default is `NULL` (fits to data automatically), can be manually adjusted as `c(min, max)`; same for `ylim` and `zlim`.
#' @return 3D scatterplot of all sets of trees in metric space.
#'
#' @noRd
#'
#' @importFrom scatterplot3d scatterplot3d
#' @importFrom scales alpha

ScatterMetricsSets <- function(empMetrics, sets=1, setsize=NA, separate=FALSE, setnames=NA, colours=c("black", "red"), transparencyEmp=0.8, pch=16, cex.symbols=0.8, main="Metrics Set", angle=-230, xlim=NULL, ylim=NULL, zlim=NULL) {
    if (is.na(setsize)) {
        setsize <- rep((nrow(empMetrics$metrics)/sets), times=sets)
    } else if (length(setsize) != sets) {
        setsize <- rep(setsize, times=sets)
    }
    plotcols <- c()
    for (set in 1:sets){
        plotcols <- c(plotcols, rep(colours[set], times=setsize[set]))
    }
    if (separate == FALSE) {
    scatterplot3d(
        empMetrics$metrics[, "Asymmetry_St"], 
        empMetrics$metrics[, "Peakedness_St"], 
        empMetrics$metrics[, "Princ_Eigenv_St"],
        xlab="Skewness", ylab="Kurtosis", zlab="Princ. Eigenvalue", pch=pch, cex.symbols=cex.symbols, highlight.3d=FALSE, color=alpha(plotcols, transparencyEmp), type="p", main=main, angle=angle, xlim=xlim, ylim=ylim, zlim=zlim)
    } else if (separate == TRUE) {
        stopset <- 0
        for (sepset in 1:sets) {
            startset <- 1 + stopset
            stopset <- stopset + setsize[sepset]
            # print(paste(models[sepset], startset, stopset, plotcols[startset], plotcols[stopset], sep=" "))
            scatterplot3d(
                empMetrics$metrics[c(startset:stopset), "Asymmetry_St"], 
                empMetrics$metrics[c(startset:stopset), "Peakedness_St"], 
                empMetrics$metrics[c(startset:stopset), "Princ_Eigenv_St"],
                xlab="Skewness", ylab="Kurtosis", zlab="Princ. Eigenvalue", pch=pch, cex.symbols=cex.symbols, highlight.3d=FALSE, 
                color=alpha(plotcols[startset:stopset], transparencyEmp), 
                type="p", main=paste(setnames[sepset], main, sep=" "), angle=angle, xlim=xlim, ylim=ylim, zlim=zlim)
        }
    }
}
