% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotPvalMetricsPDF.R
\name{plotPvalMetricsPDF}
\alias{plotPvalMetricsPDF}
\title{Plot p-values on PDF for sets of tree metrics}
\usage{
plotPvalMetricsPDF(empMetrics, simMetrics, set = NULL)
}
\arguments{
\item{empMetrics}{Metrics of empirical tree or set of trees; output of \code{GetTreeMetrics} or formatted the same way.}

\item{simMetrics}{Metrics of sets of simulated trees; output of \code{GetTreeMetrics} or formatted the same way.}

\item{set}{Numerical index for which of the sets of pairs of empirical and simulated metrics to be plotted; default NULL will plot all sets.}
}
\value{
An array of plots.
}
\description{
Creates plots of p-values on their corresponding probability density function, based on sets of simulated and empirical distributions of tree metrics,
}