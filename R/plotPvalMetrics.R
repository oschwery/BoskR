# plots for p-pValues
plotPvalMetrics <- function(pmetrics, set) {
  par(mfrow=c(length(set), ncol(pmetrics$pValues)/2))
  for (i in 1:(ncol(pmetrics$pValues)/2)) {
    plot(pmetrics$ECDs[[set]][[i]], xlab= colnames(pmetrics$pValues)[i])
    abline(v=pmetrics$pValues[set, i])
    abline(h=pmetrics$pValues[set, i+(ncol(pmetrics$pValues)/2)])
  }
}
