\name{PvalMetrics}
\alias{PvalMetrics}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
PvalMetrics(empMetrics, simMetrics)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{empMetrics}{
%%     ~~Describe \code{empMetrics} here~~
}
  \item{simMetrics}{
%%     ~~Describe \code{simMetrics} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (empMetrics, simMetrics) 
{
    dists <- list()
    targetmetrics <- c("Princ_Eigenv_St", "Asymmetry_St", "Peakedness1_St")
    pval <- matrix(nrow = nrow(empMetrics$metrics), ncol = 2 * 
        length(targetmetrics))
    for (j in 1:nrow(empMetrics$metrics)) {
        dist <- list()
        for (i in 1:length(targetmetrics)) {
            dist[[i]] <- ecdf(simMetrics[[j]]$metrics[, targetmetrics[i]])
            pval[j, i] <- empMetrics[1]$metrics[j, targetmetrics[i]]
            pval[j, i + length(targetmetrics)] <- dist[[i]](c(empMetrics[1]$metrics[j, 
                targetmetrics[i]]))
        }
        names(dist) <- targetmetrics
        dists[[j]] <- dist
    }
    colnames(pval) <- c(targetmetrics, paste("p-Value ", targetmetrics, 
        sep = ""))
    pvals <- list(ECDs = dists, pValues = pval)
    pvals
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
