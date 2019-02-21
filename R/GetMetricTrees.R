#' Simulate trees based on empirical estimations or set parameters
#'
#' `GetMetricTrees` simulates trees under a given model based on either parameter estimates from empirical trees or pre-set parameters.
#'
#' The function will simulate a number of trees based on either the parameters inferred from one or several empirical trees (given through `empParams` if `empirical_start=TRUE`), or user-specified parameters (if `empirical_start=FALSE`)
#'
#' @param trset Integer indicating the tree set to be evaluated
#' @param empirical_start `TRUE` to use parameters estimated from empirical trees, `FALSE` to use user-specified ones
#' @param current_case Integer indicating which tree set tree belongs to (used for labeling purposes).
#' @param empParams Nested list object with tree parameters as inferred through `GetParams` from one or several empirical trees
#' @param current_method Method to be used for simulation, either `"BD", "TimeD-BD", "DD", "CD", "TraitD"` for birth-death, time-dependent birth-death, diversity dependent, clade dependent, or trait dependent diversification respectively.
#' @param N Number of taxa
#' @param Numbsim1 Number of trees to simulate per each
#' @param Lambda Speciation rate
#' @param Mu Extinction rate
#' @param l Speciation rate
#' @param a Extinction fracion (Mu/Lambda)
#' @param LambdaFun Function for speciation rate
#' @param MuFun Function for extinction rate
#' @param TreeAge Stem age of tree
#' @param BiSSEpars Parameters from BiSSE
#' @param tree Phylogeny
#' @return A list of trees of class multiPhylo
GetMetricTrees <- function(trset=trset, empirical_start=FALSE, current_case=NULL, empParams=empParams, current_method, N=NULL, Numbsim1, Lambda, Mu, l=NULL, a=NULL, LambdaFun=NULL, MuFun=NULL, TreeAge=NULL, BiSSEpars=NULL, tree=NULL) {
  # use empirical trees (or not)
  if (empirical_start == TRUE) {
    N <- empParams[[trset]]$N
    Lambda <- empParams[[trset]]$Lambda
    Mu <- empParams[[trset]]$Mu
    K <- empParams[[trset]]$K
    l <- empParams[[trset]]$l
    a <- empParams[[trset]]$a
    LambdaFun <- empParams[[trset]]$LambdaFun
    MuFun <- empParams[[trset]]$MuFun
    TreeAge <-empParams[[trset]]$TreeAge
    BiSSEpars <- empParams[[trset]]$BiSSEpars
  }
  if (is.null(K)) {
    K <- N
  }
  if (is.null(l)) {
    l <- Lambda
  }
  if (length(a) == 0) {
    a <- Mu/Lambda
  }
  if (is.null(LambdaFun)) {
    LambdaFun <- function(t) l * exp(a*t)
  }
  if (is.null(MuFun)) {
    MuFun <- Mu  #as.numeric(empirical_solution[1,]$mu0)
  }
  if (is.null(TreeAge) & !is.null(tree)) {
    TreeAge <- max(branching.times(tree))
  }
  if (is.null(TreeAge) & is.null(tree)) {
    TreeAge <- (N * 0.75)
  }
  if (is.null(BiSSEpars) & current_method == "TraitD") {  # this part doesn't make that much sense at this point...
    BiSSEpars <- c(median(as.numeric(empirical_solution$lambda0), na.rm=TRUE), median(as.numeric(empirical_solution$lambda1), na.rm=TRUE), median(as.numeric(empirical_solution$mu0), na.rm=TRUE), median(as.numeric(empirical_solution$mu1), na.rm=TRUE), median(as.numeric(empirical_solution$q01), na.rm=TRUE), median(as.numeric(empirical_solution$q10), na.rm=TRUE))
  }
# display used metrics
print(trset)
if (empirical_start==TRUE) {
  print(names(emptrees[trset]))
}
print(paste("N",N, sep=": "))
print(paste("Numbsim1", Numbsim1, sep=": "))
print(paste("Lambda", Lambda, sep=": "))
print(paste("Mu", Mu, sep=": "))
print(paste("K", K, sep=": "))
print(paste("l", l, sep=": "))
print(paste("a", a, sep=": "))
print(paste("MuFun", MuFun, sep=": "))
print(paste("TreeAge", TreeAge, sep=": "))
print(paste("BiSSEpars", BiSSEpars, sep=": "))
##############
# pick right approach
  if (is.na(N)) {
    trees <- NA
  } else if (!is.na(N)) {
    if (current_method == "BD") {
      trees <- try(sim.bd.age(age=TreeAge, numbsim=Numbsim1, lambda=Lambda, mu=Mu, frac=1, mrca=TRUE, complete=FALSE, K=0), FALSE)
    } else if (current_method == "TimeD-BD") {
      trees <- try(tess.sim.age(n=Numbsim1, age=TreeAge, lambda=LambdaFun, mu=MuFun, massExtinctionTimes = c(), massExtinctionSurvivalProbabilities = c(), samplingProbability = 1, samplingStrategy = "uniform", maxTaxa = Inf, MRCA = TRUE), FALSE)
    } else if (current_method == "DD") {
      if (empirical_start == TRUE) {
        trees <- try(DDtreeSim(numbsim=Numbsim1, lambda=Lambda, mu=Mu, K=K,  age=TreeAge, ddmodel=1), FALSE)
      } else if (empirical_start == FALSE) {
        print(paste("start sim trset ", trset, " lam ", lams, " mu ", mus, sep=""))
        trees <- DDtreeSim(numbsim=Numbsim1, lambda=Lambda, mu=Mu, K=N,  age=TreeAge, ddmodel=1)
        print(paste("done sim trset ", trset, " lam ", lams, " mu ", mus, sep=""))
      }
    } else if (current_method == "CD") {
      trees <- try(Newfunction(n=N, numbsim=Numbsim1, lambda=Lambda, mu=Mu, frac=1, complete=FALSE, stochsampling=TRUE), FALSE)
    } else if (current_method == "TraitD") {
      trees <- list()
      for (i in 1:Numbsim1) {
        repeat {
          tree <- list(try(tree.bisse(BiSSEpars, max.taxa=N, max.t=Inf, include.extinct=FALSE, x0=0), FALSE))
          if (sum(tree[[1]]$tip.state) != 0 & sum(tree[[1]]$tip.state) != length(tree[[1]]$tip.label)) {
            trees[[i]] <- tree[[1]]
            class(trees) <- "multiPhylo"
            break
          }
        }
      }
    }
  }
  if (class(trees) != "try-error" && !is.na(trees)) {
    class(trees) <- "multiPhylo"
  }
  prefix <- "Tree"
  suffix <- seq(1:length(trees))
  names(trees) <- paste(prefix, suffix, sep="_")
  trees
}


#' Get several DD trees out of DDD
#'
#' Internal function used by `GetMetricTrees` to simulate several trees under DD for each supplied tree or tree set.
#'
#' The function uses `dd_sim` from the package `DDD`.
#'
#' @param numbsim Number of simulations per tree
#' @param lambda Speciation rate
#' @param mu Extinction rate
#' @param K Carrying capacity
#' @param age Crown age for simulated tree
#' @param ddmodel Model for dependence of diversification rates on K, see `dd_sim` description.
#' @return A set of simulated trees
DDtreeSim <- function(numbsim, lambda, mu, K, age, ddmodel) {
  outtrees <- list()
  for (rounds in 1:numbsim) {
    onetree <- c()
    treenumber <- as.character(rounds)
    onetree <- dd_sim(c(lambda, mu, K), age, ddmodel)
    print(paste("done tree ", rounds, " of ", numbsim, sep=""))
    outtrees[[treenumber]] <- onetree$tes
  }
  outtrees
}
