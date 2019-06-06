#' Get diversification parameters from trees
#'
#' `GetTreeParams` estimates parameters from a supplied tree or tree set, which can subsequently be used as input for tree simulations using `GetMetricTrees`.
#'
#' The function wraps around the internal `GetParams`, and uses either
#'
#' @param trees Tree or set of trees, list or multiPhylo-object, or list of tree sets
#' @param current_method_est String specifying the method to be used to estimate the parameters. Can be `"BD", "TimeD-BD", "DD", "CD", "TraitD"` for birth-death, time-dependent birth-death, diversity dependent, clade dependent, or trait dependent diversification respectively.
#' @return A nested list of parameter estimates for every tree in `trees`, or every tree in each tree set therein respectively.
#'
#' @export

GetTreeParams <- function(trees, current_method_est) {
  outparams <- c()
  if (!is.list(trees[[1]][[1]][[1]])) {
    outparams <- GetParams(trees, current_method_est)
  } else if (is.list(trees[[1]][[1]][[1]])) {
    for (k in 1:length(trees$metricTreeSet)) {
      paramname <- c()
      paramsies <- c()
      print(names(trees$metricTreeSet[k]))
      if (is.na(trees$metricTreeSet[[k]])) {
        paramsies <- NA
      } else {
        paramsies <- try(GetParams(trees$metricTreeSet[[k]], current_method_est), FALSE)
      }
      paramname <- paste("params", names(trees$metricTreeSet[k]), sep="_")
      outparams[[paramname]] <- paramsies
    }
  }
  outparams
}


#' Internal to get diversification parameters from trees
#'
#' `GetParams` estimates parameters from a supplied tree or tree set, within the wrapper of 'GetTreeParams', which allows to use it for lists of tree sets (e.g. simulated trees based on a set of trees). The output can subsequently be used as input for tree simulations using `GetMetricTrees`.
#'
#' The function uses either
#'
#' @param emptrees Set of (probably empirical) phylogenies, list or multiPhylo-object.
#' @param current_method_est String specifying the method to be used to estimate the parameters. Can be `"BD", "TimeD-BD", "DD", "CD", "TraitD"` for birth-death, time-dependent birth-death, diversity dependent, clade dependent, or trait dependent diversification respectively.
#' @return A nested list of parameter estimates for every tree in `emptrees`.
#'
#' @noRd
#'
#' @import ape

GetParams <- function(emptrees, current_method_est) {
  params <- c()
  for (j in 1:length(emptrees)) {
    tree <- emptrees[j]
    current_name <- names(emptrees[j])
    # rate estimations empirical tree
    if (current_method_est == "BD") {
      empirical_solution <- try(BDredux(tree), FALSE)
    } else if (strsplit(x=current_method_est, split="_")[[1]][1] == "Time") {
      empirical_solution <- try(TimeDepBD(tree, current_method_est), FALSE)
    } else if (current_method_est == "DD") {
      empirical_solution <- try(DDredux(tree), FALSE)
    } else if (current_method_est == "CD") {
      empirical_solution <- try(Newfunction(tree), FALSE)  # WRITE THOSE
    } else if (current_method_est == "TraitD") {
      empirical_solution <- try(RateEstTraitDep(tree), FALSE)
    }
    #setting params
    if (class(empirical_solution) == "try-error") {
      N <- NA
      Lambda <- NA
      Mu <- NA
      K <- NA
      l <- NA
      a <- NA
      LambdaFun <- NA
      MuFun <- NA
      TreeAge <- NA
      BiSSEpars <- c(NA, NA, NA, NA, NA, NA)
      lik <- NA
      AIC <- NA
      params[[j]] <- list(N=N, Lambda=Lambda, Mu=Mu, K=K, l=l, a=a, LambdaFun=LambdaFun, MuFun=MuFun, TreeAge=TreeAge, BiSSEpars=BiSSEpars, lik=lik, AIC=AIC)
    } else {
      N <- length(tree[[1]]$tip.label)
      Lambda <- median(as.numeric(empirical_solution$lambda0), na.rm=TRUE)
      Mu <- median(as.numeric(empirical_solution$mu0), na.rm=TRUE)
      K <- as.numeric(empirical_solution[1,]$lambda1)
      l <- as.numeric(empirical_solution[1,]$lambda0)
      a <- as.numeric(empirical_solution[1,]$"a, k, etc")
      LambdaFun <- function(t) l * exp(a*t)
      MuFun <- as.numeric(empirical_solution[1,]$mu0)
      TreeAge <- max(branching.times(tree[[1]]))
      BiSSEpars <- c(median(as.numeric(empirical_solution$lambda0), na.rm=TRUE), median(as.numeric(empirical_solution$lambda1), na.rm=TRUE), median(as.numeric(empirical_solution$mu0), na.rm=TRUE), median(as.numeric(empirical_solution$mu1), na.rm=TRUE), median(as.numeric(empirical_solution$q01), na.rm=TRUE), median(as.numeric(empirical_solution$q10), na.rm=TRUE))
      lik <- as.numeric(empirical_solution[1,]$lik)
      AIC <- as.numeric(empirical_solution[1,]$AIC)
      params[[j]] <- list(N=N, Lambda=Lambda, Mu=Mu, K=K, l=l, a=a, LambdaFun=LambdaFun, MuFun=MuFun, TreeAge=TreeAge, BiSSEpars=BiSSEpars, lik=lik, AIC=AIC)
    }
  }
  params
}


#' Get BD diversification parameters from trees
#'
#' Internal function used by `GetParams` to estimate BD parameters from a supplied tree or tree set.
#'
#' The function uses `birthdeath` from the package `ape`.
#'
#' @param treeset Set of (probably empirical) phylogenies, list or multiPhylo-object.
#' @return A dataframe with BD parameters.
#'
#' @noRd
#'
#' @import ape

BDredux <- function(treeset) {
  outmatrix <- matrix(data=NA, nrow=length(treeset), ncol=11, dimnames=list(c(), c("Model", "Tree", "Method", "lambda0", "mu0", "lambda1", "mu1", "d/b (epsilon)", "b-d (r)", "lnLik", "AIC")))
  # ape
  ape_result <- list()
  for (i in 1:length(treeset)) {
    ape_result[[i]] <- list(birthdeath(treeset[[i]]))
  }
  # fill result into matrix
  for (i in 1:length(treeset)) {
    outmatrix[i,] <- c("bd", paste(deparse(substitute(treeset)),i, sep=" "), "ape_birthdeath", (ape_result[[i]][[1]]$para[2])/(1-(ape_result[[i]][[1]]$para[1])), ((ape_result[[i]][[1]]$para[2])*(ape_result[[i]][[1]]$para[1]))/(1-(ape_result[[i]][[1]]$para[1])), NA, NA, (ape_result[[i]][[1]]$para[1]), (ape_result[[i]][[1]]$para[2]), ((ape_result[[i]][[1]]$dev)/(-2)), ((-2)*((ape_result[[i]][[1]]$dev)/(-2))+(2*length(ape_result[[i]][[1]]$para))))
  }
  outframe <- as.data.frame(outmatrix, stringsAsFactors=FALSE)
  return(outframe=outframe)
}


#' Get TimeD-BD diversification parameters from trees
#'
#' Internal function used by `GetParams` to estimate TimeD-BD parameters from a supplied tree or tree set.
#'
#' The function uses `fit_bd` from the package `RPANDA`.
#'
#' @param treeset Set of (probably empirical) phylogenies, list or multiPhylo-object.
#' @return A dataframe with TimeD-BD parameters.
#'
#' @noRd
#'
#' @import ape
#' @importFrom picante node.age
#' @importFrom RPANDA fit_bd

TimeDepBD <- function(treeset, current_method_est) {
  outmatrix <- matrix(data=NA, nrow=length(treeset), ncol=11, dimnames=list(c(), c("Model", "Tree", "Method", "lambda0", "mu0", "lambda1", "mu1", "a, k, etc", "2nd a, k, etc", "lnLik", "AIC")))
  # RPANDA
  RPANDA_result <- list()
  for (i in 1:length(treeset)) {
    tot_time <- max(node.age(treeset[[i]])$ages)  # get max time for crown age
    # set Lambda formula and settings
    f.lamb <- c()
    lamb_par<-c()
    cst.lamb <- c()
    expo.lamb <- c()
    if (strsplit(x=current_method_est, split="_")[[1]][2] == "const") {
      f.lamb <- function(t,y){y[1]}
      lamb_par<-c(0.09)
      cst.lamb <- TRUE
      expo.lamb <- FALSE
    } else if (strsplit(x=current_method_est, split="_")[[1]][2] == "lin") {
      f.lamb <- function(t,y){y[1] + y[2] * t}
      lamb_par<-c(0.09, 0.001)
      cst.lamb <- FALSE
      expo.lamb <- FALSE
    } else if (strsplit(x=current_method_est, split="_")[[1]][2] == "exp") {
      f.lamb <- function(t,y){y[1] * exp(y[2] * t)}
      lamb_par <- c(0.05, 0.01)
      cst.lamb <- FALSE
      expo.lamb <- TRUE
    }
    # set mu formula and settings
    f.mu <- c()
    mu_par<-c()
    fix.mu <- c()
    cst.mu <- c()
    expo.mu <- c()
    if (strsplit(x=current_method_est, split="_")[[1]][3] == "PB") {
      f.mu <- function(t,y){0}
      mu_par<-c()
      fix.mu <- TRUE
      cst.mu <- FALSE
      expo.mu <- FALSE
    } else if (strsplit(x=current_method_est, split="_")[[1]][3] == "const") {
      f.mu <- function(t,y){y[1]}
      mu_par <- c(0.005)
      fix.mu <- FALSE
      cst.mu <- TRUE
      expo.mu <- FALSE
    } else if (strsplit(x=current_method_est, split="_")[[1]][3] == "lin") {
      f.mu <- function(t,y){y[1] + y[2] * t}
      mu_par <- c(0.005, 0.0001)
      fix.mu <- FALSE
      cst.mu <- FALSE
      expo.mu <- FALSE
    } else if (strsplit(x=current_method_est, split="_")[[1]][3] == "exp") {
      f.mu <- function(t,y){y[1] * exp(y[2] * t)}
      mu_par <- c(0.005, 0.0001)
      fix.mu <- FALSE
      cst.mu <- FALSE
      expo.mu <- TRUE
    }
    # run the model
    RPANDA_result[[i]] <- list(fit_bd(treeset[[i]], tot_time, f.lamb, f.mu, lamb_par, mu_par, f=1, meth = "Nelder-Mead", cst.lamb, cst.mu, expo.lamb, expo.mu, fix.mu, dt=1e-3, cond="crown"))
    RPANDA_result[[i]]$model <- current_method_est
    if (is.null(RPANDA_result[[i]]$mu_par)) {
      RPANDA_result[[i]]$mu_par <- NA
    }
  }
  # fill result into matrix
  for (i in 1:length(treeset)) {
    outmatrix[(i+(0*length(treeset))),] <- c(RPANDA_result[[i]]$model, paste(deparse(substitute(treeset)),i, sep=" "), paste("RPANDA fit_bd", "lambda", strsplit(x=current_method_est, split="_")[[1]][2], "mu", strsplit(x=current_method_est, split="_")[[1]][3] sep=" "), RPANDA_result[[i]]$lambda_par[1],  RPANDA_result[[i]]$mu_par[1], NA, NA, RPANDA_result[[i]]$lambda_par[2], RPANDA_result[[i]]$mu_par[2], RPANDA_result[[i]]$LH, RPANDA_result[[i]]$aicc)
  }
  outframe <- as.data.frame(outmatrix, stringsAsFactors=FALSE)
  return(outframe=outframe)
}

#' Get DD diversification parameters from trees
#'
#' Internal function used by `GetParams` to estimate DD parameters from a supplied tree or tree set.
#'
#' The function uses `dd_ML` from the package `DDD`.
#'
#' @param treeset Set of (probably empirical) phylogenies, list or multiPhylo-object.
#' @return A dataframe with TimeD-BD parameters.
#'
#' @noRd
#'
#' @import ape
#' @importFrom laser getBtimes
#' @importFrom DDD dd_ML

DDredux <- function(treeset) {
  outmatrix <- matrix(data=NA, nrow=length(treeset), ncol=11, dimnames=list(c(), c("Model", "Tree", "Method", "lambda0", "mu0", "lambda1", "mu1", "d/b (epsilon)", "b-d (r)", "lnLik", "AIC")))
    # DDD1
    DDD1_result <- list()
    for (i in 1:length(treeset)) {
      branching_times <- getBtimes(string=write.tree(treeset[[i]]))
      DDD1_result[[i]] <- list(dd_ML(branching_times, ddmodel = 1, cond = 1, soc = 2))
    }
    # fill result into matrix
    for (i in 1:length(treeset)) {
      outmatrix[i,] <- c("dd", paste(deparse(substitute(treeset)),i, sep=" "), "DDD1_dd_ML ddmodel=1 cond=1", DDD1_result[[i]][[1]]$lambda, DDD1_result[[i]][[1]]$mu, DDD1_result[[i]][[1]]$K, NA, (DDD1_result[[i]][[1]]$mu/DDD1_result[[i]][[1]]$lambda), (DDD1_result[[i]][[1]]$lambda-DDD1_result[[i]][[1]]$mu), DDD1_result[[i]][[1]]$loglik, (-2*(DDD1_result[[i]][[1]]$loglik)+(2*length(DDD1_result[[i]][[1]]$df))))
    }
  outframe <- as.data.frame(outmatrix, stringsAsFactors=FALSE)
  return(outframe=outframe)
}


#' Get TraitD diversification parameters from trees
#'
#' Internal function used by `GetParams` to estimate TraitD parameters from a supplied tree or tree set.
#'
#' The function uses `make.bisse` and `find.mle` from the package `diversitree`.
#'
#' @param treeset Set of (probably empirical) phylogenies, list or multiPhylo-object.
#' @return A dataframe with TimeD-BD parameters.
#'
#' @noRd
#'
#' @import ape
#' @import diversitree

RateEstTraitDep <- function(treeset) {
  outmatrix <- matrix(data=NA, nrow=length(treeset), ncol=11, dimnames=list(c(), c("Model", "Tree", "Method", "lambda0", "mu0", "lambda1", "mu1", "q01", "q10", "lnLik", "AIC")))
  # diversitree2  THERE ARE WARNINGS FOR THE CODE, possibly related to solution close to start values
  divtree_result <- list()
  for (i in 1:length(treeset)) {
    lik <- make.bisse(treeset[[i]], treeset[[i]]$tip.state)
    lik.m<- constrain(lik,  mu0 ~ mu1)
    p <- starting.point.bisse(treeset[[i]])
    divtree_result[[i]] <- list(find.mle(lik.m, p[argnames(lik.m)]))
    lik <- c()
    lik.m <- c()
    print(paste(i, "has passed", sep=" "))
  }
  # fill result into matrix
  for (i in 1:length(treeset)) {
    outmatrix[i,] <- c("trait-dep bd", paste(deparse(substitute(treeset)),i, sep=" "), "diversitree_make.bisse", divtree_result[[i]][[1]]$par.full[1], divtree_result[[i]][[1]]$par.full[3], divtree_result[[i]][[1]]$par.full[2], divtree_result[[i]][[1]]$par.full[4], divtree_result[[i]][[1]]$par.full[5], divtree_result[[i]][[1]]$par.full[6], divtree_result[[i]][[1]]$lnLik, (-2*(divtree_result[[i]][[1]]$lnLik)+(2*length(divtree_result[[i]][[1]]$par))))
  }
  outframe <- as.data.frame(outmatrix, stringsAsFactors=FALSE)
  result.list <- list(divtree_result=divtree_result)
  save(result.list, file=paste("/Users/", whatsMyName, "/Documents/", pathToVictory, "/PhD-Projects/DivModelTesting/Outputs_", StartDate,  "/", current_name, "_", current_method, "/","result_list_", current_name, "_", deparse(substitute(treeset)),".Rdata", sep=""))
    save(outframe, file=paste("/Users/", whatsMyName, "/Documents/", pathToVictory, "/PhD-Projects/DivModelTesting/Outputs_", StartDate,  "/", current_name, "_", current_method, "/","result_outframe_", current_name, "_", deparse(substitute(treeset)),".Rdata", sep=""))
  return(outframe=outframe)
}
