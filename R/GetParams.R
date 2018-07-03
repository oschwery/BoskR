# Function to get Parameters from trees
GetParams <- function(emptrees, current_method_est) {
  params <- c()
  # load main analysis script
  if (current_method_est == "BD") {
    source(paste("/Users/", whatsMyName, "/Documents/", pathToVictory, "/PhD-Projects/DivModelTesting/TreeSpaces/BDredux.R", sep=""))
#  source(paste("/Users/", whatsMyName, "/Documents/", pathToVictory, "/PhD-Projects/DivModelTesting/AllTanja41TreeSet.R", sep=""))
  } else if (current_method_est == "TimeD-BD") {
    source(paste("/Users/", whatsMyName, "/Documents/", pathToVictory, "/PhD-Projects/DivModelTesting/TreeSpaces/TimeD-BDredux.R", sep=""))
#    source(paste("/Users/", whatsMyName, "/Documents/", pathToVictory, "/PhD-Projects/DivModelTesting/RateEstTimeDep.R", sep=""))
  } else if (current_method_est == "DD") {
    source(paste("/Users/", whatsMyName, "/Documents/", pathToVictory, "/PhD-Projects/DivModelTesting/TreeSpaces/DDredux.R", sep=""))
#    source(paste("/Users/", whatsMyName, "/Documents/", pathToVictory, "/PhD-Projects/DivModelTesting/Newscript.R", sep=""))
  } else if (current_method_est == "CD") {
    source(paste("/Users/", whatsMyName, "/Documents/", pathToVictory, "/PhD-Projects/DivModelTesting/Newscript.R", sep=""))
  } else if (current_method_est == "TraitD") {
    source(paste("/Users/", whatsMyName, "/Documents/", pathToVictory, "/PhD-Projects/DivModelTesting/RateEstTraitDep.R", sep=""))
  }
  for (j in 1:length(emptrees)) {
    tree <- emptrees[j]
    current_name <- names(emptrees[j])
    # rate estimations empirical tree
    if (current_method_est == "BD") {
      empirical_solution <- try(BDredux(tree), FALSE)
#      empirical_solution <- AllTanja41TreeSet(tree, FALSE)
    } else if (current_method_est == "TimeD-BD") {
      empirical_solution <- try(TimeDBDredux(tree), FALSE)
#      empirical_solution <- RateEstTimeDep(tree)
    } else if (current_method_est == "DD") {
      empirical_solution <- try(DDredux(tree), FALSE)
#      empirical_solution <- Newfunction(tree)  # WRITE THOSE
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
      params[[j]] <- list(N=N, Lambda=Lambda, Mu=Mu, K=K, l=l, a=a, LambdaFun=LambdaFun, MuFun=MuFun, TreeAge=TreeAge, BiSSEpars=BiSSEpars)
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
      params[[j]] <- list(N=N, Lambda=Lambda, Mu=Mu, K=K, l=l, a=a, LambdaFun=LambdaFun, MuFun=MuFun, TreeAge=TreeAge, BiSSEpars=BiSSEpars)
    }
  }
  params
}
