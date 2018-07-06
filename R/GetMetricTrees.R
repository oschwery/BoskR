################################################
# Simulate trees based on empirical estimations
################################################
GetMetricTrees <- function(trset=trset, empirical_start=FALSE, current_case=NULL, empParams=empParams, current_method, current_method_est=current_method, N=NULL, Numbsim1, Lambda, Mu, l=NULL, a=NULL, LambdaFun=NULL, MuFun=NULL, TreeAge=NULL, BiSSEpars=NULL, tree=NULL) {
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

# display used metrics (only leave in for now as safety net... drop into a text file later)

#print(paste("Tree Used", current_case[trset], sep=": "))
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
#print(paste("LambdaFun", LambdaFun, sep=": "))
print(paste("MuFun", MuFun, sep=": "))
print(paste("TreeAge", TreeAge, sep=": "))
print(paste("BiSSEpars", BiSSEpars, sep=": "))
##############
#print("works until print all new params")
#try(print(class(trees)))
#try(print(trees))
##############
# pick right approach\
  if (is.na(N)) {
    trees <- NA
    ##############
#    print("works until is.na(N)")
#    print(class(trees))
#    print(trees)
    ##############

  } else if (!is.na(N)) {
    if (current_method == "BD") {
      trees <- try(sim.bd.taxa(n=N, numbsim=Numbsim1, lambda=Lambda, mu=Mu, frac=1, complete=FALSE, stochsampling=TRUE), FALSE)
    } else if (current_method == "TimeD-BD") {
      trees <- try(tess.sim.taxa.age(n=Numbsim1, nTaxa=N, age=TreeAge, lambda=LambdaFun, mu=MuFun), FALSE)
      #########
#      print("does it do anything in here?")
      #########
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
        #while (sum(trees[[i]]$tip.state) == 0 | sum(trees[[i]]$tip.state) == length(trees[[i]]$tip.label)) {
        repeat {
          tree <- list(try(tree.bisse(BiSSEpars, max.taxa=N, max.t=Inf, include.extinct=FALSE, x0=0), FALSE))
          if (sum(tree[[1]]$tip.state) != 0 & sum(tree[[1]]$tip.state) != length(tree[[1]]$tip.label)) {
            trees[[i]] <- tree[[1]]
            class(trees) <- "multiPhylo"
            break
          }
        }
        #}
      }
      #trees <- trees(BiSSEpars, type=c("bisse"), n=Numbsim1, max.taxa=N, max.t=Inf, include.extinct=FALSE, x0=0)
    }
  }
  ##############
#  print("worksuntil after loops of different methods")
  ##############
#print(class(trees))
#print(is.na(trees))
  if (class(trees) != "try-error" && !is.na(trees)) {
    class(trees) <- "multiPhylo"
    ##############
#    print("works turning trees to multiphylo")
#    print(class(trees))
    ##############
  }
  ##############
#  print("works until after multiPhylo")
  ##############

  prefix <- "Tree"
  suffix <- seq(1:length(trees))
  names(trees) <- paste(prefix, suffix, sep="_")
  ##############
#  print("works assigning trees")
  ##############

  trees
}
