# script to combine all trees (or add them)
# Function to calculate tree characteristics in "treespace"
CombineEmpiricalTrees <- function(set, location=NULL) {
  emptrees <- list()
  names <- c()
  Ntax <- c()
  Age <- c()
  if (set == "initial") {
    emptreesALL <- LoadEmpiricalTrees("ALL")  # load package trees
    for (i in 1:length(emptreesALL)) {
      emptrees[i] <- emptreesALL[[i]]$tree
      names <- c(names, emptreesALL[[i]]$current_name)
      Ntax <- c(Ntax, length(emptreesALL[[i]][[1]][[1]]$tip.label))
      Age <- c(Age, max(branching.times(emptreesALL[[i]][[1]][[1]])))
    }
    emptreesOTL <- LoadEmpiricalTrees("otl")
    for (i in 1:length(emptreesOTL)) {
      emptrees[length(emptreesALL)+i] <- emptreesOTL[[i]][[1]]
      names <- c(names, emptreesOTL[[i]][[2]])
      Ntax <- c(Ntax, length(emptreesOTL[[i]][[1]][[1]]$tip.label))
      Age <- c(Age, max(branching.times(emptreesOTL[[i]][[1]][[1]])))
    }
    prefix <- "Emp"
    suffix <- seq(1:length(emptrees))
    identifiers <- paste(prefix, suffix, sep="_")
    # print(length(identifiers))
    # print(length(Ntax))
    # print(length(Age))
    # print(length(names))
    names(emptrees) <- identifiers
    treeindex <- data.frame(TreeID=identifiers, Ntax=Ntax, Age=Age, Source=names)
    save(emptrees, file=paste("/Users/", whatsMyName, "/Documents/", pathToVictory, "/PhD-Projects/DivModelTesting/TreeSpaces/EmptreesCombined.Rdata", sep=""))
    save(treeindex, file=paste("/Users/", whatsMyName, "/Documents/", pathToVictory, "/PhD-Projects/DivModelTesting/TreeSpaces/Treeindex.Rdata", sep=""))
  } else if (set =="addon") {
    #to be continued
  }
  return(paste("Saved ", length(emptrees), " trees to file", sep=""))
}
