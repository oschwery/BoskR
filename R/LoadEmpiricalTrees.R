# Function to calculate tree characteristics in "treespace"
LoadEmpiricalTrees <- function(current_case) {
  returntree <- c()
  if (length(current_case) > 1) {
    for (i in 1:length(current_case)) {
      notreturntree <- LoadEmpiricalTreesSingle(current_case[i])
      returntree[[i]] <- notreturntree
    }
  } else if (current_case == "combo") {
    load(file=paste("/Users/", whatsMyName, "/Documents/", pathToVictory, "/PhD-Projects/DivModelTesting/TreeSpaces/EmptreesCombined.Rdata", sep=""))
    returntree <- emptrees
  } else if (current_case == "ALL") {
    for (i in 1:11) {  # 11 needs to be updated if more trees added...
      notreturntree <- LoadEmpiricalTreesSingle(i)
      returntree[[i]] <- notreturntree
    }
  } else if (current_case == "otl") {
    load(file=paste("/Users/", whatsMyName, "/Documents/", pathToVictory, "/PhD-Projects/DivModelTesting/TreeSpaces/goodsizetrees.Rdata", sep=""))
    returntree <- goodsizetrees
  } else if (length(current_case) == 1) {
    notreturntree <- LoadEmpiricalTreesSingle(current_case)
    returntree[[1]] <- notreturntree
  }
  return(returntree)
}

LoadEmpiricalTreesSingle <- function(current_case) {
  if (current_case == 1) {
    library(BAMMtools)
    data(whales)
    tree <- list(whales)
    current_name <- "whales"
  } else if (current_case == 2) {
    library(MonoPhy)
    data(Ericactree)
    tree <- list(Ericactree)
    current_name <- "Ericac77"
  } else if (current_case == 3) {
    Ericac450 <- read.tree(paste("/Users/", whatsMyName, "/Documents/", pathToVictory, "/PhD-Projects/DivModelTesting/EmpiricalTrees/Ericaceae_MCC_with_outgroup.tre", sep=""))
    tree <- list(Ericac450)
    current_name <- "Ericac450"
  } else if (current_case == 4) {
    Poales <- read.tree(paste("/Users/", whatsMyName, "/Documents/", pathToVictory, "/PhD-Projects/DivModelTesting/EmpiricalTrees/MCC_Poales.tre", sep=""))
    tree <- list(Poales)
    current_name <- "Poales"
  } else if (current_case == 5) {
    Fagales <- read.tree(paste("/Users/", whatsMyName, "/Documents/", pathToVictory, "/PhD-Projects/DivModelTesting/EmpiricalTrees/MCC_Fagales.tre", sep=""))
    tree <- list(Fagales)
    current_name <- "Fagales"
  } else if (current_case == 6) {
    Angios <- read.nexus(paste("/Users/", whatsMyName, "/Documents/", pathToVictory, "/PhD-Projects/DivModelTesting/EmpiricalTrees/floral_1.nex", sep=""))
    tree <- list(Angios)
    current_name <- "Angios"
  } else if (current_case == 7) {
    Liverworths <- read.nexus(paste("/Users/", whatsMyName, "/Documents/", pathToVictory, "/PhD-Projects/DivModelTesting/EmpiricalTrees/T75503.nex", sep=""))
    tree <- list(Liverworths)
    current_name <- "Liverworths"
  } else if (current_case == 8) {
    library(RPANDA)
    data(Balaenopteridae)
    tree <- list(Balaenopteridae)
    current_name <- "Balaenopteridae"
  } else if (current_case == 9) {
    library(RPANDA)
    data(Calomys)
    tree <- list(Calomys)
    current_name <- "Calomys"
  }else if (current_case == 10) {
    library(RPANDA)
    data(Phocoenidae)
    tree <- list(Phocoenidae)
    current_name <- "Phocoenidae"
  }else if (current_case == 11) {
    library(RPANDA)
    data(Phyllostomidae)
    tree <- list(Phyllostomidae)
    current_name <- "Phyllostomidae"
  }
return(list(tree=tree, current_name=current_name))
}
# Tree selection catalogue
# current_case --> tree  --> add provenance later?
#  1 whales
#  2 Ericac77
#  3 Ericac450
#  4 Poales
#  5 Fagales
#  6 Angios
#  7 Liverworths
#  8 Balaenopteridae
#  9 Calomys
# 10 Phocoenidae
# 11 Phyllostomidae
CorrUltramet <- function(emptrees) {
  for (i in 1:length(emptrees)) {
#    tree <- emptrees[[i]]$tree[[1]]
    tree <- emptrees[[i]]
    nnls <- c()
    if (is.ultrametric(tree) == FALSE) {
#      print(paste("not ultrametric", emptrees[[i]]$current_name, sep=" "))
      print(paste("not ultrametric", names(emptrees[i]), sep=" "))
      nnls<-nnls.tree(cophenetic(tree),tree,rooted=TRUE)
      if (is.ultrametric(nnls) == TRUE) {
#        emptrees[[i]]$tree[[1]] <- nnls
        emptrees[[i]] <- nnls
#        print(paste("fixed", emptrees[[i]]$current_name, sep=" "))
        print(paste("fixed", names(emptrees[i]), sep=" "))
      } else if (is.ultrametric(nnls) == FALSE) {
#        print(paste("Still not ultrametric", emptrees[[i]]$current_name, sep=" "))
        print(paste("Still not ultrametric", names(emptrees[i]), sep=" "))
      }
    }
    tree <- c()
  }
  emptrees
}
# #make tree ultrametric again (rounding errors)
# nnls<-nnls.tree(cophenetic(tree),tree,rooted=TRUE)
# ## check
# is.ultrametric(nnls)
# tree <- nnls
CorrZerobranch <- function(emptrees) {
  for (i in 1:length(emptrees)) {
#    tree <- emptrees[[i]]$tree[[1]]
    tree <- emptrees[[i]]
    nnls <- c()
    if (is.binary(tree) == FALSE) {
      print(paste("not binary", names(emptrees[i]), sep=" "))
      newtree <- multi2di(tree, random=TRUE)
      if (is.binary(newtree) == TRUE) {
        emptrees[[i]] <- newtree
        print(paste("fixed", names(emptrees[i]), sep=" "))
      } else if (is.binary(newtree) == FALSE) {
        print(paste("Still not binary", names(emptrees[i]), sep=" "))
      }
    }
    tree <- c()
  }
  emptrees
}
# reorder all to cladewise, bcz this suddenly causes a problem
ReorderCladewise <- function(emptrees) {
  for (i in 1:length(emptrees)) {
    emptrees[[i]] <- reorder.phylo(emptrees[[i]],"cladewise")
  }
  emptrees
}
