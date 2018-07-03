# function to get several DD trees out of DDD
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
