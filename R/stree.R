stree<-function(data, family = c("likelihood", "log-rank",
                "kaplan-meiser distance", "adaptive normalization",
                "global normalization")){
  cl <- match.call()
  if(missing(data)) stop("'data' argument is required")
  if(missing(family)) family <- "likelihood"

  srule <- match(family, c("likelihood", "log-rank",
            "kaplan-meiser distance", "adaptive normalization",
            "global normalization"), nomatch = 0)
  if(!srule) stop("wrong 'family'")
  colname <- rep(1, ncol(data))

  index <- which(sapply(data, class)=="factor")
  for(i in index){
    colname[i] <- 2
  }
  colname[1] <- -2 # response
  colname[2] <- -1 # censoring indicator (1=death, 0=alive)
  output <- cstree(data, colname, srule)
  out <- trim_stree(output)
  names(out$dt) <- 1:out$nnd
  names(out$pt) <- 1:out$nnd
  names(out$spv) <- 1:out$nnd
  colnames(out$spvl) <- 1:out$nnd
  names(out$ncases) <- 1:out$nnd
  names(out$death_catg) <- 1:out$nnd
  names(out$median) <- 1:out$nnd
  names(out$colname) <- 1:length(out$colname)
  out <- c(out, family = list(family), call = list(cl), learning.data = list(data))
  class(out) <- "stree"
  out
}
