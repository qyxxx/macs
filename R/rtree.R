rtree <- function(data, method = c("entropy", "gini"), alpha = 0.01,
                  cost = NULL)
{
  cl <- match.call()
  if(missing(data)) stop("'data' argument is required")
  if(missing(method)) method <- "entropy"

  if(method=="entropy"){
    srule <- 1
  }else if(method=="gini"){
    srule <- 2
  }else{
    stop("wrong 'method'")
  }
  index <- which(sapply(data, class)=="factor")
  colname <- rep(1, ncol(data))
  colname[1] <- -1
  for(i in index){
    colname[i] <- length(unique(data[, i]))
  }
  if(is.null(cost)) cost <- rep(1, length(unique(data[, 1])))
  else if(length(cost)!=length(unique(data[, 1]))){
    stop("the length of 'cost' does not match the input data")
  }
  output <- cctree(data, colname, srule, alpha)
  p_value <- pchisq(output$chi[1:output$nnd], length(unique(data[, 1])) - 1)
  out <- trim_ctree(output)
  classification <- classify(out, cost)
  out$nodeclass <- classification$nodeclass
  out$p_value <- p_value
  out$cost <- cost
  names(out$cost) <- 1:length(cost)
  names(out$dt) <- 1:out$nnd
  names(out$pt) <- 1:out$nnd
  names(out$spv) <- 1:out$nnd
  colnames(out$spvl) <- 1:out$nnd
  names(out$p_value) <- 1:out$nnd
  colnames(out$final_counts) <- paste0("count_", 1:ncol(out$final_counts))
  out$final_counts <- data.frame(out$final_counts)
  names(out$varcatg) <- 1:length(out$varcatg)
  names(out$nodeclass) <- 1:out$nnd
  names(out$p_value) <- 1:out$nnd
  out <- c(out, method = list(method), call = list(cl), learning.data = list(data))
  class(out) <- "ctree"
  out
}
