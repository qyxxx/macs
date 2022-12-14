classify <- function(fit, cost){
  nnd <- fit$nnd
  final_counts <- fit$final_counts
  dt <- fit$dt

  num <- sum(final_counts[1, ])
  nodeclass <- rep(0, nnd)
  p_value <- rep(0, nnd)
  leavenum <- rep(1, nnd)
  R_s <- matrix(0, nnd,2)
  for(i in 1:nnd){
    index <- which.max(cost*final_counts[i,])
    R_s[i,1] <- sum(cost*final_counts[i,])-cost[index]*final_counts[i,index]
    nodeclass[i] <- index
  }
  for(i in nnd:1){
    if(dt[i]>0){
      leavenum[i] = leavenum[dt[i]+1] + leavenum[dt[i]+2]
      R_s[i, 2] <- R_s[dt[i]+1, 2] + R_s[dt[i]+2, 2]
    }else{
      R_s[i, 2] <- R_s[i, 1]
    }
  }
  for(i in 1:nnd){
    if(dt[i]>0){
      p_value[i] <- (R_s[i, 1]-R_s[i, 2])/(leavenum[i]-1)
    }
  }
  out <- list(nodeclass = nodeclass, p_value = p_value)
  out
}
