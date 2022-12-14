sort_covariate <- function(model){
  id_all <- paste("person", model$data[, 1])
  data_processed <- model$data_processed
  num <- data_processed[, ncol(data_processed)-2]
  ntime <- model$ntime
  cov_matrix <- model$cov_matrix
  z <- NULL
  id <- NULL
  j <- 1
  for(i in 1:length(ntime)){
    id <- c(id, id_all[j])
    n <- num[j:(j+ntime[i]-1)]
    r <- rbind(n, cov_matrix[[i]])
    if(dim(r)[2]!=1){
      r <- r[, order(r[1, ])]
      r <- r[2:dim(r)[1], ]
      r <- cbind(n, r)
      r <- r[order(r[, 1]), ]
      r <- r[, 2:dim(r)[2]]
      z <- c(z, list(r))
      j <- j + ntime[i]
    }
    else{
      r <- r[-1,]
      r <- cbind(n,r)
      r <- r[,-1]
      z = c(z, list(r))
      j <- j + ntime[i]
    }
  }
  names(z) <- id
  z
}
