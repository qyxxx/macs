upfun <- function(x, tau){
  return(ifelse((x-tau)>0, x-tau, 0))
}

process <- function(out){
  finTerm <- out$finTerm
  beta <- out$beta
  nterm <- length(beta)
  varList <- out$varList
  fitted <- out$fitted
  data <- out$data
  x <- data[, 2:(ncol(data)-1)]

  if(nterm>2){
    f <- substring(fitted, 1:nchar(fitted), 1:nchar(fitted))
    col <- NULL
    i <- 1
    while(i<=length(f)){
      if(f[i]=="x"){
        j <- i
        while(f[j]!=" "&j<=length(f)) j <- j + 1
        col <- c(col, paste(f[i:(j-1)], collapse = ""))
        i <- j
      }
      else i <- i + 1
    }
    z <- matrix(1, nrow = nrow(x), ncol = length(beta) - 2)
    n <- 1
    for(i in 3:nterm){
      j <- finTerm[i]
      k <- varList[j+1, 1]
      for (l in 1:k) {
        if(varList[j+1, 3*l] > 0.){
          z[, n] <- z[, n]*upfun(x[, varList[j+1, 3*l-1]], varList[j+1, 3*l+1])
        }
        else
          z[, n] <- z[, n]*x[, varList[j+1, 3*l-1]]
      }
      n <- n + 1
    }
    colnames(z) <- col
    z <- data.frame(y = data[, ncol(data)], z)
  }
  else z <- data.frame(y = data[, ncol(data)]) # intercept only
  z
}
