enforce_test <- function(object, number){
  if(!is.null(number)){
    formula <- object$fitted
    f <- substring(formula, 1:nchar(formula), 1:nchar(formula))
    x <- character(0)
    number <- as.character(number)
    j <- 0
    for(i in 1:length(f)){
      if(f[i]=="{"&i<=length(f)){
        j <- j+1
        x[j] <- f[i+1]
      }
    }

    k <- rep(0,length(number))

    for(i in 1:length(number)) {
      for(j in 1:length(x)) {
        if(number[i]==x[j]) {
          k[i] <- 1 #forced coviriate was included in the final model
          break
        }
      }
    }

    if(length(which(k==0))==0) cat("all forced covariates were included in the final model.\n")
    else {
      k <- number[which(k==0)]
      k <- as.numeric(k) #k save for the covariates number

      data <- as.matrix(object$model)
      orig <- object$learning.data

      if(length(k)==1) noenforce <- paste0("x_{", k[1], "}")
      else{
        for(i in 2:length(k)) noenforce <- paste0("x_{", k[1], "}", "x_{", k[i], "}")
      }
      cat(paste("covariates which were forced but not added in the final model are", noenforce), ".\n")
    }
  }
}
