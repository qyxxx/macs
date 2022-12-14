#' @export
forecast <- function(object, ...) UseMethod("forecast")

forecast.ctree <- function(object, newdata, ...){
  if(!inherits(object, "ctree")) stop("not a legitimate \" ctree \" object")
  if(missing(newdata) || is.null(newdata)){
    newdata <- object$learning.data[, -1]
  }
  varcatg <- object$varcatg
  if(ncol(newdata)!=(length(varcatg)-1))
    stop("the column number of input data does not match")
  index <- ifelse(sum(varcatg>1), which(varcatg>1)-1, 0)
  index_new <- ifelse(sum(sapply(newdata, class)=="factor"), which(sapply(newdata, class)=="factor"), 0)
  if(sum(index_new!=index)) stop("wrong input data: check the 'factor' variable")
  if(index!=0)
  {
    for(i in index){
      newdata[, i] <- as.integer(newdata[, i])
    }
  }
  dt <- object$dt
  spv <- object$spv
  spvl <- object$spvl
  nodeclass <- object$nodeclass

  spv0 <- rep(0, length(spv))
  for(i in 1:length(spv)){
    if(!is.na(spv[i])){
      len <- nchar(spv[i])
      spv0[i] <- as.numeric(substr(spv[i], 2, len))
    }
  }

  spvl0 <- list()
  for(i in 1:ncol(spvl)){
    if(is.na(spvl[i])) spvl0 <- c(spvl0, list(NA))
    else{
      if(is.character(spvl[, i]))
      {
        val <- substring(spvl[i], 1:nchar(spvl[i]), 1:nchar(spvl[i]))
        n <- NULL
        j <- 2
        while(j<length(val)){
          for(k in j:length(val)){
            if(val[k]==","|val[k]=="}"){
              n <- c(n, as.numeric(paste(val[j:(k-1)], collapse = "")))
              break
            }
          }
          j <- k+1
        }
        spvl0 <- c(spvl0, list(n))
      }
      else spvl0 <- c(spvl0, list(spvl[, i]))
    }
  }

  node <- NULL
  predictor <- NULL
  for(i in 1:nrow(newdata)){
    temp <- 1
    while(dt[temp]!=0){
      if(varcatg[spv0[temp]+1]>1){
        if(newdata[i, spv0[temp]] %in% spvl0[[temp]])
          temp <- dt[temp]+2
        else
          temp <- dt[temp]+1
      }else{
        if(newdata[i, spv0[temp]]>spvl0[[temp]])
          temp <- dt[temp]+2
        else
          temp <- dt[temp]+1
      }
    }
    node <- c(node, temp)
    predictor <- c(predictor, nodeclass[temp])
  }
  z <- data.frame(node = node, predictor = predictor)
  rownames(z) <- rownames(newdata)
  class(z) <- "forecast.ctree"
  z
}

forecast.stree <- function(object, newdata, ...){
  if(!inherits(object, "stree")) stop("not a legitimate \" stree \" object")
  if(missing(newdata) || is.null(newdata)){
    newdata <- object$learning.data[, c(-1, -2)]
  }
  if(ncol(newdata)!=(length(object$colname)-2))
    stop("the column number of input data does not match")
  index <- ifelse(sum(object$colname>1), which(object$colname>1)-2, 0)
  index_new <- ifelse(sum(sapply(newdata, class)=="factor"), which(sapply(newdata, class)=="factor"), 0)
  if(sum(index_new!=index)) stop("wrong input data: check the 'factor' variable")
  for(i in index) newdata[, i] <- as.integer(newdata[, i])
  dt <- object$dt
  spv <- object$spv
  spvl <- object$spvl
  median <- object$median
  colname <- object$colname
  spv0 <- NULL
  for(i in 1:length(spv)){
    if(is.na(spv[i])) spv0 <- c(spv0, 0)
    else{
      spv0 <- c(spv0, as.numeric(substring(spv[i], 2, nchar(spv[i]))))
    }
  }

  spvl0 <- list()
  for(i in 1:ncol(spvl)){
    if(is.na(spvl[i])) spvl0 <- c(spvl0, list(NA))
    else{
      if(is.character(spvl[, i]))
      {
        val <- substring(spvl[i], 1:nchar(spvl[i]), 1:nchar(spvl[i]))
        n <- NULL
        j <- 2
        while(j<length(val)){
          for(k in j:length(val)){
            if(val[k]==","|val[k]=="}"){
              n <- c(n, as.numeric(paste(val[j:(k-1)], collapse = "")))
              break
            }
          }
          j <- k+1
        }
        spvl0 <- c(spvl0, list(n))
      }
      else spvl0 <- c(spvl0, list(spvl[, i]))
    }
  }
  node <- NULL
  predictor <- NULL
  for(i in 1:nrow(newdata)){
    temp <- 1
    while(dt[temp]!=0){
      if(colname[spv0[temp]+2]>1){
        if(newdata[i, spv0[temp]] %in% spvl0[[temp]])
          temp <- dt[temp]+2
        else
          temp <- dt[temp]+1
      }
      else{
        if(newdata[i, spv0[temp]] > spvl0[[temp]])
          temp <- dt[temp]+2
        else
          temp <- dt[temp]+1
      }
    }
    node <- c(node, temp)
    predictor <- c(predictor, median[temp])
  }
  z <- data.frame(node = node, predictor = predictor)
  rownames(z) <- rownames(newdata)
  class(z) <- "forecast.stree"
  z
}

forecast.masal <- function(object, newdata = NULL, col = NULL, ...){
  if(!inherits(object, "masal")) stop("not a legitimate \" masal \" object")
  if(!is.null(col)){
    if(is.null(newdata)) stop("when 'col' given, 'newdata' is required")
    s_index <- which(col=="s")# sort
    if(length(s_index)!=1)
      stop("input one sort-category in 'col'")
    t_index <- which(col=="t")# time
    if(length(t_index)!=1)
      stop("input one time-category in 'col'")
    c_index <- which(col=="c")# covariate
    if(length(c_index)==0)
      stop("input covaraite-category in 'col'")
    covariate <- data.frame(newdata)
    covariate[, 1:(ncol(newdata)-2)] <- newdata[, c_index]
    covariate[, ncol(newdata)-1] <- newdata[, s_index]
    covariate[, ncol(newdata)] <- newdata[, t_index]
  }else{
    covariate <- newdata
  }
  if(is.null(covariate)) covariate <- object$learning.data[, c(-1, -ncol(object$learning.data))]
  if(ncol(covariate)!=(ncol(object$learning.data)-2))
    stop("the column number of input data does not match")
  fitted <- object$fitted
  coefficients <- object$coefficients
  cov_matrix <- object$cov_matrix
  f <- substring(fitted, 1:nchar(fitted), 1:nchar(fitted))

  f <- unlist(strsplit(fitted, split = "\\s+"))[-1] # terms (intercept dropped)
  f <- f[nchar(f)>1] # drop "+" & "-"
  if(!length(f)){ # intercept only
    data_const <- rep(1, nrow(covariate))
    masaled.values <- drop(data_const*coefficients)
  }
  else {
    sf <- strsplit(f, split = "") # list
    z <- NULL
    for(i in 1:length(f)){ # number of terms
      fi <- f[i]
      sfi <- sf[[i]]
      lbrace <- which(sfi=="{") + 1
      rbrace <- which(sfi=="}") - 1
      n <- as.numeric(substring(fi, lbrace, rbrace))
      if("(" %in% sfi){
        rbracket <- which(sfi==")") - 1
        tau <- substring(fi, rbrace + 2, rbracket)
        tauf <- substring(tau, 1, 1)
        tau[tauf=="("|tauf==""] <- NA
        tau <- as.numeric(tau)
      }
      else tau <- rep(NA, length(n))
      xi <- 1
      for(j in 1:length(n)) {
        if(!is.na(tau[j])) xi <- xi*ifelse((covariate[, n[j]] + tau[j])>0, covariate[, n[j]] + tau[j], 0)
        else xi <- xi*covariate[, n[j]]
      }
      z <- cbind(z, xi)
    }
    data_const <- cbind(rep(1, nrow(z)), z)
    masaled.values <- drop(data_const %*% coefficients)
  }
  names(masaled.values) <- rownames(covariate)
  class(masaled.values) <- "forecast.masal"
  masaled.values
}
