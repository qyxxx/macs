#geno is a data.frame with its genes' name
tarv <- function(pheno, geno, formula, method = c("entropy", "gini"),
                  family = c("binomial", "gaussian"),
                  direction = c("both", "positive", "negative"),
                  alpha = 0.01, cost = NULL){
  cl <- match.call()
  if(missing(pheno)) stop("'pheno' argument is required")
  if(missing(geno)) stop("'geno' argument is required")
  if(missing(formula)) stop("'formula' argument is required")
  if(match(0, cost, nomatch = 0)) stop("element of 'cost' must be nonzero")
  if(missing(method)) method <- "entropy"
  if(missing(family)) family <- "binomial"
  if(missing(direction)) direction <- "both"

  if(method=="entropy"){
    srule <- 1
  }else if(method=="gini"){
    srule <- 2
  }else{
    stop("wrong 'method'")
  }
  if(!match(family, c("binomial", "gaussian"), nomatch = 0)) {
    stop("wrong 'family'")
  }
  if(!match(direction, c("both", "positive", "negative"), nomatch = 0)) {
    stop("wrong 'direction'")
  }
  if(nrow(pheno) != ncol(geno)) stop("nrow(pheno) != ncol(geno)")
  direction <- match.arg(direction)
  Z <- TARV_marginal_effect(pheno, geno, formula, family)
  annotation <- rownames(geno)
  X <- TARV_transform(geno, annotation, Z, direction)
  data <- cbind(pheno,X)
  colname <- rep(1, ncol(data))
  colname[1] <- -1
  index <- which(sapply(data, class)=="factor")
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
  colnames(out$final_counts) <- paste0("count_", 1:ncol(out$final_counts))
  out$final_counts <- data.frame(out$final_counts)
  names(out$varcatg) <- 1:length(out$varcatg)
  names(out$nodeclass) <- 1:out$nnd
  names(out$p_value) <- 1:out$nnd
  out <- c(out, family = list(family), direction = list(direction), method = list(method),
           formula = list(formula), call = list(cl), learning.data = list(data))
  class(out) <- "ctree"
  out
}
