masal <- function(data, col = NULL, maxterm = 20, maxinteraction = 2,
                  varhat = c(1, 1, 1), forced_num = NULL, maxiter = 2){
  cl <- match.call()
  if(missing(data)) stop("'data' argument is required")
  num <- forced_num
  if(!is.null(col)){
    s_index <- which(col=="s")# sort
    if(length(s_index)!=1)
      stop("input one sort-category in 'col'")
    i_index <- which(col=="i")# ID
    if(length(i_index)!=1)
      stop("input one ID-category in 'col'")
    t_index <- which(col=="t")# time
    if(length(t_index)!=1)
      stop("input one time-category in 'col'")
    r_index <- which(col=="r")# response
    if(length(r_index)!=1)
      stop("input one response-category in 'col'")
    c_index <- which(col=="c")# covariate
    if(length(c_index)==0)
      stop("input covariate-category in 'col'")
    data0 <- data.frame(nrow(data), ncol(data))
    data0[, 1] <- data[, i_index]
    data0[, 2:(ncol(data)-3)] <- data[, c_index]
    data0[, ncol(data)-2] <- data[, s_index]
    data0[, ncol(data)-1] <- data[, t_index]
    data0[, ncol(data)] <- data[, r_index]
  }else{
    data0 <- data
  }
  if(length(varhat)!=3) stop("the length of varhat is not 3")
  if(is.null(forced_num)){
    yes <- 0
    num <- 0
  }else{
    yes <- 1
  }
  output <- cmasal(data0, maxterm, maxinteraction, maxiter, varhat, yes, num)
  rand_effects <- data.frame(U = output$rand_effects[, 1],
                             V1 = output$rand_effects[, 2],
                             V2 = output$rand_effects[, 3])
  ntime <- output$ntime
  Mobs <- max(ntime)
  nobs <- nrow(data0)
  if(sum(varhat)) regression <- 1
  else regression <- 0
  fitted <- get_formula(output)
  output <- c(output, fitted = list(fitted), data = list(data0))
  cov_matrix <- sort_covariate(output)

  Sigma_inv <- matrix(0, nrow = sum(ntime), ncol = sum(ntime))
  Sigma_inv[1:ntime[1], 1:ntime[1]] <- solve(cov_matrix[[1]])
  for(i in 2:length(ntime)){
    begin <- sum(ntime[1:(i-1)])+1
    end <- sum(ntime[1:i])
    Sigma_inv[begin:end, begin:end] <- solve(cov_matrix[[i]])
  }
  S_inv <- chol(Sigma_inv)

  data_masaled <- process(output)

  data_const <- cbind(data_masaled[, 1], rep(1, nobs), data_masaled[, -1])
  colnames(data_const) <- c(colnames(data_masaled)[1], "Intercept", colnames(data_masaled)[-1])
  data_weighted <- as.data.frame(S_inv %*% as.matrix(data_const))
  out <- lm(y ~ . - 1, data_weighted)
  masaled.values <- drop(as.matrix(data_const[, -1]) %*% out$coefficients)
  masaled.residuals <- drop(data_const[, 1] - masaled.values)
  names(masaled.values) <- names(masaled.residuals) <- rownames(data)

  out$call <- NULL

  out <- c(out, fitted = list(fitted), masaled.values = list(masaled.values),
           masaled.residuals = list(masaled.residuals),
           sigma = list(output$sigma), cov_matrix = list(cov_matrix),
           regression = list(regression), call = list(cl), nsub = list(output$nsub),
           nobs = list(nobs), ntime = list(ntime), Mobs = list(Mobs),
           learning.data = list(data0))
  enforce_test(out, forced_num)
  class(out) <- c("masal", "lm")
  out
}
