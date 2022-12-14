get_formula <- function(model){
  finTerm <- model$finTerm
  beta <- model$beta
  varList <- model$varList
  beta1 <- round(beta, digits = 3)
  varList1 <- round(varList, digits = 3)
  nterm <- length(beta1)
  output <- as.character(beta1[2])
  if(nterm>2){
    for(i1 in 3:nterm){
      i <- finTerm[i1];
      if(beta[i1]<0) output <- paste(output, " - ", abs(beta1[i1]), sep = "")
      else output <- paste(output, " + ", beta1[i1], sep = "")
      k <- varList[i+1, 1]
      for (j in 1:k) {
        if(varList[i+1, 3*j]> 0.){
          output <- paste(output,"(x_{", varList1[i+1, 3*j-1], "}", sep = "")
          if(varList[i+1, 3*j+1] > 0)
            output <- paste(output,"-", varList1[i+1, 3*j+1], ")^+", sep = "")
          else
            output <- paste(output,"+", -varList1[i+1, 3*j+1], ")^+", sep = "")
        }
        else
          output <- paste(output,"x_{", varList1[i+1, 3*j-1], "}", sep = "")
      }
    }
  }
  else warning("no variable was selected")
  return(output)
}
