trim_ctree <- function(object){
  nnd <- object$nnd
  dt <- object$dt[1:nnd]
  pt <- object$pt[1:nnd]
  spv <- object$spv[1:nnd]
  spvl <- object$spvl[1:nnd]
  final_counts <- object$final_counts[1:nnd, ]
  varcatg <- object$varcatg

  spvl[dt==0] <- NA
  threshold <- data.frame(t(spvl))
  for(i in 1:nnd){
    if(varcatg[spv[i] + 2]>1 & !is.na(spvl[i])){
      sp <- as.character(spvl[i])
      if(nchar(sp)<=1) threshold[, i] <- paste0("{", sp, "}")
      else {
        single <- substring(sp, 1:nchar(sp), 1:nchar(sp))
        threshold[, i] <- paste0("{", single[1])
        for(j in 2:nchar(sp)) threshold[, i] <- paste(threshold[, i], single[j], sep = ",")
        threshold[, i] <- paste0(threshold[, i], "}")
      }
    }
  }

  variable <- paste("x", spv + 1, sep = "")
  variable[dt==0] <- NA
  return(list(nnd = nnd, dt = dt, pt = pt, spv = variable, spvl = threshold,
              final_counts = final_counts, varcatg = varcatg))
}
