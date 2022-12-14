trim_stree <- function(object){
  nnd <- object$nnd
  dt <- object$dt
  pt <- object$pt
  maxnodes <- length(dt)
  spv <- object$spv
  spvl <- object$spvl
  ncases <- object$ncases
  death_catg <- object$death_catg
  median <- object$median
  colname <- object$colname
  dtemp <- dt
  dtemp[dtemp>0] <- seq(1, nnd-2, 2)
  i <- 2
  j <- 0
  while(i<=maxnodes) {
    if(dt[pt[i]+1]==0) { # avoid printing pruned node C-Y
      dtemp[(i-j):maxnodes] <- c(dtemp[(i-j+2):maxnodes], c(0, 0))
      spv[(i-j):maxnodes] <- c(spv[(i-j+2):maxnodes], c(0, 0))
      spvl[(i-j):maxnodes] <- c(spvl[(i-j+2):maxnodes], c(0, 0))
      ncases[(i-j):maxnodes] <- c(ncases[(i-j+2):maxnodes], c(0, 0))
      death_catg[(i-j):maxnodes] <- c(death_catg[(i-j+2):maxnodes], c(0, 0))
      median[(i-j):maxnodes] <- c(median[(i-j+2):maxnodes], c(0, 0))
      i <- i + 2
      j <- j + 2
    }
    else i <- i + 1
  }
  dt <- dtemp[1:nnd]
  pt <- pt[1:nnd]
  spv <- spv[1:nnd]
  spvl <- spvl[1:nnd]
  ncases <- ncases[1:nnd]
  death_catg <- death_catg[1:nnd]
  median <- median[1:nnd]
  for(i in 1:nnd){
    if(dt[i]!=0) {
      pt[dt[i]+1] <- i-1
      pt[dt[i]+2] <- i-1
    }
  }

  spvl[dt==0] <- NA
  threshold <- data.frame(t(spvl))
  for(i in 1:nnd){
    if(colname[spv[i] + 3]>1 & !is.na(spvl[i])){
      sp <- as.character(spvl[i])
      if(nchar(sp)<=1) threshold[, i] <- paste0("{", spvl[i]-1, "}")
      else {
        single <- as.character(as.numeric(substring(sp, 1:nchar(sp), 1:nchar(sp)))-1)
        threshold[, i] <- paste0("{", single[1])
        for(j in 2:nchar(sp)) threshold[, i] <- paste(threshold[, i], single[j], sep = ",")
        threshold[, i] <- paste0(threshold[, i], "}")
      }
    }
  }

  variable <- paste("x", spv + 1, sep = "")
  variable[dt==0] <- NA
  return(list(nnd = nnd, dt = dt, pt = pt, spv = variable, spvl = threshold,
              ncases = ncases, death_catg = death_catg, median = median,
              colname = colname))
}
