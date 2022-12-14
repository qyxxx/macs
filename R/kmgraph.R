kmgraph <- function(object, nodecode = 0){
  median <- object$median
  ncases <- object$ncases
  nnd <- object$nnd
  dt <- object$dt
  pt <- object$pt
  data <- object$learning.data
  leafindex <- which(dt==0)
  node <- forecast.stree(object)$node
  pro <- list()
  xlimright <- NULL
  xlimleft <- NULL
  numberofcol = ceiling(length(leafindex)/18)
  for (k in 1:length(leafindex)) {
    datatemp <- data[which(node == leafindex[k]),]
    censordata <- datatemp[which(datatemp[,2]==1),]
    censor <- sort(unique(censordata[,1]))
    if(length(censor)!=0){
      xlimright <- c(xlimright, max(censor))
      xlimleft <- c(xlimleft, min(censor))
      y <- datatemp[,1]
      d <- c()
      for (i in 1:length(censor)) {
        d[i] <- freq(censordata[,1],censor[i])

      }
      riskset <- rep(0,length(censor))
      for (i in 1:length(censor)) {
        for (j in 1:nrow(datatemp)) {
          if(y[j]>=censor[i])
            riskset[i] <- riskset[i] + 1
        }
      }
      ratio <- 1-(d/riskset)
      product <- c()
      for (i in 1:length(censor)) {
        if(i == 1) product[i] <- ratio[i]
        else product[i] <- product[i-1]*ratio[i]
      }
    }
    else product <- 0
    pro <- c(pro,list(product))
  }
  par_1 <- par(bg = "white", xaxt = "s", yaxt = "s", xaxs = "i", yaxs = "i", bty = "n", ann = TRUE,
      mgp = c(3, 1, 0), las = 1, cex = 1,
      mar = c(5.1, 4.1, 4.1, 2.1), oma = c(0, 0, 0, 0))
  on.exit(par(par_1))
  rb <- rainbow(length(leafindex))
  if(!nodecode){
    xleft <- floor(min(xlimleft))
    xright <- ceiling(max(xlimright))
    for (k in 1:length(leafindex)) {
      if(k<(length(pro)+1))
        product <- pro[[k]]
      datatemp <- data[which(node == leafindex[k]),]
      censordata <- datatemp[which(datatemp[,2]==1),]
      censor <- sort(unique(censordata[,1]))
      if(length(censor)!=0){
        color <- rb[k]
        if(k==1){
          plot(censor[1:length(censor)], product[1:length(censor)],
               xlim = c(xleft,xright+1+numberofcol), ylim = c(0,1),xlab = "Time",
               ylab = "Survival", main = "Kaplan - Meier curve",
               type = "n")
        }
        lines(c(0,censor[1]),c(1,1),col = color)
        lines(c(censor[1],censor[1]),c(1,product[1]),col = color)
        if(length(censor)>1){
          for (i in 1:(length(censor)-1)) {
            lines(c(censor[i],censor[i+1]),c(product[i],product[i]),col = color)
            lines(c(censor[i+1],censor[i+1]),c(product[i],product[i+1]),col= color)
          }
        }
        lines(c(censor[length(censor)],xright),c(product[length(censor)],product[length(censor)]),col = color)
      }
    }
    legend(xright+0.2,1,c(1:length(leafindex)),lty = rep(1,length(leafindex)), lwd =  rep(1,length(leafindex)),col = rb,ncol = numberofcol,cex = 0.5)
  }
  else if(nodecode < (1+length(leafindex))){
    if(sum(pro[[nodecode]]) == 0)
    {
      stop("this node has no censored data")
    }
    x1 <- xlimleft[nodecode]
    x2 <- xlimright[nodecode]
    product <- pro[[nodecode]]
    datatemp <- data[which(node == leafindex[nodecode]),]
    censordata <- datatemp[which(datatemp[,2]==1),]
    censor <- sort(unique(censordata[,1]))
    if(length(censor)!=0){
      color <- rb[nodecode]
      #par(no.readonly = FALSE)
      par_2 <- par(mar = c(4,4,2,2))
      on.exit(par(par_2))
      if(x1 == x2) {
        x2 <- 2*x1 - floor(x1)
        x1 <- floor(x1)
      }
      plot(censor[1:length(censor)],product[1:length(censor)],
           xlim = c(x1,x2),ylim = c(0,1),xlab = "Time",
           ylab = "Survival", main = "Kaplan - Meier curve",
           type = "n")
      lines(c(0,censor[1]),c(1,1),col = color)
      lines(c(censor[1],censor[1]),c(1,product[1]),col = color)
      if(length(censor)>1){
        for (i in 1:(length(censor)-1)) {
          lines(c(censor[i],censor[i+1]),c(product[i],product[i]),col = color)
          lines(c(censor[i+1],censor[i+1]),c(product[i],product[i+1]),col= color)
        }
      }
      lines(c(censor[length(censor)],x2),c(product[length(censor)],product[length(censor)]),col = color)
    }
  }
  else{
    stop("code number is wrong")
  }
}

freq <- function(vector,number){
  j = 0
  for (i in 1:length(vector)) {
    if(vector[i] == number)
      j = j+1
  }
  return(j)
}
