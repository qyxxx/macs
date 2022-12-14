plot.stree <- function(x, shape = 1, main = "stree", ...){
  # draw a tree
  #####################################################################
  # x:    output of cstree(a list)                            #
  # nnd:          number of nodes                                     #
  # dt:           array of left child                                 #
  # pt:           array of parent                                     #
  # spv:          array of categorical variables                      #
  # spvl:         array of threshold                                  #
  # ncases:       number of samples contained in the node             #
  # death_catg:
  # median:
  # colname:
  #####################################################################
  if(!inherits(x,"stree")) stop("Not a legitimate \" stree \" object")
  nnd <- x$nnd
  dt <- x$dt
  pt <- x$pt
  colname <- x$colname
  spv <- x$spv
  spvnum <- as.numeric(substring(x$spv, 2, nchar(x$spv)))
  spvl <- x$spvl

  threshold <- NULL
  for(i in 1:nnd){
    if(!is.na(spvl[, i])){
      if(colname[spvnum[i] + 2]>1) threshold <- c(threshold, spvl[, i])
      else threshold <- c(threshold, paste(">", round(spvl[, i], 2), sep = ""))
    }
    else threshold <- c(threshold, spvl[, i])
  }
  ncases <- x$ncases
  death_catg <- x$death_catg
  median <- round(x$median, 1)
  colname <- x$colname
  if(shape == 1){
    layer = c(1,2,2)
    if(nnd>3){
      for (i in 4:nnd) {
        count = 2
        j = i
        while (!(pt[j] == 0)) {
          j <- pt[j]+1
          count <- count +1
        }
        layer [i] <- count
      }
    }
    total <- layer[nnd]
    #Compute the index of leaf and non-leaf
    k<-1
    p<-1
    leafindex = c()
    notleafindex = c()
    for (i in 1:nnd) {
      if(dt[i] == 0){
        leafindex[k] = i-1
        k=k+1
      }
      else{
        notleafindex[p] = i-1
        p=p+1
      }
    }
    #Compute weight for each node
    weight = rep(0,nnd)
    weightfun <- function(x){
      if(dt[x+1]==0){
        return(3.2)
      }
      else{
        return(weightfun(dt[x+1])+weightfun(dt[x+1]+1))
      }
    }
    for (i in 1:length(leafindex)) {
      weight[leafindex[i]+1] = 1
    }
    for(i in 1:length(notleafindex)){
      weight[notleafindex[i]+1] = weightfun(notleafindex[i])
    }
    #edgelength indicates the length between brothers
    edgelength = c()
    for (i in 1:nnd) {
      edgelength[i] = weight[pt[i]+1]/2
    }
    #compute the coor of x and y
    coorx=c()
    coory=c()
    for (i in 1:nnd) {
      if(i==1){
        coorx[i] = 50
        coory[i] = (total-layer[i]+1)*10-5
      }
      else if(dt[pt[i]+1]+1 == i){
        coorx[i] = coorx[pt[i]+1] - edgelength[i]
        coory[i] = (total-layer[i]+1)*10-5
      }
      else{
        coorx[i] = coorx[pt[i]+1] + edgelength[i]
        coory[i] = (total-layer[i]+1)*10-5
      }
    }
    #get the index of children
    finddes <- function(number){
      des = c(number+1)
      k=2
      for (i in (number+1):nnd) {
        j = i
        while (!(pt[j] == number)) {
          j = pt[j]+1
          if(j == 1) break
        }
        if(!(j==1)){
          des[k] = i
          k=k+1
        }
      }
      return(des - 1)
    }
    #to judge whether Node Number's two sub-tree ovelap
    adjust <- function(number){
      leftdescent <- finddes(dt[number+1])
      rightdescent <- finddes(dt[number+1]+1)
      layerofleftdes = c()
      layerofrightdes = c()
      for(i in 1:length(leftdescent)){
        layerofleftdes[i] = layer[leftdescent[i]+1]
      }
      for(i in 1:length(rightdescent)){
        layerofrightdes[i] = layer[rightdescent[i]+1]
      }
      #Get the indexs of layers they both have
      index = intersect(layerofleftdes,layerofrightdes)
      #if(is.na(index) == TRUE){
      #  dd = c(-10000000)
      #  return(dd)
      #}
      dd=c()
      for (i in 1:length(index)) {
        k=1
        cr = c()
        for (j in 1:length(rightdescent)) {
          if(layer[rightdescent[j]+1] == index[i]) {
            cr[k] = coorx[rightdescent[j]+1]
            k=k+1
          }
        }
        k=1
        cl=c()
        for (j in 1:length(leftdescent)) {
          if(layer[leftdescent[j]+1] == index[i]) {
            cl[k] = coorx[leftdescent[j]+1]
            k=k+1
          }
        }
        dd[i] = max(cl) - min(cr)
      }
      return(dd)
    }
    judge=-100
    while (!(judge==0))
    {
      primaryx = coorx
      primaryy = coory
      for (i in 1:length(notleafindex)) {
        dd = adjust(notleafindex[i])
        node = notleafindex[i]
        j = node+1
        if(max(dd)>-6){
          edgelength[dt[node+1]+1] = edgelength[dt[node+1]+1] + max(abs(dd))/2 + 1
          edgelength[dt[node+1]+2] = edgelength[dt[node+1]+2] + max(abs(dd))/2 + 1
          while(!(pt[j]==0)){
            j = pt[j]+1
            edgelength[dt[j]+1] = edgelength[dt[j]+1] + max(abs(dd))/2 + 1
            edgelength[dt[j]+2] = edgelength[dt[j]+2] + max(abs(dd))/2 + 1
          }
          edgelength[2] = edgelength[2] + max(abs(dd))/2 + 1
          edgelength[3] = edgelength[3] + max(abs(dd))/2 + 1
        }
      }
      for (i in 1:length(notleafindex)) {
        node = notleafindex[i]
        j = node+1
        dd = adjust(notleafindex[i])
        if(max(dd)+7<0 && max(dd)>-10000000){
          edgelength[dt[j]+1] =edgelength[dt[j]+1]-(-7-max(dd))/2
          edgelength[dt[j]+2] =edgelength[dt[j]+2]-(-7-max(dd))/2
        }
      }
      coorx=c()
      coory=c()
      for (i in 1:nnd) {
        if(i==1){
          coorx[i] = 50
          coory[i] = (total-layer[i]+1)*10-5
        }
        else if(dt[pt[i]+1]+1 == i){
          coorx[i] = coorx[pt[i]+1] - edgelength[i]
          coory[i] = (total-layer[i]+1)*10-5
        }
        else{
          coorx[i] = coorx[pt[i]+1] + edgelength[i]
          coory[i] = (total-layer[i]+1)*10-5
        }
      }

      judge = sum(abs(primaryx-coorx)+abs(primaryy-coory))
    }
    i=1
    ##################################################################
    coorx=c()
    coory=c()
    for (i in 1:nnd) {
      if(i==1){
        coorx[i] = 50
        coory[i] = (total-layer[i]+1)*10-5
      }
      else if(dt[pt[i]+1]+1 == i){
        coorx[i] = coorx[pt[i]+1] - edgelength[i]
        coory[i] = (total-layer[i]+1)*10-5
      }
      else{
        coorx[i] = coorx[pt[i]+1] + edgelength[i]
        coory[i] = (total-layer[i]+1)*10-5
      }
    }
    #Get the scale of picture
    height = total * 10
    i=1
    while (!(dt[i]==0)) i = dt[i]+1
    left = i
    i=1
    while (!(dt[i]==0)) i = dt[i]+2
    right = i
    #################
    h=6/total
    coorx<-(coorx-50)

    par_1 <- par(bg = "white", xaxt = "n", yaxt = "n", bty = "n", cex = 1,
        ann = FALSE, mar = c(0, 0, 0, 0), oma = c(0, 0, 2, 0))
    on.exit(par(par_1))
    plot(300,300,xlim = c(coorx[left]-5,coorx[right]+5),ylim = c(-5,total*10),axes =FALSE,ann = F,cex=1,asp=1)
    t<-seq(from=0,to=(2*pi),len=1000)
    for(i in 1:nnd)
    {
      if(!(dt[i]==0)) lines(3*sin(t)+coorx[i],3*cos(t)+coory[i],type="l")
      else rect(coorx[i]-3,coory[i]-3,coorx[i]+3,coory[i]+3)
      if(dt[i]!=0)
      {
        arrows(coorx[i]-1.5,coory[i]-1.5*sqrt(3),coorx[dt[i]+1],coory[dt[i]+1]+3,length = 0.05*h)
        arrows(coorx[i]+1.5,coory[i]-1.5*sqrt(3),coorx[dt[i]+2],coory[dt[i]+1]+3,length = 0.05*h)
      }
    }
    fcount = paste(ncases, "\n", death_catg, sep = "")
    for (i in 1:nnd) {
      text(coorx[i], coory[i], fcount[i],cex= h)
      text(coorx[i], coory[i]-3.6, median[i],cex= 0.8*h, col = "red")
      if(!(dt[i]==0)){
        text(coorx[i], (coory[i])-5, spv[i], cex= 0.8*h)
        text((coorx[i]+coorx[dt[i]+2]+1.5)/2, (coory[i]+coory[dt[i]+2]+(1-sqrt(3)/2)*3)/2,threshold[i], cex=h)
      }
    }
  }
  else{
    depthv <- 1
    for(i in 2:nnd){
      d <- 1
      j <- i
      while(pt[j]!=0){
        d <- d + 1
        j <- pt[j] + 1
      }
      depthv <- c(depthv, d + 1)
    }
    height <- 100
    x <- rep(0, nnd)
    x <- calculate_x(1, c(0, height), dt, x)
    x[1] <- height/2
    for(i in seq(nnd, 3, -2)){
      x[pt[i]+1] <- (x[i]+x[i-1])/2
    }

    par_2 <- par(bg = "white", xaxt = "n", yaxt = "n", bty = "n", cex = 1,
        ann = FALSE, mar = c(0, 0, 0, 0), oma = c(0, 0, 2, 0))# c(bottom, left, top, right)
    on.exit(par(par_2))
    plot(0, xlim = c(0, height), ylim=c(-5, height + 5), type = "n", asp = 1)
    draw_large_stree(nnd, dt, pt, spv, threshold, ncases, death_catg, median, depthv, x, height)
  }
  title(main, cex.main = 1.5, font.main = 1, outer = TRUE)
}

