calculate_x <- function(i, ra, dt, x){# calculate abscissa of nnd nodes
  #####################################################################
  # i:  ordinal number of nodes                                       #
  # ra: range of ith node                                             #
  # dt: array of left child                                           #
  # x:  array of abscissa of nnd nodes(have already initialized to 0) #
  #####################################################################
  # if leaf
  if(dt[i]==0) return(x)
  else {
    # range ratio of left child of ith node
    ratio <- (children(dt[i]+1, dt)+1)/(children(dt[i]+1, dt)+children(dt[i]+2, dt)+2)
    x[dt[i]+1] <- ra[1]+diff(ra)*ratio/2# abscissa of left child of ith node
    x[dt[i]+2] <- ra[2]-diff(ra)*(1-ratio)/2# abscissa of right child of ith node
    # recursive
    x <- calculate_x(dt[i]+1, c(ra[1], ra[1]+diff(ra)*ratio), dt, x)
    x <- calculate_x(dt[i]+2, c(ra[1]+diff(ra)*ratio, ra[2]), dt, x)
    return(x)
  }
}
