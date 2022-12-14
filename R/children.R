children <- function(i, dt){# number of children of node i
  #####################################################################
  # i:  ordinal number of nodes                                       #
  # dt: array of left child                                           #
  #####################################################################

  if(dt[i]==0) return(0)
  else {
    # recursive
    return(children(dt[i]+1, dt)+children(dt[i]+2, dt)+2)
  }
}
