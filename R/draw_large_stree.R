draw_large_stree <- function(nnd, dt, pt, spv, spvl, ncases, death_catg, median, depthv, x, height){
  depth <- max(depthv)
  n <- 6/depth
  intervaly <- height/(depth-1)
  points(x[1], height, pch = 1, cex = n)
  text(x[1] - 3*n, height + 2*n, labels = spv[1], cex = n)
  text(x[1] + 4*n, height + 5*n, labels = ncases[1], cex = n)
  text(x[1] + 4*n, height + 2*n, labels = death_catg[1], cex = n)
  text(x[1], height - 4*n, labels = median[1], cex = 3*n/4, col = "red")
  for(i in 2:nnd){
    x_i <- x[i]
    y_i <- height-intervaly*(depthv[i]-1)
    if(dt[i]==0) {
      points(x_i, y_i, pch = 0, cex = n)
      text(x_i + 3*n, y_i, labels = ncases[i], cex = 3*n/4)
      text(x_i + 3*n, y_i - 2*n, labels = death_catg[i], cex = 3*n/4)
      text(x_i - n, y_i - 2*n, labels = median[i], cex = 3*n/4, col = "red")
    }
    else {
      points(x_i, y_i, pch = 1, cex = n)
      text(x_i - 3*n, y_i, labels = spv[i], cex = n)
      text(x_i + 4*n, y_i + 3*n, labels = ncases[i], cex = n)
      text(x_i + 4*n, y_i, labels = death_catg[i], cex = n)
      text(x_i, y_i - 4*n, labels = median[i], cex = 3*n/4, col = "red")
    }
    lines(c(x[pt[i]+1], x_i), c(y_i+intervaly, y_i))
    if(!sum(i-1==dt)) text((x_i+x[pt[i]+1])/2 + 3*n, y_i+intervaly/2 + 3*n/2,
                           labels = spvl[pt[i] + 1], cex = n)
  }
}
