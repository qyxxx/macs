draw_large_ctree <- function(nnd, dt, pt, spv, spvl, final_counts, nodeclass, depthv, x, height){
  depth <- max(depthv)
  n <- 8/depth
  intervaly <- height/(depth-1)
  points(x[1], height, pch = 1, cex = n)
  text(x[1] - 3*n, height, labels = spv[1], cex = n)
  text(x[1] + 4*n, height + 3*n, labels = final_counts[1, 1], cex = n)
  text(x[1] + 4*n, height, labels = final_counts[1, 2], cex = n)
  for(i in 2:nnd){
    x_i <- x[i]
    y_i <- height-intervaly*(depthv[i]-1)
    if(dt[i]==0) {
      points(x_i, y_i, pch = 0, cex = n)
      text(x_i, y_i - 2*n, labels = final_counts[i, 1], cex = 0.8*n)
      text(x_i, y_i - 4*n, labels = final_counts[i, 2], cex = 0.8*n)
      text(x_i, y_i - 8*n, labels = nodeclass[i], cex = 1.2*n, col = "red")
    }
    else {
      points(x_i, y_i, pch = 1, cex = n)
      text(x_i - 3*n, y_i, labels = spv[i], cex = n)
      text(x_i + 4*n, y_i + 3*n, labels = final_counts[i, 1], cex = n)
      text(x_i + 4*n, y_i, labels = final_counts[i, 2], cex = n)
    }
    lines(c(x[pt[i]+1], x_i), c(y_i+intervaly, y_i))
    if(!sum(i-1==dt)) text((x_i+x[pt[i]+1])/2 + 3*n, y_i+intervaly/2 + 3*n/2,
                           labels = spvl[pt[i] + 1], cex = n)
  }
}
