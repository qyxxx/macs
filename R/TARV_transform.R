TARV_transform <- function(geno, annotation, Z, direction = c("both", "positive", "negative"))
{
  direction <- match.arg(direction)
  if(nrow(geno) != length(annotation)) stop("Error: nrow(geno) != length(annotation)\n")
  if(nrow(geno) != length(Z)) stop("Error: nrow(geno) != length(Z)\n")
  gid <- unique(annotation)
  if(direction %in% c("positive", "both")) {
    X.pos <- matrix(NA, ncol(geno), length(gid))
    colnames(X.pos) <- paste(gid, "+", sep = "")
  }
  if(direction %in% c("negative", "both")) {
    X.neg <- matrix(NA, ncol(geno), length(gid))
    colnames(X.neg) <- paste(gid, "-", sep = "")
  }
  for(i in 1:length(gid)) {
    g <- gid[i]
    gindex <- annotation == g
    z <- Z[gindex]
    if(direction %in% c("positive", "both")) {
      order.pos <- order(-z)
      index.pos <- which(gindex)[order.pos]
      X.pos[, i] <- apply(rbind(geno[index.pos, ] > 0, TRUE), 2, function(x) which(x)[1])
    }
    if(direction %in% c("negative", "both")) {
      order.neg <- order(z)
      index.neg <- which(gindex)[order.neg]
      X.neg[, i] <- apply(rbind(geno[index.neg, ] > 0, TRUE), 2, function(x) which(x)[1])
    }
  }
  if(direction == "both") {
    X <- cbind(X.pos, X.neg)
  } else if(direction == "positive") {
    X <- X.pos
  } else if(direction == "negative") {
    X <- X.neg
  }
  return(X)
}
