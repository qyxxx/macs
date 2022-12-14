\name{plot.ctree}
\alias{plot.ctree}
\concept{classification tree}
\title{
  Plot function for a ctree object
}
\description{
  Draw a tree structure plot for \code{ctree} object. The plot includes main parameters of RTREE model.
}
\usage{
  \method{plot}{ctree}(x, shape = 1, main = "ctree", ...)
}
\arguments{
  \item{x}{\code{ctree} object, typically result of \code{\link[macs:tarv]{tarv}} and \code{\link[macs:rtree]{rtree}}.}
  \item{shape}{has two options: 1 or 2. Determine the shape of tree where '1' uses circle and square to denote nodes while '2' uses point to represent nodes.}
  \item{main}{title to \code{ctree} plot.}
  \item{...}{Further arguments to be passed to \code{plot.ctree}}
}

\details{
  This function is a method for the generic function plot, for objects of class \code{ctree}.
}

\references{
Zhang, H. and Singer, B. (1999), Recursive partitioning in the health sciences, Springer Verlag.
}

\seealso{
  \code{\link[macs:tarv]{tarv}}, \code{\link[macs:rtree]{rtree}}
}

\examples{
# load "macs" library
library("macs")

# construct RTREE model
set.seed(1234)
data <- data.frame(y = sample(1:3, 1000, replace = TRUE),
                   n = as.factor(sample(1:4, 1000, replace = TRUE,
                                        prob = c(0.1, 0.3, 0.2, 0.4))),
                   o1 = sample(1:50, 1000, replace = TRUE),
                   o2 = sample(1:30, 1000, replace = TRUE),
                   o3 = sample(1:10, 1000, replace = TRUE),
                   o4 = sample(1:60, 1000, replace = TRUE),
                   o5 = sample(1:20, 1000, replace = TRUE),
                   o6 = sample(1:40, 1000, replace = TRUE))
head(data)
result <- rtree(data, method = "entropy", alpha = 0.01, cost = c(1, 1, 1))
result

# plot result of rtree
plot(result, main = "example of ctree plot")
plot(result, shape = 2)
}
