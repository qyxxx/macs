\name{plot.stree}
\alias{plot.stree}
\concept{survival tree}
\title{
  Plot function for a stree object
}
\description{
  Draw a tree structure plot for \code{stree} object. The plot includes main parameters of STREE model.
}
\usage{
  \method{plot}{stree}(x, shape = 1, main = "stree", ...)
}
\arguments{
  \item{x}{\code{stree} object, typically result of \code{\link[macs:stree]{stree}}.}
  \item{shape}{has two options: 1 or 2. Determine the shape of tree where '1' uses circle and square to denote nodes while '2' uses point to represent nodes.}
  \item{main}{title to \code{stree} plot.}
  \item{...}{Further arguments to be passed to \code{plot.stree}}
}

\details{
  This function is a method for the generic function plot, for objects of class \code{stree}.
}

\references{
Zhang, H. and Singer, B. (1999), Recursive partitioning in the health sciences, Springer Verlag.
}

\seealso{
  \code{\link[macs:stree]{stree}}
}

\examples{
# load "macs" library
library("macs")

# construct STREE model
set.seed(1234)
data <- data.frame(r = runif(1000, 1, 10), c = rbinom(1000, 1, 0.5),
                   n = as.factor(sample(1:4, 1000, replace = TRUE,
                                        prob = c(0.1, 0.3, 0.2, 0.4))),
                   o1 = sample(1:50, 1000, replace = TRUE),
                   o2 = sample(1:30, 1000, replace = TRUE),
                   o3 = sample(1:10, 1000, replace = TRUE),
                   o4 = sample(1:60, 1000, replace = TRUE),
                   o5 = sample(1:20, 1000, replace = TRUE),
                   o6 = sample(1:40, 1000, replace = TRUE))
head(data)
result <- stree(data, family = "likelihood")

# plot result of stree
plot(result, main = "example of stree plot")
plot(result, shape = 2)
}
