\name{rtree}
\alias{rtree}
\title{
  Tree-based analysis of rare variants to analyze genomic data.
}
\description{
  Fit a \code{rtree} model
}

\usage{
rtree(data, method = c("entropy", "gini"), alpha = 0.01,cost = NULL)
}
\arguments{
  \item{data}{inputs \code{data.frame}. The first column of the data refers to the response  and the other columns are the covariates.}
  \item{method}{defines the splitting rule, that is, 'entropy' for the entropy criterion and 'gini' for the Gini criterion.}
  \item{alpha}{determines the significance threshold by which \code{\link[macs:rtree]{rtree}} prunes the tree based on the chi-square test of independence between the distributions of the response in any two daughter nodes.}
  \item{cost}{defines misclassification cost of each category. It is a vector whose dimension is the same as the levels of response variable. Default value is a vector with all elements 1.}
}

\details{
The object rtree function produces can be used to run \code{\link[macs:plot.ctree]{plot.ctree}} and \code{\link[macs:forecast.ctree]{forecast.ctree}} function.
}

\value{
  \item{nnd}{the total number of nodes in the tree.}
  \item{dt}{
  the sequence number of a left daughter node for each internal node.
  }
  \item{pt}{
  the sequence number of the parent node for any daughter node.
  }
  \item{spv}{
  the splitting variable used to split a given node.
  }
  \item{spvl}{
  the cut-off value of the splitting variable above.
  }
  \item{final_counts}{
  the table that contains the number of observations in each node.
  }
  \item{varcatg}{
  a numerical indicator for the category of each variable. Value '-1' points to the response variable, '1' to oridinal variables, an integer greater than 1 to a nominal variable with the number of levels equal to the integer.
  }
  \item{nodeclass}{
  the class membership of a terminal node which depends on the choice of the misclassification cost.
  }
  \item{p_value}{the p-value of the chi-square test performed at each internal node. It forms the basis to prune the offspring nodes of any internal node. More details in Recursive Partitioning and Applications [Zhang and Singer].
  }
  \item{call}{the call by which this object is generated.}
  \item{learning.data}{the data that are actually used in \code{\link[macs:rtree]{rtree}} and \code{\link[macs:tarv]{tarv}}.}
}

\references{
Zhang, H. and Singer, B. (1999) \emph{Recursive partitioning in the health sciences.} Springer Verlag.
}

\seealso{
  \code{\link[macs:plot.ctree]{plot.ctree}}, \code{\link[macs:forecast.ctree]{forecast.ctree}}
}

\examples{
library(macs)
set.seed(1234)
data <- data.frame(y = sample(1:3, 1000, replace = TRUE),
                   n = sample(1:4, 1000, replace = TRUE,
                                        prob = c(0.1, 0.3, 0.2, 0.4)),
                   o1 = sample(1:50, 1000, replace = TRUE),
                   o2 = sample(1:30, 1000, replace = TRUE),
                   o3 = sample(1:10, 1000, replace = TRUE),
                   o4 = sample(1:60, 1000, replace = TRUE),
                   o5 = sample(1:20, 1000, replace = TRUE),
                   o6 = sample(1:40, 1000, replace = TRUE))
data[,2] <- as.factor(data[,2])
result <- rtree(data, method = "entropy", alpha = 0.01, cost = NULL)
plot(result)
}

