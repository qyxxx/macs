\name{kmgraph}

\alias{kmgraph}

\title{
  Kaplan-Meier curves plot
}

\description{
Function \code{kmgraph} is able to generate Kaplan Meier curves based on censored data.
}

\usage{
kmgraph(object, nodecode = 0)
}

\arguments{
\item{object}{a 'stree' S3 object, i.e. the result returned by stree function.}
\item{nodecode}{default value 0 will return all the Kaplan-Meier curves of leaf nodes in the same graph. Users can also choose to display certain curve of a node by changing the value to other positive interger. Specifically, the leaf nodes are ordered from up to bottom, from left to right. Nodecode equals to k will return the k-th leaf node's Kaplan Meier curve.}
}
\details{
This function is a visualization method for class "stree". Given a \code{stree} object, \code{kmgraph} will use the data by which the model is generated, together with the fitted model to draw Kaplan-Meier curves of all or single node in one plot.
}
\seealso{
  \code{\link[macs:stree]{stree}}
}
\examples{
test <- data.frame(runif(1000, 1, 10), rbinom(1000, 1, 0.5),
                 sample(1:4, 1000, replace = TRUE, prob = c(0.1, 0.3, 0.2, 0.4)),
                 sample(1:50, 1000, replace = TRUE),
                 sample(1:30, 1000, replace = TRUE),
                 sample(1:10, 1000, replace = TRUE),
                 sample(1:60, 1000, replace = TRUE),
                 sample(1:20, 1000, replace = TRUE),
                 sample(1:40, 1000, replace = TRUE))
test[,3] <- as.factor(test[,3])
object <- stree(test,"likelihood")
kmgraph(object)
}
