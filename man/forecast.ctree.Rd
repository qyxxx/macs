\name{forecast.ctree}
\alias{forecast.ctree}

\concept{classification tree}

\title{
  Predict method for RTREE Model Fits
}

\description{
  Predicted values based on RTREE model object.
}

\usage{
  \method{forecast}{ctree}(object, newdata, ...)
}

\arguments{

\item{object}{Object of class inheriting from \code{"ctree"}}

\item{newdata}{An optional data frame in which to look for variables with
  which to predict. If omitted, the fitted values are used.}

\item{...}{Further arguments to be passed to \code{forecast.ctree}}
}

\details{
  \code{forecast.ctree} produces predicted values, obtained by evaluating
  the function \code{\link[macs:tarv]{tarv}} and \code{\link[macs:rtree]{rtree}} in the frame \code{newdata}.
}

\value{

\item{predictor}{result of prediction for each row.}

\item{node}{each value matches the ordinal number of terminal node.}
}

\seealso{
  The model fitting function \code{\link[macs:tarv]{tarv}} and \code{\link[macs:rtree]{rtree}}
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

# predictions
forecast(result, data[1:10, 2:8])
}
