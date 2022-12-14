\name{forecast.stree}
\alias{forecast.stree}

\concept{survival tree}
\title{
  Predict method for STREE Model Fits
}
\description{
  Predicted values based on RTREE model object.
}
\usage{
  \method{forecast}{stree}(object, newdata, ...)
}
\arguments{
  \item{object}{Object of class inheriting from \code{"stree"}}
  \item{newdata}{An optional data frame in which to look for variables with
  which to predict. If omitted, the fitted values are used.}
  \item{...}{Further arguments to be passed to \code{forecast.stree}}
}
\details{
  \code{forecast.stree} produces predicted values, obtained by evaluating
  the function \code{\link[macs:stree]{stree}} in the frame \code{newdata}.
}
\value{
  \item{predictor}{result of prediction for each row.}
  \item{node}{each value matches the ordinal number of terminal node.}
}
\seealso{
  The model fitting function \code{\link[macs:stree]{stree}}
}
\examples{
# load "macs" library
library("macs")
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
result <-  stree(data, family = "likelihood")
result

# predictions
forecast(result, data[1:10, 3:9])
}
