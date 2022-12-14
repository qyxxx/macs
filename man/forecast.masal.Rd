\name{forecast.masal}
\alias{forecast.masal}

\concept{regression}
\title{
  Predict method for MASAL Model Fits
}
\description{
  Predicted values based on MASAL model object.
}
\usage{
  \method{forecast}{masal}(object, newdata = NULL, col = NULL, ...)
}
\arguments{
\item{object}{Object of class inheriting from \code{"masal"}}
\item{newdata}{An optional data frame in which to look for variables with
    which to predict. If omitted, the fitted values are used.}
\item{col}{A character vector indicating each column in the \code{data} : \code{c}(covariate), \code{t}(time), \code{s}(sort). If omitted, the fitted values are used.}
\item{...}{Further arguments to be passed to \code{forecast.masal}}
}
\details{
  \code{forecast.masal} produces predicted values, obtained by evaluating
  the regression function in the frame \code{newdata}.
}
\value{
  \code{forecast.masal} produces a vector of predictions for all rows.
}
\seealso{
  The model fitting function \code{\link[macs:masal]{masal}}.
}

\examples{
\donttest{
library("macs")
set.seed(123)
data <- data.frame(ID = sort(rep(1:100, 5)),
                    X1 = runif(500, min = 0, max = 1),
                    X2 = runif(500, min = 0, max = 2),
                    X3 = runif(500, min = 0, max = 4),
                    X4 = rnorm(500, mean = 0, sd = 1),
                    X5 = rnorm(500, mean = 0, sd = 2),
                    X6 = rnorm(500, mean = 0, sd = 3),
                    time = rep(1:5, 100),
                    Y = rnorm(500, mean = 40, sd = 10))
result <- masal(data, col = NULL, maxterm = 20, maxinteraction = 2,
 varhat = c(1, 1, 1), forced_num = NULL, maxiter = 2)
forecast(result)
}
}
