\name{masal}
\alias{masal}
\title{Recursive partitioning based masal}

\description{
  Fit the regression function based on multivariate adaptive splines.
}

\usage{
  masal(data, col = NULL, maxterm = 20, maxinteraction = 2,
                varhat = c(1, 1, 1), forced_num = NULL, maxiter = 2)
}

\arguments{
\item{data}{a data.frame including one column of ID of the observations, covariates, time, and response.
}
\item{col}{a charater vector indicating each column in the \code{data}: \code{c}(covariate), \code{i}(ID), \code{r}(response), \code{t}(time). If NULL, column 1 of \code{data} should be the ID number, columns 2 to the third last from the end should be the covariates, the second last from the end should be time, and the response should be in the last column. The default is \code{NULL}.
}
\item{maxterm}{the maximum number of terms allowed in the model during the forward step. The default is \code{20}.
}
\item{maxinteraction}{the maximum number of interactive predictors allowed in the model. The default is \code{2}.
}
\item{varhat}{a vector of length three indicating whether any of the three random effects should be considered. A value of `1' or `0' indicates whether the corresponding term should be included or not. The default is \code{c(1, 1, 1)}.
}
\item{forced_num}{the number of the covariates that a user wishes to force into the \code{masal}. The default is \code{NULL}.
}
\item{maxiter}{the maximum number of iterations between the surface estimation and covariance estimation. The default is 2.}
}
\details{
For cross-sectional data, function \code{masal} fits a nonparametric surface for a continuous response variable. For longitudinal data, \code{masal} provides a surface that depends on time and selected covariates, and their interactions are allowed.  The regression surface is based on the multivariate adaptive splines.

The data should include the identification (ID) number that links possibly multiple observations from the same unit of observation, the covariates, the time of each observation, and the response variables. For categorical covariates, dummy variables should be created before input.
}
\value{
\item{fitted}{the fitted fixed-effect function f.
}
\item{masaled.values}{the fitted values with the formula obtained by \code{masal}.
}
\item{masaled.residuals}{the residuals from \code{masal}.
}
\item{sigma}{a vector of 4 containing the estimated parameters related to the random effects.}
\item{cov_matrix}{the estimated covariance matrix depending on \code{sigma}.
}
\item{regression}{an indicator for whether some of the three random effects in the model has been included. A value of `1' means yes and `0' otherwise.
}
\item{call}{the call by which this object is generated.
}
\item{nsub}{the number of subjects contained by input \code{data}.
}
\item{nobs}{the number of observations, i.e., the number of rows of \code{data}.}

\item{ntime}{a vector that indicates the number of observations for each subject of \code{data}.
}
\item{Mobs}{the maximum number of \code{ntime}.
}
\item{learning.data}{the data that are actually used in \code{masal}. They are not necessarily the same as original \code{data} since some processing may be conducted.}
}

\references{
Zhang, H. (1997) \emph{Multivariate adaptive splines for analysis of longitudinal data}.

Zhang, H. and Singer, B. (1999) \emph{Recursive partitioning in the health sciences.} Springer Verlag.
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
plot(result)
forecast(result)
summary(result)
}
}