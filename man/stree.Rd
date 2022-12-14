\name{stree}
\alias{stree}
\title{
  Recursive partitioning based survival tree
}

\description{
  Fit a \code{stree} model for censored data
}
\usage{stree(data, family)
}

\arguments{

\item{data}{input data, including a time response, a censoring indicator (‘1’ for acensored time and ‘0’ for an observed time), and at least one nominal or ordinal variable, and typically more covariates of either type.}

\item{family}{the splitting rule for an internal node in a survival tree. It includes five choices: `likelihood', `log-rank', `kaplan-meier distance', `adaptive normalization', and `global normalization'}
}


\details{
\code{stree} uses five rules which can be chosen by users to create a survival tree based on censored data. Data format is required to include at least one column of response variable, one column of censored variable, one column of nominal or ordinal variable. Given data and chosen family, it will produce an 'stree' object that can be used for \code{\link[macs:plot.stree]{plot.stree}}, \code{\link[macs:kmgraph]{kmgraph}} and other further analysis like prediction. The 'stree' object will include necessary properties of the tree.
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

\item{ncases}{
  the number of observations in each node.
  }

\item{death_catg}{
  the number of deaths in each node.
  }

\item{median}{
  the median survival time in each node.
  }

\item{colname}{
  a numerical indicator for the category of each variable. Value ‘-2' points to the response variable, ‘1' to oridinal variables, ‘2' to a nominal variable, ‘-1' to censoring indicator.
  }

\item{call}{the call by which this object is generated.}

\item{learning.data}{the data that are actually used in \code{\link[macs:stree]{stree}}.}

}


\references{
 Zhang, H. and Singer, B. (1999) \emph{Recursive partitioning in the health sciences.} Springer Verlag.
}


\seealso{
  \code{\link[macs:plot.stree]{plot.stree}}, \code{\link[macs:forecast.stree]{forecast.stree}}, \code{\link[macs:kmgraph]{kmgraph}}}



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
 result <- stree(test, "likelihood")
plot(result)
}

