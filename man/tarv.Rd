\name{tarv}
\alias{tarv}
\title{
  Recursive partitioning based classification trees for risk profile and diagnosis
}
\description{
  Fit a \code{ctree} model
}
\usage{
tarv(pheno, geno, formula, method = c("entropy", "gini"),
                  family = c("binomial", "gaussian"),
                  direction = c("both", "positive", "negative"),
                  alpha = 0.01, cost = NULL)
}
\arguments{
  \item{pheno}{defines the phenotypic data, serving as the response variable.}
  \item{geno}{defines the genetic variants.}
  \item{formula}{specifies the model structure such as '\code{disease~sex+race}'.}
  \item{method}{is the same as that in \code{\link[macs:rtree]{rtree}}.}
  \item{family}{specifies the regression method: '\code{binomial}' for the logistic regression for a binary phenotype and '\code{gaussian}' for the linear regression with a normally distributed phenotype.}
  \item{direction}{offers the option for the inclusion of genetic markers. When \code{tarv} processes the original genetic markers, it may rank their univariate effects according to their t-values in three ways, namely t (\code{positive}), -t (\code{negative}), and |t| (\code{both}). We recommend \code{both} to avoid any potential confusion which means that more statistically significant genetic markers would be considered favorably regardless of their directions of the effect on the phenotype.}
  \item{alpha}{is the same as that in \code{\link[macs:rtree]{rtree}}.}
  \item{cost}{is the same as that in \code{\link[macs:rtree]{rtree}}.}
}

\details{
'pheno' and 'geno' are processed via TRAV[Song and Zhang 2014] to create new covariates and then these new virables are used to create a tree construction. So it has a similar input and output form with rtree. The object derives from rtree and ctree function are both 'ctree' object, and can be used to plot.ctree and predict.ctree.
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
Zhang, H. and Singer, B. (1999), Recursive partitioning in the health sciences, Springer Verlag.

Song C. and Zhang H.(2014), Tree-based Analysis of Rare Variants Identifying Risk Modifying Variants in CTNNA2 and CNTNAP2 for Alcohol Addiction.
}

\seealso{
  \code{\link[macs:plot.ctree]{plot.ctree}}, \code{\link[macs:forecast.ctree]{forecast.ctree}}, \code{\link[macs:rtree]{rtree}}
}

\examples{
library(macs)
set.seed(123)
sex <- rbinom(1000, 1, 0.5)
race <- sample(1:3, 1000, replace = TRUE, prob = c(0.4, 0.4, 0.2))
gg <- replicate(100, rbinom(1000, 1, runif(1, 0.005, 0.05)))
annotation <- paste("gene", rep(1:5, each = 20), sep = "")
causal <- rbinom(40, 1, 0.8)
x1 <- rowSums(gg[, 1:20][, causal[1:20] > 0]) > 0
x2 <- rowSums(gg[, 21:40][, causal[21:40] > 0]) > 0
xb <- sex * 0.2 + (race == 2) * 0.2 + x1 * 0.6 + x2 * 0.8 - 0.7
r <- rbinom(1000, 1, exp(xb)/(1+exp(xb)))
pheno <- data.frame(disease = r, sex = sex, race = as.factor(race))
geno <- t(gg)
rownames(geno) <- annotation

result <- tarv(pheno, geno, formula = "disease~sex+race",
                 method = "entropy", family = "binomial",
                 direction = "both", alpha = 0.01, cost = c(1, 1))
plot(result)
}

