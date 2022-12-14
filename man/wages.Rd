\name{wages}
\alias{wages}
\docType{data}
\title{Wages of male dropouts from high school.}
\description{
The dataset originally comprised of 888 male dropouts from high school and their wages at different times. Our goal is to model the longitudinal data to track the wages of them. Each dropout has several observations. Exper is the duration between the observation and his first day. For example, dropout 1 has 8 observations; the first three data entries refer to the first two years of work (at 0.015, 0.715 and 1.734 years). Wages were listed on an hourly base. The wages dataset contains 11 potential explanatory variables. There are 123 dropouts who have 8 observations, and 127 have 9 observations. The maximum number of observations of dropouts is 13 while 38 dropouts only had 1 observation. These dropouts with 1 observation are deleted because they provide less information about the wage change of the person. So, the dataset in the package only comprised of 850 male dropouts.
}
\usage{data("wages")}
\format{A data.frame containing 6363 observations from 888 male dropouts.}
\references{
Singer, Judith D and Willett, John B and Willett, John B and others (2003) \emph{Applied longitudinal data analysis: Modeling change and event occurrence}. Oxford university press.
}
\keyword{datasets}
