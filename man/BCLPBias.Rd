\name{BCLPBias}
\alias{BCLPBias}
\title{
Bias for Bias-Corrected Local Polynomial Regression
}
\description{
Confidence interval estimates for bias in local polynomial regression.
}
\usage{
BCLPBias(xy,k1,k2,h,h2,numgrid=401,alpha=.95)
}
\arguments{
  \item{xy}{A data frame, whose first column must be the explanatory variable
and whose second column must be the response variable.
}
  \item{k1}{degree of local polynomial used in curve estimator.}
  \item{k2}{degree of local polynomial used in bias estimator.}
  \item{h}{bandwidth for regression estimator.}
  \item{h2}{bandwidth for bias estimator.}
  \item{numgrid}{number of gridpoints used in the curve estimator.}
  \item{alpha}{nominal confidence level.}
}
\value{A list containing the confidence interval limits, pointwise 
estimates of bias, standard deviation of bias, curve estimate, standard
deviation of curve estimate, and approximate confidence limits for curve
estimates and corresponding bias-corrected estimates. }
\author{
W. John Braun and Wenkai Ma
}
\keyword{models}
