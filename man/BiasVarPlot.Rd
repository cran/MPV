\name{BiasVarPlot}
\alias{BiasVarPlot}
\title{
Local Polynomial Bias and Variability 
}
\description{
Graphs of confidence interval estimates for bias and standard
deviation of in local polynomial regression curve estimates.
}
\usage{
BiasVarPlot(data, k1=1, k2=2, h, h2, output=FALSE, g, layout=TRUE)
}
\arguments{
  \item{data}{A data frame, whose first column must be the explanatory variable
and whose second column must be the response variable.
}
  \item{k1}{degree of local polynomial used in curve estimator.}
  \item{k2}{degree of local polynomial used in bias estimator.}
  \item{h}{bandwidth for regression estimator.}
  \item{h2}{bandwidth for bias estimator.}
  \item{output}{if true, numeric output is printed to the console window.}
  \item{g}{the target function, if known (for use in simulations).}
  \item{layout}{if true, a 2x1 layout of plots is sent to the graphics device.}
}
\value{A list containing the confidence interval limits, pointwise 
estimates of bias, standard deviation of bias, curve estimate, standard
deviation of curve estimate, and approximate confidence limits for curve
estimates. Graphs of the curve estimate confidence limits and the bias
confidence limits. }
\author{
W. John Braun and Wenkai Ma
}
\keyword{graphics}
