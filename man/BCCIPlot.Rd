\name{BCCIPlot}
\alias{BCCIPlot}
\title{
Confidence Intervals for Bias Corrected Local Regression
}
\description{
Graphs of confidence interval estimates for bias and standard
deviation of in bias-corrected local polynomial regression curve estimates.
}
\usage{
BCCIPlot(data, k1=1, k2=2, h, h2, output, g, layout, incl.biasplot, plotdata)
}
\arguments{
  \item{data}{A data frame, whose first column must be the explanatory variable
and whose second column must be the response variable.
}
  \item{k1}{degree of local polynomial used in curve estimator.}
  \item{k2}{degree of local polynomial used in bias estimator.}
  \item{h}{bandwidth for regression estimator.}
  \item{h2}{bandwidth for bias estimator.}
  \item{output}{if TRUE, numeric output is printed to the console window.}
  \item{g}{the target function, if known (for use in simulations).}
  \item{layout}{if TRUE, a 2x1 layout of plots is sent to the graphics device.}
  \item{incl.biasplot}{if TRUE, the confidence intervals for the bias of the 
uncorrected estimate are plotted.}
  \item{plotdata}{if TRUE, the data points are plotted as a scatter plot.}
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
