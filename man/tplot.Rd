\name{tplot}
\alias{tplot}
\title{
Graphical t Test for Regression
}
\description{
This function analyzes regression data graphically.  It 
allows visualization of the usual t-tests for individual
regression coefficients.
}
\usage{
tplot(X, y, plotIt=TRUE, type="hist", includeIntercept=TRUE)
}
\arguments{
  \item{X}{The design matrix.
}
  \item{y}{A numeric vector containing the response.}
  \item{plotIt}{Logical: if TRUE, a graph is drawn.}
  \item{type}{"QQ" or "hist"}
  \item{includeIntercept}{Logical: if TRUE, the intercept effect
is plotted; otherwise, it is omitted from the plot.}
}
\value{A QQ-plot or a histogram and rugplot, or
a list if plotIt=FALSE}
\author{
W. John Braun
}
\examples{
# Jojoba oil data set
X <- p4.18[,-4]
y <- p4.18[,4]
tplot(X, y, type="hist", includeIntercept=FALSE)
title("Tests for Individual Coefficients in the Jojoba Oil Regression")
# Simulated data set where none of the predictors are in the true model:
set.seed(4571)
Z <- matrix(rnorm(400), ncol=10)
A <- matrix(rnorm(81), ncol=9)
simdata <- data.frame(Z[,1], crossprod(t(Z[,-1]),A))
names(simdata) <- c("y", paste("x", 1:9, sep=""))
X <- simdata[,-1]
y <- simdata[,1]
tplot(X, y, type="hist", includeIntercept=FALSE)
title("Tests for Individual Coefficients for the Simulated Data Set")
# NFL Data set:
X <- table.b1[,-1]
y <- table.b1[,1]
tplot(X, y, type="hist", includeIntercept=FALSE)
title("Tests for Individual Coefficients for the NFL Data Set")
# Simulated Data set where x8 is the only predictor in the true model:
X <- pathoeg[,-10]
y <- pathoeg[,10]
par(mfrow=c(2,2))
tplot(X, y)
tplot(X, y, type="QQ")
}
\keyword{graphics}
