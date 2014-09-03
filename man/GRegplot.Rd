\name{GRegplot}
\alias{GRegplot}
\title{
Graphical Regression Plot 
}
\description{
This function analyzes regression data graphically.  It 
allows visualization of the usual F-test for significance
of regression.
}
\usage{
GRegplot(X, y, sortTrt=FALSE, includeIntercept=TRUE, type=c("hist", "dot"))
}
\arguments{
  \item{X}{The design matrix.
}
  \item{y}{A numeric vector containing the response.}
  \item{sortTrt}{Logical: if TRUE, an attempt is made
at sorting the predictor effects in descending order.}
  \item{includeIntercept}{Logical: if TRUE, the intercept effect
is plotted; otherwise, it is omitted from the plot.}
  \item{type}{Character: hist, for histogram; dot, for stripchart}
}
\value{A histogram or dotplot and rugplot}
\author{
W. John Braun
}
\source{Braun, W.J. 2014. Visualization of Evidence in 
Regression Analysis with the QR Decomposition. Preprint.}
\examples{
# Example 1
X <- p4.18[,-4]
y <- p4.18[,4]
GRegplot(X, y, includeIntercept=FALSE)
title("Evidence of Regression in the Jojoba Oil Data")
# Example 2
set.seed(4571)
Z <- matrix(rnorm(400), ncol=10)
A <- matrix(rnorm(81), ncol=9)
simdata <- data.frame(Z[,1], crossprod(t(Z[,-1]),A))
names(simdata) <- c("y", paste("x", 1:9, sep=""))
GRegplot(simdata[,-1], simdata[,1], includeIntercept=FALSE)
title("Evidence of Regression in Simulated Data Set")
# Example 3
GRegplot(table.b1[,-1], table.b1[,1], includeIntercept=FALSE)
title("Evidence of Regression in NFL Data Set")
# An example where stepwise AIC selects the complement
# of the set of variables that are actually in the true model:
X <- pathoeg[,-10]
y <- pathoeg[,10]
par(mfrow=c(2,1))
GRegplot(X, y)
GRegplot(X, y, sortTrt=TRUE)
X <- table.b1[,-1]  # NFL data
y <- table.b1[,1]
GRegplot(X, y)
}
\keyword{graphics}
