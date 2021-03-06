\name{GFplot}
\alias{GFplot}
\title{
Graphical F Plot for Significance in Regression
}
\description{
This function analyzes regression data graphically.  It 
allows visualization of the usual F-test for significance
of regression.
}
\usage{
GFplot(X, y, plotIt=TRUE, sortTrt=FALSE, type="hist", includeIntercept=TRUE, labels=FALSE)
}
\arguments{
  \item{X}{The design matrix.
}
  \item{y}{A numeric vector containing the response.}
  \item{plotIt}{Logical: if TRUE, a graph is drawn.}
  \item{sortTrt}{Logical: if TRUE, an attempt is made
at sorting the predictor effects in descending order.}
  \item{type}{"QQ" or "hist"}
  \item{includeIntercept}{Logical: if TRUE, the intercept effect
is plotted; otherwise, it is omitted from the plot.}
  \item{labels}{logical: if TRUE, names of predictor variables
are used as labels; otherwise, the design matrix column numbers
are used as labels}
}
\value{A QQ-plot or a histogram and rugplot, or
a list if plotIt=FALSE}
\author{
W. John Braun
}
\source{Braun, W.J. 2013. Regression Analysis and the QR Decomposition. Preprint.}
\examples{
# Example 1
X <- p4.18[,-4]
y <- p4.18[,4]
GFplot(X, y, type="hist", includeIntercept=FALSE)
title("Evidence of Regression in the Jojoba Oil Data")
# Example 2
set.seed(4571)
Z <- matrix(rnorm(400), ncol=10)
A <- matrix(rnorm(81), ncol=9)
simdata <- data.frame(Z[,1], crossprod(t(Z[,-1]),A))
names(simdata) <- c("y", paste("x", 1:9, sep=""))
GFplot(simdata[,-1], simdata[,1], type="hist", includeIntercept=FALSE)
title("Evidence of Regression in Simulated Data Set")
# Example 3
GFplot(table.b1[,-1], table.b1[,1], type="hist", includeIntercept=FALSE)
title("Evidence of Regression in NFL Data Set")
# An example where stepwise AIC selects the complement
# of the set of variables that are actually in the true model:
X <- pathoeg[,-10]
y <- pathoeg[,10]
par(mfrow=c(2,2))
GFplot(X, y)
GFplot(X, y, sortTrt=TRUE)
GFplot(X, y, type="QQ")
GFplot(X, y, sortTrt=TRUE, type="QQ")
X <- table.b1[,-1]  # NFL data
y <- table.b1[,1]
GFplot(X, y)

}
\keyword{graphics}
