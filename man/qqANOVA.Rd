\name{qqANOVA}
\alias{qqANOVA}
\title{
QQ Plot for Analysis of Variance
}
\description{ This function is used to display the weight of the 
evidence against null main effects in data coming from a 1 factor 
design, using a QQ plot.  In practice this method is often called via 
the function GANOVA.
}
\usage{
qqANOVA(x, y, plot.it = TRUE, xlab = deparse(substitute(x)),
    ylab = deparse(substitute(y)), ...)
}
\arguments{
  \item{x}{numeric vector of errors}
  \item{y}{numeric vector of scaled responses}
  \item{plot.it}{logical vector indicating whether to plot
or not}
  \item{xlab}{character, x-axis label}
  \item{ylab}{character, y-axis label}
  \item{...}{any other arguments for the plot function}
}
\value{A QQ plot is drawn.}
\author{
W. John Braun
}
\keyword{graphics}
