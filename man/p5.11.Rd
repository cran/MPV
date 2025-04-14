\name{p5.11}
\alias{p5.11}
\title{Data Set for Problem 5-11}
\usage{data(p5.11)}
\description{
The \code{p5.11} data frame has 8 observations on an 
experiment with a catapult.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{x1}{hook}
    \item{x2}{arm length}
    \item{x3}{start angle}
    \item{x4}{stop angle}
    \item{yi1}{response 1}
    \item{yi2}{response 2}
    \item{yi3}{response 3}
  }
}
\source{
 Montgomery, D.C., Peck, E.A., and Vining, C.G. (2001)
Introduction to Linear Regression Analysis. 3rd Edition, John Wiley and Sons.
}
\examples{
data(p5.11)
attach(p5.11)
ybar.i <- apply(p5.11[,5:7], 1, mean)
sd.i <- apply(p5.11[,5:7], 1, sd)
y.lm <- lm(ybar.i ~ x1 + x2 + x3 + x4)
plot(y.lm, which=1)
detach(p5.11)
}
\keyword{datasets}
