\name{p5.12}
\alias{p5.12}
\title{Data Set for Problem 5-12}
\usage{data(p5.12)}
\description{
The \code{p5.12} data frame has 27 observations on 9
variables.  
}
\format{
  This data frame contains the following columns:
  \describe{
\item{i}{a numeric vector}
    \item{xi}{a numeric vector}
    \item{x2}{a numeric vector}
    \item{x3}{a numeric vector}
    \item{yi1}{response 1}
    \item{yi2}{response 2}
    \item{yi3}{response 3}
    \item{in211.1.gif}{a numeric vector}
    \item{si}{a numeric vector}
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
