\name{p5.10}
\alias{p5.10}
\title{Data Set for Problem 5-10}
\usage{data(p5.10)}
\description{
The \code{p5.10} data frame has 27 observations on the
effect of three factors on a printing machine's ability
to apply coloring inks on package labels.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{x1}{speed}
    \item{x2}{pressure}
    \item{x3}{distance}
    \item{yi1}{response 1}
    \item{yi2}{response 2}
    \item{yi3}{response 3}
    \item{ybar.i}{average response}
    \item{si}{standard deviation of the 3 responses}
  }
}
\source{
 Montgomery, D.C., Peck, E.A., and Vining, C.G. (2001)
Introduction to Linear Regression Analysis. 3rd Edition, John Wiley and Sons.
}
\examples{
data(p5.10)
attach(p5.10)
y.lm <- lm(ybar.i ~ x1 + x2 + x3)
plot(y.lm, which=1)
detach(p5.10)
}
\keyword{datasets}
