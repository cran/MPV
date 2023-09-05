\name{p5.12}
\alias{p5.12}
\title{Data Set for Problem 5-12}
\usage{data(p5.12)}
\description{
The \code{p5.12} data frame has 27 observations on 3 variables, with
responses replicated 3 times.  Averages and standard deviations are calculated
for each level of the experimental design.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{i}{numeric, experimental run number}
    \item{xi}{numeric}
    \item{x2}{numeric}
    \item{x3}{numeric}
    \item{yi1}{response 1}
    \item{yi2}{response 2}
    \item{yi3}{response 3}
    \item{ybari}{average of 3 responses at ith level}
    \item{si}{standard deviation of 3 responses at ith level}
  }
}
\source{
 Montgomery, D.C., Peck, E.A., and Vining, C.G. (2021)
Introduction to Linear Regression Analysis. 6th Edition, John Wiley and Sons.
}
\references{
Vining, G. and Myers, R. (1990) "Combining Taguchi and Response Surface Philosophies: A Dual Response Approach,"
Journal of Quality Technology, 22, 15-22. 
}
\examples{
y.lm <- lm(ybari ~ xi + x2 + x3, data = p5.12)
plot(y.lm, which=1)
}
\keyword{datasets}
