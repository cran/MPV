\name{lof.lm}
\alias{lof.lm}
\title{Lack of Fit Test}
\description{
Tests lack of fit for simple linear regression.
}
\usage{
lof.lm(lm.obj, approx = FALSE, call.plot=TRUE)
}
\arguments{
  \item{lm.obj}{A \code{lm} object }
  \item{approx}{binary, if set to TRUE, it gives an approximation
when no replicates are available  }
  \item{call.plot}{binary}
}
\value{
F-test for lack of fit.
}
\author{W.J. Braun}
\seealso{\code{lm}}

\examples{
data(p4.18)
y.lm <- lm(y ~ x1, data=p4.18)
lof.lm(y.lm)
}
\keyword{models}
