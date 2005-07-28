\name{PRESS}
\alias{PRESS}
\title{PRESS statistic}
\description{
Computation of Allen's PRESS statistic for an lm object.
}
\usage{
PRESS(x)
}
\arguments{
  \item{x}{An \code{lm} object }
  }
\value{
Allen's PRESS statistic.
}
\author{W.J. Braun}
\seealso{\code{lm}}

\examples{
data(p4.18)
attach(p4.18)
y.lm <- lm(y ~ x1 + I(x1^2))
PRESS(y.lm)
detach(p4.18)
}
\keyword{models}
