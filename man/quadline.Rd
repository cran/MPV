\name{quadline}
\alias{quadline}
\title{Quadratic Overlay}
\description{
Overlays a quadratic curve to a fitted quadratic model.
}
\usage{
quadline(lm.obj, ...)
}
\arguments{
  \item{lm.obj}{A \code{lm} object (a quadratic fit) }
  \item{...}{Other arguments to the \code{lines} function; e.g. \code{col}}
  }
\value{
The function superimposes a quadratic curve onto an existing scatterplot.
}
\author{W.J. Braun}
\seealso{\code{lm}}

\examples{
data(p4.18)
attach(p4.18)
y.lm <- lm(y ~ x1 + I(x1^2))
plot(x1, y)
quadline(y.lm)
detach(p4.18)
}
\keyword{models}
