\name{vif}
\alias{vif}
\title{Variance Inflation Factors}
\description{
Variance inflation factors are computed for the standard
errors of linear model coefficient estimates.
}
\usage{
vif(obj, digits=5)
}
\arguments{
  \item{obj}{ A \code{lm} object }
  \item{digits}{ Number of digits}
}
\value{
A vector of variance inflation factors corresponding to 
the coefficient estimates given in the \code{lm} object.
}
\author{J.H. Maindonald}
\seealso{\code{lm}}

\examples{
data(p4.20)
y.lm <- lm(y ~ acidtemp+acidconc+watertemp+sulfconc+amtbl, data=p4.20)
vif(y.lm)
}
\keyword{models}
