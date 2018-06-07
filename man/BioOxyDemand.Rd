\name{BioOxyDemand}
\alias{BioOxyDemand}
\title{Biochemical Oxygen Demand}
\usage{data(BioOxyDemand)}
\description{
The \code{BioOxyDemand} data frame has 14 rows and 2 columns.
}
\format{
This data frame contains the following columns:
\describe{
\item{x}{
a numeric vector
}
\item{y}{
a numeric vector
}
}
}
\source{
  Devore, J. L. (2000) \emph{Probability and Statistics for Engineering and the Sciences (5th ed)}, Duxbury
}
\examples{
plot(BioOxyDemand)
summary(lm(y ~ x, data = BioOxyDemand))
}
\keyword{datasets}

