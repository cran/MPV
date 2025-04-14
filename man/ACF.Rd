\name{ACF}
\alias{ACF}
\title{Aberrant Crypt Foci in Rat Colons}
\usage{ACF}
\description{
Numbers of aberrant crypt foci (ACF) in 
colons of 66 rats subjected to a various numbers of 
dose of the carcinogen azoxymethane (AOM), sacrificed
at 3 different times. 
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{INJ}{The number of carcinogen injections}
    \item{T}{Time of sacrifice, in weeks following injection
of AOM}
    \item{COUNT}{The number of ACF observed in each rat colon}
  }
}
\source{
 Ranjana P. Bird, Faculty of Human Ecology, University of Manitoba,
Winnipeg, Canada.
}
\references{
E.A. McLellan, A. Medline and R.P. Bird.  Dose response and
proliferative characteristics of aberrant crypt foci: putative
preneoplastic lesions in rat colon.  Carcinogenesis, 12(11): 2093-2098, 
1991.
}
\examples{
sapply(split(ACF$COUNT,ACF$T),var)
}
\keyword{datasets}
