\name{lesions}
\alias{lesions}
\title{Lesions in Rat Colons}
\usage{lesions}
\description{
Numbers of aberrant crypt foci (ACF) in each of 
six cross-sectional regions of the 
 colons of 66 rats subjected to varying 
doses of the carcinogen azoxymethane (AOM), sacrificed
at 3 different times. 
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{T}{Incubation time factor, levels: 6, 12 and 18 weeks}
    \item{INJ}{Number of injections}
    \item{SECT}{Section of colon, a factor with levels 1 through 6, where 1
denotes the proximal end of the colon and 6 denotes the distal end}  
    \item{RAT}{Label for animal within a particular T-INJ factor level 
    combination}
    \item{ACF.Total}{Total number of ACF lesions in a section of a rat's colon}
    \item{ACF.total.mult}{Sum of ACF multiplicities for a section of a rat's colon}
    \item{id}{Identifier for each of the 66 rats.}
}
}
\source{
 Ranjana P. Bird, University of Northern British Columbia,
Prince George, Canada.
}
\references{
E.A. McLellan, A. Medline and R.P. Bird.  Dose response and
proliferative characteristics of aberrant crypt foci: putative
preneoplastic lesions in rat colon.  Carcinogenesis, 12(11): 2093-2098, 
1991.
}
\examples{
summary(lesions)
ACF.All <- aggregate(ACF.Total ~  id + INJ + T, FUN=sum, data = lesions)
lesions.glm <- glm(ACF.Total ~ INJ * T, data = ACF.All, family=poisson)
summary(lesions.glm)
lesions.qp <- glm(ACF.Total ~ INJ * T, data = ACF.All, family=quasipoisson)
summary(lesions.qp)
lesions.noInt <- glm(ACF.Total ~ INJ + T, data = ACF.All, family=quasipoisson)
summary(lesions.noInt)
}
\keyword{datasets}
