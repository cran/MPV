\name{Juliet}
\alias{Juliet}
\title{Juliet}
\description{
\code{Juliet} has 28 rows and 9 columns. The data is of the input and output of the Spirit Still "Juliet" from Endless Summer Distillery. It is suggested to split the data by the Batch factor for ease of use.
}
\usage{Juliet}
\format{
  The data frame contains the following 9 columns. 
  \describe{
    \item{\code{Batch}}{a Factor determing how many times the volume has been through the still.}
    \item{\code{Vol1}}{Volume in litres, initial}
    \item{\code{P1}}{Percent alcohol present, initial}
    \item{\code{LAA1}}{Litres Absolute Alcohol initial, \code{Vol1*P1}}
    \item{\code{Vol2}}{Volume in litres, final}
    \item{\code{P2}}{Percent alcohol present, final}
    \item{\code{LAA2}}{Litres Absolute Alcohol final, \code{Vol2*P2}}
    \item{\code{Yield}}{Percent yield obtained, \code{LAA2/LAA1}}
    \item{\code{Date}}{Character, Date of run}
  }
}
\details{
The purpose of this information is to determine the optimal initial volume and percentage. The information is broken down by \code{Batch}. A batch factor 1 means that it 
is the first time the liquid has gone through the spirit still. The first run through the still should have the most loss due to the "heads" and "tails".
Literature states that the first run through a spirit still should yield 70 percent. 
A batch factor 2 means that it is the second time the liquid has gone through the spirit still. 
A batch factor 3 means that it is the third time or more that the liquid has gone through the spirit still. 
Each subsequent distillation should result in a higher yield, never to exceed 95 percent. 
}
\source{
Charisse Woods, Endless Summer Distillery, (2015).
}
\examples{
summary(Juliet)

#Split apart the Batch factor for easier use.
juliet<-split(Juliet,Juliet$Batch)
juliet1<-juliet$'1'
juliet2<-juliet$'2'
juliet3<-juliet$'3'

plot(LAA1~LAA2,data=Juliet)
plot(LAA1~LAA2,data=juliet1)

}
\keyword{datasets}
