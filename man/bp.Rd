\name{bp}
\alias{bp}
\title{Blood Pressure Measurements on a Single Adult Male}
\description{
Systolic and diastolic blood pressure measurement readings
were taken on a 56-year-old male over a 39 day period, sometimes
in the mornings (AM) and sometimes in the evening (PM). 
Varying number of replicate
measurements were taken at each time point.  
}
\usage{bp}
\format{
  A data frame with 121 observations on the following 4 variables.
  \describe{
\item{\code{TimeofDay}}{factor with levels AM and PM}
\item{\code{Date}}{numeric}
\item{\code{Systolic}}{numeric}
\item{\code{Diastolic}}{numeric}
  }
}
\examples{
require(lattice)
xyplot(Date ~ Diastolic|TimeofDay, groups=cut(Systolic, c(0, 130, 140,
   200)), data = bp, col=c(3, 1, 2), pch=16)
matplot(bp[, c(3, 4)], type="l", lwd=2, ylab="Pressure")
n <- nrow(bp)
abline(v=(1:n)[bp[,1]=="PM"]-.5, col="grey")
abline(v=(1:n)[bp[,1]=="PM"], col="grey")
abline(v=(1:n)[bp[,1]=="PM"]+.5, col="grey")
bp.stk <- stack(bp, c("Systolic", "Diastolic"))
bp.tmp <- rbind(bp[,1:2], bp[,1:2])
bp.stk <- cbind(bp.tmp, bp.stk)
names(bp.stk) <- c("TimeofDay", "Date", "Pressure", "Type")
reps <- NULL
for (j in rle(paste(bp.stk$Date, bp.stk$TimeofDay))$lengths) reps <- c(reps, (1:j))
bp.stk$Rep <- reps
xyplot(Pressure ~ I(Date+Rep/24)|TimeofDay, groups=Type, data = bp.stk, xlab="Date", pch=16)
}
\keyword{datasets}
