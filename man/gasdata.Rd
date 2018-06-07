\name{gasdata}
\alias{gasdata}
\docType{data}
\title{Natural Gas Consumption in a Single-Family Residence}
\description{
This data frame contains the average monthly volume
of natural gas used in the furnace of a 1600 square foot
house located in London, Ontario, for each month from 2006 until 2011. 
It also contains the average temperature 
for each month, and a measure of degree days.  Insulation was 
added to the roof on one occasions, the walls were insulated on a second
occasion, 
and the mid-efficiency furnace was replaced with a high-efficiency
furnace on a third occasion. 

}
\usage{data("gasdata")}
\format{
  A data frame with 70 observations on the following 9 variables.
  \describe{
\item{\code{month}}{numeric 1=January, 12=December}
\item{\code{degreedays}}{numeric, Celsius}
\item{\code{cubicmetres}}{total volume of gas used in a month}
\item{\code{dailyusage}}{average amount of gas used per day}
\item{\code{temp}}{average temperature in Celsius}
\item{\code{year}}{numeric}
\item{\code{I1}}{indicator that  roof insulation is present}
\item{\code{I2}}{indicator that wasll insulation is present}
\item{\code{I3}}{indicator that high efficiency furnace is present}
  }
}
\keyword{datasets}
