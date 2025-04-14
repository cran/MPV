\name{soilstudy}
\alias{soilstudy}
\docType{data}
\title{Soil Moisture Data}
\description{
Percent soil moisture measurements at 26 different locations in a forest in
southwestern British Columbia.  Some of the locations were in stands that had
been thinned.   
}

\usage{data("soilstudy")}
\format{
  A data frame with 26 observations on the following 3 variables.
  \describe{
    \item{\code{location}}{character vector identifying forest stand}
    \item{\code{moisture}}{numeric vector, percentage moisture content}
    \item{\code{treatment}}{character vector identifying fuel treatment: thinned or unthinned}
  }
}
\source{Millikin, R.L., Braun, W.J., Alexander, M.E., Fani, S. (2024), 
The Impact of Fuel Thinning on the Microclimate in Coastal Rainforest Stands of Southwestern British Columbia, Canada.
Fire. Vol 7(8), 2024, pp 285-309. 
}
\keyword{datasets}
