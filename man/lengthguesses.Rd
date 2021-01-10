\name{lengthguesses}
\alias{lengthguesses}
\title{Length Guesses Data}
\usage{data(lengthguesses)}
\description{
The \code{lengthguesses} list consists of 2 numeric vectors, one
giving the metric-converted length guesses (in feet) of an auditorium
whose actual length (in meters) was 13.1m, and the other containing
the length guesses of 69 others (in meters).
}
\format{
  This list contains the following columns:
  \describe{
    \item{imperial}{a numeric vector of 69 student guesses as to the
length of an auditorium using the imperial system, converted to meters.}
    \item{metric}{a numeric vector of 44 student guesses as to the length of an
auditorium using the metric system.}
  }
}
\source{
Hills, M. and the M345 Course Team (1986) M345 Statistical Methods, Unit 1:
Data, distributions and uncertainty, Milton Keynes: The Open University.
Tables 2.1 and 2.4. 
}
\references{
Hand, D.J., Daly, F., Lunn, A.D., McConway, K.J. and Ostrowski, E. (1994)
A Handbook of Small Data Sets. Boca Raton: Chapman & Hall/CRC.  
}
\examples{
with(lengthguesses, t.test(imperial, metric))
}
\keyword{datasets}
